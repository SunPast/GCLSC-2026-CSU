"""BERT-Set-VAE for antibody sequences (v5 final).

Architecture (literature-grounded, minimal, robust):
  Encoder: BERT bidirectional (6 layers) → Set Transformer pooling (16 learnable
    queries × cross-attention) → mu/logvar → z (D-dim)
  Decoder: ResNet (8 Conv1D blocks), FiLM-conditioned on z, per-position
    independent decoding.  No self-attention, no teacher forcing, no
    autoregression.  Cannot ignore z — must use it to reconstruct.

Key design choices:
  - Set Transformer pooling (Lee et al. 2019): 16 learnable queries extract
    diverse positional information, avoiding FR-dominated mean pooling.
  - ResNet FiLM decoder: z injects via FiLM (γ, β) in every ResBlock.
    Per-position decoding — no position-to-position interaction.
    Posterior collapse prevented architecturally.
  - No disentanglement (TC is not estimated).  Attribute conditioning is
    deferred to Stage 2 latent diffusion.
  - β-VAE loss: recon_CE + β·KL(q(z|x)||N(0,I)), fixed β=0.1.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict, Optional, Tuple

import torch
import torch.nn as nn
import torch.nn.functional as F


# =============================================================================
# Config
# =============================================================================

@dataclass
class BERTSetVAEConfig:
    vocab_size: int = 26
    max_len: int = 140

    # Encoder
    d_model: int = 512
    n_enc_layers: int = 6
    n_heads: int = 8
    d_ff: int = 2048
    dropout: float = 0.1

    # Set Transformer pooling
    n_queries: int = 16

    # Latent
    latent_dim: int = 512  # D
    beta: float = 0.1       # KL weight

    # ResNet decoder
    n_dec_blocks: int = 8
    dec_hidden: int = 256

    # Special tokens
    pad_id: int = 0
    eos_id: int = 3
    h_chain_id: int = 4
    l_chain_id: int = 5


# =============================================================================
# Set Transformer Pooling (Lee et al. 2019)
# =============================================================================

class SetTransformerPooling(nn.Module):
    """Learnable query vectors attend to encoder outputs via cross-attention.

    N learnable queries → cross-attn on encoder outputs → N pooled vectors.
    Flatten → linear → output_dim.  More informative than mean pooling.
    """

    def __init__(self, d_model: int, n_queries: int, output_dim: int,
                 n_heads: int = 8):
        super().__init__()
        self.queries = nn.Parameter(torch.randn(n_queries, d_model) * 0.02)
        self.cross_attn = nn.MultiheadAttention(
            d_model, n_heads, batch_first=True,
        )
        self.norm = nn.LayerNorm(d_model, eps=1e-6)
        self.proj = nn.Linear(n_queries * d_model, output_dim)

    def forward(self, enc_out: torch.Tensor,
                mask: torch.Tensor) -> torch.Tensor:
        """enc_out: (B, T, E); mask: (B, T) bool. → (B, output_dim)."""
        B = enc_out.size(0)
        q = self.queries.unsqueeze(0).expand(B, -1, -1)  # (B, N, E)
        # Cross-attention: queries attend to encoder tokens
        attended, _ = self.cross_attn(
            q, enc_out, enc_out,
            key_padding_mask=~mask,
        )
        attended = self.norm(attended + q)  # residual
        return self.proj(attended.reshape(B, -1))


# =============================================================================
# FiLM (Feature-wise Linear Modulation) layer
# =============================================================================

class FiLM(nn.Module):
    """z → (γ, β) for modulating feature maps.  (Perez et al. 2018, AAAI)."""

    def __init__(self, z_dim: int, n_features: int):
        super().__init__()
        self.net = nn.Sequential(
            nn.Linear(z_dim, n_features),
            nn.GELU(),
            nn.Linear(n_features, n_features * 2),
        )

    def forward(self, z: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        """z: (B, z_dim); x: (B, T, C) → (B, T, C)."""
        params = self.net(z)  # (B, 2C)
        gamma, beta = params.chunk(2, dim=-1)
        return x * (1.0 + gamma.unsqueeze(1)) + beta.unsqueeze(1)


# =============================================================================
# ResNet FiLM Decoder Block
# =============================================================================

class ResFiLMBlock(nn.Module):
    """Conv1D residual block with FiLM conditioning from z."""

    def __init__(self, hidden: int, z_dim: int):
        super().__init__()
        self.norm1 = nn.LayerNorm(hidden, eps=1e-6)
        self.conv1 = nn.Conv1d(hidden, hidden * 2, kernel_size=1)
        self.norm2 = nn.LayerNorm(hidden * 2, eps=1e-6)
        self.conv2 = nn.Conv1d(hidden * 2, hidden, kernel_size=1)
        self.film = FiLM(z_dim, hidden)
        self.film2 = FiLM(z_dim, hidden * 2)

    def forward(self, x: torch.Tensor, z: torch.Tensor) -> torch.Tensor:
        # x: (B, T, C), z: (B, z_dim)
        residual = x
        x = self.norm1(x)
        x = self.film(z, x)  # FiLM modulation
        x = x.permute(0, 2, 1)  # (B, T, C) → (B, C, T)
        x = self.conv1(x).permute(0, 2, 1)  # → (B, T, 2C)
        x = F.gelu(x)
        x = self.norm2(x)
        x = self.film2(z, x)
        x = x.permute(0, 2, 1)
        x = self.conv2(x).permute(0, 2, 1)  # → (B, T, C)
        return x + residual


# =============================================================================
# BERT-Set-VAE Model
# =============================================================================

class BERTSetVAE(nn.Module):
    """VAE with BERT encoder, Set Transformer pooling, and ResNet-FiLM decoder.

    Encoder: BERT → Set Transformer (16 queries) → mu/logvar → z
    Decoder: ResNet (8 FiLM blocks), per-position independent, no self-attn.
    """

    def __init__(self, cfg: BERTSetVAEConfig):
        super().__init__()
        self.cfg = cfg

        # ---- Embeddings ----
        self.embed = nn.Embedding(cfg.vocab_size, cfg.d_model, padding_idx=cfg.pad_id)
        self.pos = nn.Parameter(torch.randn(cfg.max_len, cfg.d_model) * 0.02)
        self.drop = nn.Dropout(cfg.dropout)

        # ---- Encoder ----
        enc_layer = nn.TransformerEncoderLayer(
            d_model=cfg.d_model, nhead=cfg.n_heads,
            dim_feedforward=cfg.d_ff, dropout=cfg.dropout,
            activation="gelu", batch_first=True, norm_first=True,
        )
        self.encoder = nn.TransformerEncoder(enc_layer, num_layers=cfg.n_enc_layers)

        # ---- Set Transformer Pooling ----
        self.pool = SetTransformerPooling(
            cfg.d_model, cfg.n_queries, cfg.latent_dim * 2, cfg.n_heads,
        )
        # pool output → split into mu, logvar
        self.fc_mu = nn.Linear(cfg.latent_dim * 2, cfg.latent_dim)
        self.fc_logvar = nn.Linear(cfg.latent_dim * 2, cfg.latent_dim)

        # ---- Decoder (ResNet FiLM) ----
        # Initial projection: flat → per-position tokens
        self.dec_proj = nn.Linear(cfg.latent_dim, cfg.max_len * cfg.dec_hidden)
        self.dec_blocks = nn.ModuleList([
            ResFiLMBlock(cfg.dec_hidden, cfg.latent_dim)
            for _ in range(cfg.n_dec_blocks)
        ])
        self.dec_norm = nn.LayerNorm(cfg.dec_hidden, eps=1e-6)
        self.out_head = nn.Linear(cfg.dec_hidden, cfg.vocab_size)

        self._init_weights()

    def _init_weights(self):
        for m in self.modules():
            if isinstance(m, nn.Linear):
                nn.init.xavier_uniform_(m.weight)
                if m.bias is not None:
                    nn.init.constant_(m.bias, 0.0)

    # ------------------------------------------------------------------
    # Encoder
    # ------------------------------------------------------------------
    def encode(self, tokens: torch.Tensor, mask: torch.Tensor
               ) -> Tuple[torch.Tensor, torch.Tensor]:
        B, T = tokens.shape
        x = self.embed(tokens)
        x = x + self.pos[:T].unsqueeze(0)
        x = self.drop(x)
        x = self.encoder(x, src_key_padding_mask=~mask)
        pooled = self.pool(x, mask)
        mu = self.fc_mu(pooled)
        logvar = self.fc_logvar(pooled)
        # Clamp for FP16 safety: prevents NaN from exp(logvar) overflow
        # + mu clamp prevents encoder output spikes that trigger collapse
        logvar = logvar.clamp(min=-10, max=10)
        mu = mu.clamp(min=-10, max=10)
        return mu, logvar

    def reparameterize(self, mu: torch.Tensor, logvar: torch.Tensor,
                       sample: bool = True) -> torch.Tensor:
        if not sample:
            return mu
        std = torch.exp(0.5 * logvar.clamp(min=-20, max=20))
        return mu + torch.randn_like(std) * std

    # ------------------------------------------------------------------
    # Decoder
    # ------------------------------------------------------------------
    def decode(self, z: torch.Tensor) -> torch.Tensor:
        """z (B, D) → logits (B, T, V). Per-position independent."""
        B, D = z.shape
        T = self.cfg.max_len
        # Project z to per-position tokens
        x = self.dec_proj(z)  # (B, T*C)
        x = x.view(B, T, self.cfg.dec_hidden)  # (B, T, C)
        # ResNet FiLM blocks
        for blk in self.dec_blocks:
            x = blk(x, z)  # z injected at every block via FiLM
        x = self.dec_norm(x)
        return self.out_head(x)  # (B, T, V)

    # ------------------------------------------------------------------
    # Generation (greedy sampler)
    # ------------------------------------------------------------------
    @torch.no_grad()
    def generate(self, z: torch.Tensor, temperature: float = 1.0) -> torch.Tensor:
        """Decode z to token ids. temperature=1 uses argmax; else multinomial."""
        logits = self.decode(z)
        if temperature <= 0:
            raise ValueError("temperature must be > 0")
        if temperature == 1.0:
            return logits.argmax(dim=-1)
        probs = torch.softmax(logits / temperature, dim=-1)
        return torch.multinomial(probs.view(-1, probs.size(-1)), 1).view(probs.shape[:2])

    # ------------------------------------------------------------------
    # Forward
    # ------------------------------------------------------------------
    def forward(self, tokens: torch.Tensor, mask: torch.Tensor,
                beta: Optional[float] = None,
                sample: bool = True) -> Dict[str, torch.Tensor]:
        b = beta if beta is not None else self.cfg.beta

        mu, logvar = self.encode(tokens, mask)
        z = self.reparameterize(mu, logvar, sample=sample)
        logits = self.decode(z)

        # ---- Losses ----
        # Reconstruction
        recon = F.cross_entropy(
            logits.reshape(-1, self.cfg.vocab_size),
            tokens.reshape(-1),
            reduction="none",
        )
        recon = (recon * mask.reshape(-1).float()).sum() / (mask.sum() + 1e-8)

        # KL divergence
        kl = -0.5 * (1 + logvar - mu.pow(2) - logvar.exp()).sum(dim=-1).mean()

        loss = recon + b * kl
        return {
            "loss": loss, "recon": recon, "kl": kl,
            "mu": mu, "logvar": logvar, "z": z, "logits": logits,
        }
