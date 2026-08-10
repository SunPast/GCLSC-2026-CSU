"""BERT-Set-VAE v2 for antibody sequences (set-latent).

Encoder: BERT → multi-query pooling → per-token mu/logvar → z_tokens (B, Q, D).
Decoder: position cross-attention over z_tokens, FiLM, Conv1d blocks.
Loss: recon_CE + β·KL.

`kernel_size` is configurable (production freeze uses 1 to limit
decoder-side local shortcuts that can induce posterior collapse).
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
class BERTSetVAEConfigV2:
    vocab_size: int = 26
    max_len: int = 140

    # Encoder
    d_model: int = 512
    n_enc_layers: int = 6
    n_heads: int = 8
    d_ff: int = 2048
    dropout: float = 0.1

    # Set Transformer queries
    n_queries: int = 64

    # VAE latent
    latent_dim: int = 512  # per-token latent dim
    beta: float = 0.1

    # Decoder
    n_dec_blocks: int = 8
    dec_hidden: int = 256
    kernel_size: int = 1  # Conv1d kernel (1 = pointwise; prefer for anti-collapse)

    # Special tokens
    pad_id: int = 0
    eos_id: int = 3
    h_chain_id: int = 4
    l_chain_id: int = 5


# =============================================================================
# Set Transformer Multi-Query Pooling
# =============================================================================

class MultiQueryPooling(nn.Module):
    """64 learnable query vectors independently attend to encoder tokens.

    Returns (B, N, d_model) — all tokens retained for decoder cross-attention.
    """

    def __init__(self, d_model: int, n_queries: int = 64, n_heads: int = 8,
                 n_groups: int = 8):
        super().__init__()
        self.n_queries = n_queries
        self.n_groups = n_groups
        queries_per_group = n_queries // n_groups

        # Create group-wise queries: each group learns to focus on different
        # aspects of the sequence (FR framework, CDR loops, charge patterns, etc.)
        self.group_queries = nn.ParameterList([
            nn.Parameter(torch.randn(queries_per_group, d_model) * 0.02)
            for _ in range(n_groups)
        ])

        self.cross_attn = nn.MultiheadAttention(
            d_model, n_heads, batch_first=True,
        )
        self.group_norm = nn.LayerNorm(d_model, eps=1e-6)

    def forward(self, enc_out: torch.Tensor,
                mask: torch.Tensor) -> torch.Tensor:
        """enc_out: (B, T, E); mask: (B, T) → (B, 64, E)."""
        B = enc_out.size(0)
        # Concatenate all group queries
        queries = torch.cat([q.unsqueeze(0).expand(B, -1, -1)
                             for q in self.group_queries], dim=1)
        attended, _ = self.cross_attn(
            queries, enc_out, enc_out,
            key_padding_mask=~mask,
        )
        return self.group_norm(attended + queries)


# =============================================================================
# FiLM (global z → γ, β for all positions)
# =============================================================================

class GlobalFiLM(nn.Module):
    """Pool z_tokens → global z → γ, β for all positions."""

    def __init__(self, z_dim: int, z_tokens: int, n_features: int):
        super().__init__()
        input_dim = z_tokens * z_dim
        self.net = nn.Sequential(
            nn.Linear(input_dim, n_features),
            nn.GELU(),
            nn.Linear(n_features, n_features * 2),
        )

    def forward(self, z_tokens: torch.Tensor, x: torch.Tensor) -> torch.Tensor:
        """z_tokens: (B, N, D); x: (B, T, C) → (B, T, C)."""
        B = z_tokens.size(0)
        z_flat = z_tokens.reshape(B, -1)
        params = self.net(z_flat)
        gamma, beta = params.chunk(2, dim=-1)
        return x * (1.0 + gamma.unsqueeze(1)) + beta.unsqueeze(1)


# =============================================================================
# Position-Aware Cross-Attention
# =============================================================================

class PositionCrossAttn(nn.Module):
    """140 learnable position queries attend to 64 z_tokens.

    Each position learns which z tokens carry information relevant
    to that sequence position.  No position bias — implicit alignment.
    """

    def __init__(self, d_model: int, max_len: int = 140, n_heads: int = 8):
        super().__init__()
        self.pos_queries = nn.Parameter(torch.randn(max_len, d_model) * 0.02)
        self.cross_attn = nn.MultiheadAttention(
            d_model, n_heads, batch_first=True,
        )
        self.norm = nn.LayerNorm(d_model, eps=1e-6)

    def forward(self, z_tokens: torch.Tensor) -> torch.Tensor:
        """z_tokens: (B, 64, E) → (B, 140, E)."""
        B = z_tokens.size(0)
        q = self.pos_queries.unsqueeze(0).expand(B, -1, -1)
        attended, _ = self.cross_attn(q, z_tokens, z_tokens)
        return self.norm(attended + q)


# =============================================================================
# Decoder Block (FiLM + Conv1d)
# =============================================================================

class DecoderBlockV2(nn.Module):
    """Conv1D residual block with global FiLM + local context."""

    def __init__(self, hidden: int, z_dim: int, z_tokens: int,
                 kernel_size: int = 3):
        super().__init__()
        self.norm1 = nn.LayerNorm(hidden, eps=1e-6)
        self.conv1 = nn.Conv1d(hidden, hidden * 2, kernel_size=kernel_size,
                                padding=kernel_size // 2)
        self.norm2 = nn.LayerNorm(hidden * 2, eps=1e-6)
        self.conv2 = nn.Conv1d(hidden * 2, hidden, kernel_size=kernel_size,
                                padding=kernel_size // 2)
        self.film = GlobalFiLM(z_dim, z_tokens, hidden)
        self.film2 = GlobalFiLM(z_dim, z_tokens, hidden * 2)

    def forward(self, x: torch.Tensor, z_tokens: torch.Tensor) -> torch.Tensor:
        residual = x
        x = self.norm1(x)
        x = self.film(z_tokens, x)
        x = x.permute(0, 2, 1)
        x = self.conv1(x).permute(0, 2, 1)
        x = F.gelu(x)
        x = self.norm2(x)
        x = self.film2(z_tokens, x)
        x = x.permute(0, 2, 1)
        x = self.conv2(x).permute(0, 2, 1)
        return x + residual


# =============================================================================
# BERT-Set-VAE v2
# =============================================================================

class BERTSetVAEV2(nn.Module):
    """VAE v2: Token-level latent + position cross-attn decoder."""

    def __init__(self, cfg: BERTSetVAEConfigV2):
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
        self.pool = MultiQueryPooling(cfg.d_model, cfg.n_queries, cfg.n_heads)

        # ---- VAE per-token reparameterization ----
        self.fc_mu = nn.Linear(cfg.d_model, cfg.latent_dim)
        self.fc_logvar = nn.Linear(cfg.d_model, cfg.latent_dim)

        # ---- Decoder ----
        self.pos_cross_attn = PositionCrossAttn(cfg.d_model, cfg.max_len, cfg.n_heads)
        self.dec_proj = nn.Linear(cfg.d_model, cfg.dec_hidden)
        self.dec_blocks = nn.ModuleList([
            DecoderBlockV2(cfg.dec_hidden, cfg.latent_dim, cfg.n_queries,
                           cfg.kernel_size)
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
        z_raw = self.pool(x, mask)  # (B, 64, 512)
        mu = self.fc_mu(z_raw)
        logvar = self.fc_logvar(z_raw)
        logvar = logvar.clamp(min=-10, max=10)
        mu = mu.clamp(min=-10, max=10)
        return mu, logvar

    def reparameterize(self, mu: torch.Tensor, logvar: torch.Tensor,
                       sample: bool = True) -> torch.Tensor:
        if not sample:
            return mu
        std = torch.exp(0.5 * logvar)
        return mu + torch.randn_like(std) * std

    # ------------------------------------------------------------------
    # Decoder
    # ------------------------------------------------------------------
    def decode(self, z_tokens: torch.Tensor) -> torch.Tensor:
        """z_tokens: (B, 64, D) → logits (B, 140, V)."""
        B, N, D = z_tokens.shape
        # Position queries attend to z tokens
        x = self.pos_cross_attn(z_tokens)  # (B, 140, 512)
        x = self.dec_proj(x)  # (B, 140, 256)
        for blk in self.dec_blocks:
            x = blk(x, z_tokens)
        x = self.dec_norm(x)
        return self.out_head(x)

    # ------------------------------------------------------------------
    # Generation
    # ------------------------------------------------------------------
    @torch.no_grad()
    def generate(self, z_tokens: torch.Tensor,
                 temperature: float = 1.0) -> torch.Tensor:
        """Decode z_tokens to token ids. temperature=1 uses argmax; else multinomial."""
        logits = self.decode(z_tokens)
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
        z_tokens = self.reparameterize(mu, logvar, sample=sample)
        logits = self.decode(z_tokens)

        # Reconstruction
        recon = F.cross_entropy(
            logits.reshape(-1, self.cfg.vocab_size),
            tokens.reshape(-1),
            reduction="none",
        )
        recon = (recon * mask.reshape(-1).float()).sum() / (mask.sum() + 1e-8)

        # KL: 64 tokens × D dims each
        kl = -0.5 * (1 + logvar - mu.pow(2) - logvar.exp()).sum(dim=(1, 2)).mean()

        loss = recon + b * kl
        return {
            "loss": loss, "recon": recon, "kl": kl,
            "mu": mu, "logvar": logvar, "z_tokens": z_tokens, "logits": logits,
        }
