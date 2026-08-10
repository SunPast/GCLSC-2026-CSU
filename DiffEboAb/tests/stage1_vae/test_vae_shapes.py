"""Shape / smoke tests for Stage-1 VAE models (tiny configs, CPU)."""
from __future__ import annotations

import sys
from pathlib import Path

import torch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.models.bert_set_vae import BERTSetVAE, BERTSetVAEConfig
from src.models.bert_set_vae_v2 import BERTSetVAEConfigV2, BERTSetVAEV2


def _tiny_vec_cfg() -> BERTSetVAEConfig:
    return BERTSetVAEConfig(
        vocab_size=26,
        max_len=32,
        d_model=64,
        n_enc_layers=1,
        n_heads=4,
        d_ff=128,
        dropout=0.0,
        n_queries=4,
        latent_dim=16,
        n_dec_blocks=2,
        dec_hidden=32,
        beta=0.1,
    )


def _tiny_set_cfg() -> BERTSetVAEConfigV2:
    # PositionCrossAttn expects z_tokens dim == d_model (production: both 512).
    return BERTSetVAEConfigV2(
        vocab_size=26,
        max_len=32,
        d_model=64,
        n_enc_layers=1,
        n_heads=4,
        d_ff=128,
        dropout=0.0,
        n_queries=8,
        latent_dim=64,
        n_dec_blocks=2,
        dec_hidden=32,
        kernel_size=1,
        beta=0.1,
    )


def test_vec_vae_forward_shapes():
    cfg = _tiny_vec_cfg()
    model = BERTSetVAE(cfg).eval()
    B, T = 2, cfg.max_len
    tokens = torch.randint(6, cfg.vocab_size, (B, T))
    mask = torch.ones(B, T, dtype=torch.bool)
    out = model(tokens, mask, beta=0.01, sample=True)
    assert out["logits"].shape == (B, T, cfg.vocab_size)
    assert out["mu"].shape == (B, cfg.latent_dim)
    assert out["z"].shape == (B, cfg.latent_dim)
    assert torch.isfinite(out["loss"])


def test_vec_vae_generate_shape():
    cfg = _tiny_vec_cfg()
    model = BERTSetVAE(cfg).eval()
    z = torch.randn(3, cfg.latent_dim)
    ids = model.generate(z, temperature=1.0)
    assert ids.shape == (3, cfg.max_len)


def test_set_vae_forward_shapes():
    cfg = _tiny_set_cfg()
    model = BERTSetVAEV2(cfg).eval()
    B, T = 2, cfg.max_len
    tokens = torch.randint(6, cfg.vocab_size, (B, T))
    mask = torch.ones(B, T, dtype=torch.bool)
    out = model(tokens, mask, beta=0.01, sample=True)
    assert out["logits"].shape == (B, T, cfg.vocab_size)
    assert out["mu"].shape == (B, cfg.n_queries, cfg.latent_dim)
    assert out["z_tokens"].shape == (B, cfg.n_queries, cfg.latent_dim)
    assert torch.isfinite(out["loss"])


def test_set_vae_generate_shape():
    cfg = _tiny_set_cfg()
    model = BERTSetVAEV2(cfg).eval()
    z = torch.randn(3, cfg.n_queries, cfg.latent_dim)
    ids = model.generate(z, temperature=1.0)
    assert ids.shape == (3, cfg.max_len)


def test_set_vae_kernel_size_one_matches_freeze():
    cfg = _tiny_set_cfg()
    assert cfg.kernel_size == 1
    model = BERTSetVAEV2(cfg)
    # pointwise conv: weight shape [out, in, k]
    k = model.dec_blocks[0].conv1.kernel_size[0]
    assert k == 1
