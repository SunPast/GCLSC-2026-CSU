#!/usr/bin/env python3
"""Unconditional antibody generation from the Stage-1 VAE prior.

Samples z ~ N(0, I), then decodes with the frozen non-autoregressive decoder.

Usage (cwd = DiffEboAb/):
    python scripts/stage1_vae/generate.py -o outputs/uncond.fasta --n 1000
    python scripts/stage1_vae/generate.py -o out.fasta --ckpt checkpoints/stage1_vae/vec.pt --n 100

Default checkpoint slot:
    checkpoints/stage1_vae/model.pt   # → set.pt
"""
from __future__ import annotations

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import torch

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

from src.data.tokenizer import build_tokenizer
from src.models.bert_set_vae import BERTSetVAE
from src.models.bert_set_vae_v2 import BERTSetVAEConfigV2, BERTSetVAEV2


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Unconditional VAE prior sampling")
    p.add_argument(
        "-o", "--out",
        required=True,
        help="output FASTA path (required)",
    )
    p.add_argument(
        "--ckpt",
        default="checkpoints/stage1_vae/model.pt",
        help="checkpoint path (default: stable slot model.pt)",
    )
    p.add_argument("--n", type=int, default=1000, help="number of sequences")
    p.add_argument("--batch_size", type=int, default=256)
    p.add_argument("--temperature", type=float, default=1.0)
    p.add_argument("--seed", type=int, default=42)
    p.add_argument(
        "--device",
        default="cuda" if torch.cuda.is_available() else "cpu",
    )
    return p.parse_args(argv)


def load_model(ckpt_path: Path, device: torch.device):
    ckpt = torch.load(ckpt_path, map_location="cpu", weights_only=False)
    cfg = ckpt["config"]
    if isinstance(cfg, BERTSetVAEConfigV2):
        model = BERTSetVAEV2(cfg)
        kind = "BERTSetVAEV2"
    else:
        model = BERTSetVAE(cfg)
        kind = "BERTSetVAE"
    model.load_state_dict(ckpt["model"])
    model.to(device).eval()
    return model, cfg, kind, ckpt.get("step")


def sample_prior(cfg, batch: int, device: torch.device) -> torch.Tensor:
    """z ~ N(0, I). set-latent (B, Q, D) or vector latent (B, D)."""
    if isinstance(cfg, BERTSetVAEConfigV2):
        return torch.randn(batch, cfg.n_queries, cfg.latent_dim, device=device)
    return torch.randn(batch, cfg.latent_dim, device=device)


def main(argv: list[str] | None = None) -> None:
    args = parse_args(argv)
    torch.manual_seed(args.seed)
    np.random.seed(args.seed)

    device = torch.device(args.device)
    if device.type == "cuda":
        torch.backends.cudnn.benchmark = True

    ckpt_path = Path(args.ckpt)
    if not ckpt_path.exists():
        raise FileNotFoundError(
            f"Checkpoint not found: {ckpt_path}\n"
            f"Place weights in checkpoints/stage1_vae/ (see checkpoints/README.md)."
        )

    tok = build_tokenizer(max_len=140)
    print(f"Loading {ckpt_path} on {device}...")
    model, cfg, kind, step = load_model(ckpt_path, device)
    if isinstance(cfg, BERTSetVAEConfigV2):
        lat_str = f"{cfg.n_queries}×{cfg.latent_dim}"
    else:
        lat_str = f"{cfg.latent_dim}"
    print(f"  {kind}  step={step}  latent={lat_str}")

    all_seqs: list[str] = []
    t0 = time.time()
    for i in range(0, args.n, args.batch_size):
        b = min(args.batch_size, args.n - i)
        z = sample_prior(cfg, b, device)
        with torch.no_grad():
            tokens = model.generate(z, temperature=args.temperature)
        for row in tokens.cpu().numpy():
            all_seqs.append(tok.decode(row, skip_special=True, stop_at_eos=True))
        done = i + b
        if done % max(args.batch_size, 1) == 0 or done == args.n:
            print(f"  {done}/{args.n}  ({time.time() - t0:.1f}s)", flush=True)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    tag = ckpt_path.stem
    with open(out_path, "w") as f:
        for idx, seq in enumerate(all_seqs):
            f.write(f">uncond_{tag}_{idx:05d} len={len(seq)}\n{seq}\n")

    lengths = [len(s) for s in all_seqs]
    print(f"Done. {len(all_seqs)} sequences → {out_path.resolve()}")
    print(f"Sample: {all_seqs[0][:70]}...")
    print(
        f"Length: min={min(lengths)} max={max(lengths)} "
        f"mean={float(np.mean(lengths)):.1f}"
    )


if __name__ == "__main__":
    main()
