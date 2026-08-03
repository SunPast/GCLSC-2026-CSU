# Checkpoint slots 

Inference scripts **must** load weights from here — not from local training dumps.

```
checkpoints/
├── stage1_vae/
│   ├── model.pt   → set.pt     # default for generate / downstream
│   ├── set.pt                  # set-latent VAE (64×512), pipeline freeze
│   └── vec.pt                  # vector-latent VAE (single z)
├── stage2_diffusion/           # later
└── stage3_rl/                  # later
```

Naming by latent geometry 
- **vec** — single-vector latent \(z\in\mathbb{R}^{D}\)
- **set** — set-latent \(z\in\mathbb{R}^{Q\times D}\) + position cross-attn


## Default path used by code

```text
checkpoints/stage1_vae/model.pt
```
