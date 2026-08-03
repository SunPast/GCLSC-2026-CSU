# Stage 1 scripts

## Unconditional generation

`-o` / `--out` is **required**.

```bash
# from DiffEboAb/
python scripts/stage1_vae/generate.py -o outputs/uncond.fasta --n 1000
python scripts/stage1_vae/generate.py -o out.fasta --ckpt checkpoints/stage1_vae/vec.pt --n 100
```

Loads `checkpoints/stage1_vae/model.pt` by default (prior \(z\sim\mathcal{N}(0,I)\), decode).
