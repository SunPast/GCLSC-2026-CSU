#!/usr/bin/env python3
"""
Generate example input for the Module 2 pipeline.
=================================================
Writes 25 sequences into the pipeline's expected input location so a fresh
checkout has something to run:

    data/01_generated/BATCH_HER2_20260809_01/Abseqs_BATCH_HER2_20260809_01.fa

The set deliberately covers every defect class the pipeline handles: clean
sequences, non-standard residues, each of the 5 conserved-site mutations,
double-Cys, length anomalies, quality motifs (N-glyc, deamidation/isomerization
hot spots), and an all-sites-broken sequence that must FAIL. Each header carries
a `type=` label naming the intended defect.

Template positions (0-indexed, verified):
  VH_323: C23=21, W41=35, C104=95, J-motif=WGQG(100-103) -> F118=100, G119=101
  VH_169: C23=21, W41=35, C104=95, J-motif=WGQG(100-103) -> F118=100, G119=101
  VL_K39: C23=22, W41=34, C104=87, J-motif=FGQG(97-100)  -> F118=97,  G119=98
  VL_L44: C23=21, W41=35, C104=88, J-motif=FGGG(100-103) -> F118=100, G119=101

Usage:
    python scripts/generate_mock.py
"""

import os
import random

# Anchor output to the project root so the script works from any directory.
ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

random.seed(42)

# === Templates (germline-based, realistic) ===
T_VH_323 = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGST"
    "YYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDYWGQGTLVTVSS"
)  # 111 AA, C23=21, W41=35, C104=95, F118≈100(W), G119≈101

T_VH_169 = (
    "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGGIIPIFGT"
    "ANYAQKFQGRVTITADESTSTAYMELSSLRSEDTAVYYCARDYWGQGTTVTVSS"
)  # 111 AA, C23=21, W41=35, C104=95, F118≈100(W), G119≈101

T_VL_K39 = (
    "DIQMTQSPSSLSASVGDRVTITCRASQSISSYLNWYQQKPGKAPKLLIYAASSLQS"
    "GVPSRFSGSGSGTDFTLTISSLQPEDFATYYCQQSYSTPYTFGQGTKLEIK"
)  # 107 AA, C23=22, W41=34, C104=87, F118≈97(F), G119≈98

T_VL_L44 = (
    "QSVLTQPPSVSAAPGQKVTISCSGSSSNIGNNYVSWYQQLPGTAPKLLIYDNNKRPS"
    "GIPDRFSGSKSGTSATLGITGLQTGDEADYYCGTWDSSLSAVVFGGGTKLTVL"
)  # 110 AA, C23=21, W41=35, C104=88, F118≈100(F), G119≈101


def mutate_at(seq, pos, new_aa):
    """Replace AA at position pos (0-indexed)."""
    if pos >= len(seq):
        return seq  # safety
    return seq[:pos] + new_aa + seq[pos + 1:]


def mk_header(seq_id, score, length, extras=""):
    """Create FASTA header line."""
    h = f">{seq_id}|Score={score:.2f}|Batch=BATCH_HER2_20260809_01|Length={length}"
    if extras:
        h += f"|{extras}"
    return h


def fold_seq(seq):
    """Fold to 80 chars per line."""
    return "\n".join(seq[i:i + 80] for i in range(0, len(seq), 80))


def make_fasta(seq_id, seq, score, extras=""):
    return f"{mk_header(seq_id, score, len(seq), extras)}\n{fold_seq(seq)}"


def main():
    batch_id = "BATCH_HER2_20260809_01"

    def sid(n):
        return f"{batch_id}_SEQ_{n:05d}"

    entries = []
    n = 0

    # ── Group 1: Clean sequences (5) ──
    n += 1; entries.append(make_fasta(sid(n), T_VH_323, 0.92, "type=clean_VH"))
    n += 1; entries.append(make_fasta(sid(n), T_VH_169, 0.88, "type=clean_VH"))
    n += 1; entries.append(make_fasta(sid(n), T_VL_K39, 0.85, "type=clean_VL"))
    n += 1; entries.append(make_fasta(sid(n), T_VL_L44, 0.83, "type=clean_VL"))
    n += 1; entries.append(make_fasta(sid(n), T_VH_323, 0.97, "type=clean_VH_dup"))

    # ── Group 2: Non-standard AA (5) ──
    n += 1; entries.append(make_fasta(sid(n), mutate_at(T_VH_323, 70, "B"), 0.45, "type=nonstd_B_FR3"))
    n += 1; entries.append(make_fasta(sid(n), mutate_at(T_VH_323, 55, "Z"), 0.40, "type=nonstd_Z_CDR2"))
    n += 1; entries.append(make_fasta(sid(n), mutate_at(T_VH_323, 5, "J"), 0.42, "type=nonstd_J_FR1"))
    seq_X = mutate_at(mutate_at(T_VH_323, 30, "X"), 80, "X")
    n += 1; entries.append(make_fasta(sid(n), seq_X, 0.30, "type=nonstd_XX"))
    n += 1; entries.append(make_fasta(sid(n), mutate_at(T_VH_323, 40, "U"), 0.38, "type=nonstd_U_FR2"))

    # ── Group 3: Conserved position mutations (5) ──
    # C23→S (pos 21 in VH_323)
    n += 1; entries.append(make_fasta(sid(n), mutate_at(T_VH_323, 21, "S"), 0.35, "type=C23_S"))
    # C104→Y (pos 95 in VH_323)
    n += 1; entries.append(make_fasta(sid(n), mutate_at(T_VH_323, 95, "Y"), 0.33, "type=C104_Y"))
    # W41→L (pos 35 in VH_323)
    n += 1; entries.append(make_fasta(sid(n), mutate_at(T_VH_323, 35, "L"), 0.50, "type=W41_L"))
    # F118→A (pos 97 in VL_K39 = F)
    n += 1; entries.append(make_fasta(sid(n), mutate_at(T_VL_K39, 97, "A"), 0.48, "type=F118_A"))
    # G119→E (pos 101 in VH_323)
    n += 1; entries.append(make_fasta(sid(n), mutate_at(T_VH_323, 101, "E"), 0.25, "type=G119_E"))

    # ── Group 4: Severe/multiple issues (5) ──
    # Double C mutant (C23→A, C104→G in VH_323)
    seq_2c = mutate_at(mutate_at(T_VH_323, 21, "A"), 95, "G")
    n += 1; entries.append(make_fasta(sid(n), seq_2c, 0.10, "type=double_C_mut"))
    # Too short
    n += 1; entries.append(make_fasta(sid(n),
        "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGR", 0.05, "type=too_short"))
    # Too long (VH_323 + 5×GGGGS)
    n += 1; entries.append(make_fasta(sid(n), T_VH_323 + "GGGGS" * 5, 0.08, "type=too_long"))
    # Mixed: nonstd + conserved mutation
    seq_mix = mutate_at(mutate_at(mutate_at(T_VH_323, 21, "A"), 50, "B"), 75, "X")
    n += 1; entries.append(make_fasta(sid(n), seq_mix, 0.15, "type=mixed_issues"))
    # Free cysteine (odd count: introduce extra C at pos 60)
    n += 1; entries.append(make_fasta(sid(n), mutate_at(T_VH_323, 60, "C"), 0.55, "type=free_cysteine"))

    # ── Group 5: Edge cases (5) ──
    # N-glycosylation motif (N-X-S) in CDR2: mutate pos 52→N, 54→S
    seq_ngly = mutate_at(mutate_at(T_VH_323, 52, "N"), 54, "S")
    n += 1; entries.append(make_fasta(sid(n), seq_ngly, 0.65, "type=N_glyc_CDR2"))
    # DG isomerization in CDR
    seq_dg = mutate_at(mutate_at(T_VH_323, 55, "D"), 56, "G")
    n += 1; entries.append(make_fasta(sid(n), seq_dg, 0.60, "type=DG_isomerization"))
    # VL with W41→F (conservative, might still pass)
    n += 1; entries.append(make_fasta(sid(n), mutate_at(T_VL_K39, 34, "F"), 0.58, "type=W41_F_VL"))
    # VL with C23→T
    n += 1; entries.append(make_fasta(sid(n), mutate_at(T_VL_L44, 21, "T"), 0.36, "type=C23_T_VL"))
    # All conserved positions broken (should FAIL)
    seq_all = T_VH_323
    seq_all = mutate_at(seq_all, 21, "S")   # C23→S
    seq_all = mutate_at(seq_all, 35, "L")   # W41→L
    seq_all = mutate_at(seq_all, 95, "Y")   # C104→Y
    seq_all = mutate_at(seq_all, 101, "A")  # G119→A
    n += 1; entries.append(make_fasta(sid(n), seq_all, 0.02, "type=all_conserved_mut"))

    # Write output
    outdir = os.path.join(ROOT, "data", "01_generated", batch_id)
    os.makedirs(outdir, exist_ok=True)
    outpath = os.path.join(outdir, f"Abseqs_{batch_id}.fa")
    with open(outpath, "w", encoding="utf-8") as f:
        f.write("\n".join(entries) + "\n")

    print(f"Generated {n} mock sequences -> {outpath}")
    for i, e in enumerate(entries):
        hdr = e.split("\n")[0]
        # Extract type from end of header
        parts = hdr.split("|")
        extra = parts[-1] if len(parts) >= 5 else ""
        seq_id = parts[0][1:]
        score = parts[1].split("=")[1]
        length = parts[3].split("=")[1]
        print(f"  {i+1:2d}. {seq_id:40s} Score={score:>5s}  Len={length:>3s}  {extra}")


if __name__ == "__main__":
    main()
