#!/usr/bin/env python3
"""
IMGT-annotated FASTA writer for VH sequences
=============================================
Renders every input VH sequence onto the fixed IMGT V-domain frame
(positions 1-128) and pads positions with no residue using '.'.

Why a frame-shift repair step is needed
---------------------------------------
ANARCI (via pyhmmer) assigns IMGT numbering by HMM alignment. On sequences
whose N-terminus is truncated, the alignment can slip by one match state: a
controlled test that deletes the first k residues from one full-length VH
gives a correct frame for k = 14, 16, 17, 18 but a frame shifted by +1 for
k = 15 -- after which IMGT 23 reads 'S' and IMGT 104 reads 'Y' even though
both conserved cysteines are physically present.

This is a numbering artefact, not a mutation. It is repaired here by
re-anchoring the frame on the two conserved cysteines (IMGT 23 and IMGT 104,
Chothia & Lesk 1987; Lefranc et al. 2003), which must both be Cys in a
foldable V domain.

Conventions in the output
-------------------------
- '.'  = IMGT position with no residue in this sequence. This covers both
         short-loop gap positions (e.g. CDR1/FR3 length variation, which is
         normal) and genuinely absent regions (e.g. a deleted N-terminus).
- Residues that fall past IMGT 128 are appended after column 128 and counted
  in the header as Cterm_extra.

Usage:
    python scripts/annotate_imgt.py -i D2_candidates_50.fa -o D2_candidates_50_imgt.fa
"""

import argparse
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from anarci_numbering import _hit_to_state_vector, _run_scan
from anarci.schemes import number_imgt

IMGT_MAX = 128
# The conserved 1st/2nd cysteine pair sits at IMGT 23 and IMGT 104 in every
# V domain; C23 and C104 are 81 IMGT positions apart.
C23_IMGT = 23
C104_IMGT = 104
CYS_SPAN = C104_IMGT - C23_IMGT


def read_fasta(path):
    entries, header, chunks = [], None, []
    with open(path, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header:
                    entries.append((header, "".join(chunks)))
                header, chunks = line[1:], []
            else:
                chunks.append(line)
    if header:
        entries.append((header, "".join(chunks)))
    return entries


def _raw_alignment(seq):
    """ANARCI IMGT alignment as [(column, residue)] plus residues past col 128.

    `number_imgt` emits exactly 128 columns. Residues the aligner placed past
    the last IMGT match state are dropped there, so they are recovered here
    from the sequence tail (verified to be an unmodified suffix in practice).
    """
    hit = _run_scan(seq)
    if hit is None:
        return None
    state_vector = _hit_to_state_vector(hit)
    numbering, _, _ = number_imgt(state_vector, seq)
    cols = [(num, aa) for (num, _ins), aa in numbering if aa != "-"]
    used = "".join(aa for _, aa in cols)
    if not seq.startswith(used):
        # Aligner reordered residues relative to the query; not safe to repair.
        return None
    trailing = seq[len(used):]
    cols.extend((IMGT_MAX + 1 + i, aa) for i, aa in enumerate(trailing))
    return cols


def _render(cols, shift):
    """Place [(column, residue)] shifted by `shift` onto the IMGT frame."""
    width = IMGT_MAX + max(0, max(c - shift for c, _ in cols) - IMGT_MAX)
    frame = ["."] * width
    for col, aa in cols:
        new = col - shift
        if new < 1:
            return None          # a residue pushed off the N-terminus
        if new > width:
            return None          # a residue pushed past the C-terminus
        if frame[new - 1] != ".":
            return None          # two residues in one column
        frame[new - 1] = aa
    return "".join(frame)


def _cys_anchored(cols):
    """Pick the shift placing the conserved Cys pair on IMGT 23 / 104."""
    for shift in (0, 1, -1, 2, -2, 3, -3):
        frame = _render(cols, shift)
        if frame is None:
            continue
        if frame[C23_IMGT - 1] == "C" and frame[C104_IMGT - 1] == "C":
            return frame, shift
    return None, None


def annotate(seq):
    """Return (aligned_string, meta_dict) or (None, meta) if not numberable."""
    cols = _raw_alignment(seq)
    if not cols:
        return None, {"error": "ANARCI did not recognise the sequence"}

    frame, shift = _cys_anchored(cols)
    anchored = True
    if frame is None:
        # No shift puts Cys on 23/104 -- the cysteines are genuinely broken
        # (e.g. a C23 or C104 mutation). Keep ANARCI's frame unmodified.
        frame = _render(cols, 0)
        shift, anchored = 0, False

    present = [i + 1 for i, ch in enumerate(frame[:IMGT_MAX]) if ch != "."]
    first = min(present) if present else None
    meta = {
        "imgt_first": first,
        "imgt_last": max(present) if present else None,
        # Counted in IMGT position numbers, i.e. IMGT 1..first-1. This is the
        # natural unit for an IMGT annotation and is directly readable off the
        # IMGT=<span> field.
        "nterm_missing": (first - 1) if first else None,
        "cterm_extra": len(frame) - IMGT_MAX,
        "frame_shift": shift,
        "cys_anchored": anchored,
        "c23": frame[C23_IMGT - 1],
        "c104": frame[C104_IMGT - 1],
    }
    return frame, meta


def main():
    ap = argparse.ArgumentParser(description="IMGT-annotated FASTA ('.' = absent)")
    ap.add_argument("--input", "-i", required=True)
    ap.add_argument("--output", "-o", required=True)
    args = ap.parse_args()

    entries = read_fasta(args.input)
    out, stats = [], {"ok": 0, "fail": 0, "shifted": 0, "nterm": 0}
    for header, seq in entries:
        frame, meta = annotate(seq)
        sid = header.split("|")[0]
        if frame is None:
            stats["fail"] += 1
            out.append(f">{header}|IMGT=NA")
            out.append(seq)
            continue
        stats["ok"] += 1
        if meta["frame_shift"]:
            stats["shifted"] += 1
        if meta["nterm_missing"]:
            stats["nterm"] += 1
        ann = (
            f"|IMGT={meta['imgt_first']}-{meta['imgt_last']}"
            f"|Nterm_missing={meta['nterm_missing']}"
            f"|Cterm_extra={meta['cterm_extra']}"
            f"|frame_shift={meta['frame_shift']:+d}"
            f"|C23={meta['c23']}|C104={meta['c104']}"
        )
        out.append(f">{header}{ann}")
        # One unwrapped line per sequence: the IMGT frame is fixed-width, so
        # keeping it on a single line makes the columns line up in any editor.
        out.append(frame)

    with open(args.output, "w", encoding="utf-8") as f:
        f.write("\n".join(out) + "\n")

    print(f"Input  : {args.input}  ({len(entries)} sequences)")
    print(f"Output : {args.output}")
    print(f"  numbered      : {stats['ok']}")
    print(f"  frame repaired: {stats['shifted']}")
    print(f"  N-term missing: {stats['nterm']}")
    print(f"  unrecognised  : {stats['fail']}")


if __name__ == "__main__":
    main()
