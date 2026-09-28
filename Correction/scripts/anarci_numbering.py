#!/usr/bin/env python3
"""
IMGT numbering via ANARCI + pyhmmer (no hmmscan binary needed)
================================================================
ANARCI requires the HMMER3 `hmmscan` command-line binary, which is not
available on Windows. This module replaces that subprocess call with pyhmmer
(pure Python bindings to the same HMMER3 libraries) and reuses ANARCI's own
`number_imgt` scheme for the actual numbering.

The result is exact IMGT numbering that does NOT depend on sequence motifs, so
it correctly locates conserved positions (C23, W41, C104, F/W118, G119) even
when a mutation destroys the local motif (e.g. G119→E, double-Cys breaks) —
the failure mode that the motif-based detector in pipeline.py misses.

Usage:
    from anarci_numbering import imgt_number
    num = imgt_number("EVQLVES...")   # -> {23:'C', 41:'W', 104:'C', 118:'W', 119:'G', ...}
    # or None if the sequence is not a recognizable antibody V domain

License note: ANARCI is BSD-3-Clause; pyhmmer is MIT; both commercial-friendly.
"""

import os
import re

import anarci
import pyhmmer
from anarci.schemes import number_imgt
from pyhmmer.easel import Alphabet, TextSequence
from pyhmmer.plan7 import HMMPressedFile, Pipeline

_HMM_PATH = os.path.join(os.path.dirname(anarci.__file__), "dat", "HMMs", "ALL.hmm")
_ALPHABET = Alphabet.amino()


def _run_scan(seq):
    """Run pyhmmer hmmscan-equivalent for one sequence; return best VH/VL hit.

    The HMMPressedFile is re-opened per call (it is not re-entrant, and reusing
    a module-level instance hangs on the second `with`).
    """
    s = TextSequence(name="query", sequence=seq).digitize(_ALPHABET)
    pipeline = Pipeline(alphabet=_ALPHABET)
    with HMMPressedFile(_HMM_PATH) as targets:
        top = pipeline.scan_seq(s, targets)
    best = None
    for hit in top:
        if hit.evalue < 1e-10 and "_" in hit.name:
            # prefer heavy chain H, else kappa K / lambda L
            if best is None or ("_H" in hit.name and "_H" not in best.name):
                best = hit
    return best


def _hit_to_state_vector(hit):
    """Convert a pyhmmer hit into ANARCI's state_vector format."""
    al = list(hit.domains)[0].alignment
    hmm = al.hmm_sequence
    tgt = al.target_sequence
    hmm_states = list(range(1, 129))
    state_vector = []
    h = al.hmm_from - 1   # 0-based match-state index
    s = 0                 # 0-based index into the query sequence
    for i in range(len(hmm)):
        t = tgt[i]
        hm = hmm[i]
        if t != "-" and hm != "-":
            st = "m"       # match state (query residue aligned to HMM match state)
        elif t != "-" and hm == "-":
            st = "i"       # insert state (query residue is an insertion)
        elif t == "-" and hm != "-":
            st = "d"       # delete state (HMM match state with no query residue)
        else:
            continue
        if h >= len(hmm_states):
            break           # alignment extends past the IMGT model (e.g. constant region)
        state_vector.append(((hmm_states[h], st), None if st == "d" else s))
        if st == "m":
            h += 1
            s += 1
        elif st == "i":
            s += 1
        else:
            h += 1
    return state_vector


def imgt_number(seq, chain="H"):
    """
    Return IMGT numbering for an antibody V-domain sequence.

    Returns dict {imgt_position: amino_acid} for the best-scoring domain
    (heavy chain preferred), or None if the sequence is not recognized.

    Conserved positions of interest: 23 (Cys), 41 (Trp), 104 (Cys),
    118 (F/W), 119 (Gly).
    """
    if len(seq) < 60 or len(seq) > 300:
        return None
    try:
        hit = _run_scan(seq)
    except Exception:
        return None
    if hit is None:
        return None
    try:
        state_vector = _hit_to_state_vector(hit)
        numbering, start, end = number_imgt(state_vector, seq)
    except Exception:
        return None
    result = {}
    for (num, insertion), aa in numbering:
        # Only keep the primary position (ignore insertions like 111A)
        if num not in result:
            result[num] = aa
    return result


def imgt_positions(seq):
    """
    Return {IMGT_position: sequence_index} for match states (0-based index into
    the query sequence). This is what the correction pipeline needs to edit the
    right residue. Returns None if the sequence is not recognized.
    """
    if len(seq) < 60 or len(seq) > 300:
        return None
    try:
        hit = _run_scan(seq)
    except Exception:
        return None
    if hit is None:
        return None
    al = list(hit.domains)[0].alignment
    hmm = al.hmm_sequence
    tgt = al.target_sequence
    hmm_states = list(range(1, 129))
    positions = {}
    h = al.hmm_from - 1
    s = 0
    for i in range(len(hmm)):
        t = tgt[i]
        hm = hmm[i]
        if t != "-" and hm != "-":
            # match state: query residue at s maps to IMGT position hmm_states[h]
            if h >= len(hmm_states):
                break
            positions[hmm_states[h]] = s
            h += 1
            s += 1
        elif t != "-" and hm == "-":
            # insert state: query residue advances, no IMGT match position
            s += 1
        elif t == "-" and hm != "-":
            # delete state: HMM match state advances, no query residue
            h += 1
    return positions


def conserved_residues(seq):
    """
    Extract the five conserved residues used by the correction pipeline.
    Returns dict with keys C23, W41, C104, F118, G119 (values are the amino
    acid at that IMGT position, or None).
    """
    pos = imgt_positions(seq)
    if pos is None:
        return None
    return {
        "C23": seq[pos[23]] if 23 in pos else None,
        "W41": seq[pos[41]] if 41 in pos else None,
        "C104": seq[pos[104]] if 104 in pos else None,
        "F118": seq[pos[118]] if 118 in pos else None,
        "G119": seq[pos[119]] if 119 in pos else None,
    }


if __name__ == "__main__":
    # self-test: a clean VH, then the same VH with G119 mutated to Glu.
    # The mutation destroys the J-motif that motif-based detection relies on,
    # so W41/F118/G119 must still be located correctly here.
    seq = ("EVQLVESGGGLVQPGGSLRLSCAASGFTFTDYTMDWVRQAPGKGLEWVADVNPNSGGSIYNQRFKGRFTLSVDRSKNTLYLQMNSLRAEDTAVYYCARNLGPSFYFDYWGQGTLVTVSS")
    print("clean VH   :", conserved_residues(seq))
    m = re.search(r"WGQG", seq)
    p = m.start() + 1
    mut = seq[:p] + "E" + seq[p + 1:]
    print("G119->E    :", conserved_residues(mut))
