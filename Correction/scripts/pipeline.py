#!/usr/bin/env python3
"""
Module 2: Antibody VH sequence QC → correction → audit pipeline
===============================================================
Three steps per sequence:

  a. QC         — non-standard residues, the 5 conserved sites, length, quality motifs
  b. Correction — repair what is repairable (non-standard residues, conserved sites)
  c. Audit      — score the corrected sequence, decide PASS or FAIL

Usage:
    python scripts/pipeline.py -b <BatchID>                     # reads the default input path
    python scripts/pipeline.py -b <BatchID> -i <path/to/in.fa>  # explicit input
    python scripts/pipeline.py -i <path/to/in.fa>               # batch id from the file name

Input : data/01_generated/<BatchID>/Abseqs_<BatchID>.fa   (or -i)
Output: data/02_correction/<BatchID>/manifest.csv
        data/02_correction/<BatchID>/batch_correction_summary/
        data/02_correction/<BatchID>/<Seq_ID>/reports/structure_report.md

Rules and their literature basis: docs/02_规则说明_中文.md
"""

import argparse
import csv
import json
import os
import re
import sys
from datetime import datetime

# Optional: exact IMGT numbering via ANARCI (pyhmmer backend). Falls back to
# motif-based detection if unavailable.
try:
    sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
    from anarci_numbering import imgt_positions as _imgt_positions
    _HAS_ANARCI = True
except Exception:
    _HAS_ANARCI = False

# ============================================================
# Configuration
# ============================================================
ROOT_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

# Standard proteinogenic amino acids
STANDARD_AA = set("ACDEFGHIKLMNPQRSTVWY")

# Non-standard → standard mapping (literature-based)
# Ref: IUPAC-IUB (1984) Eur. J. Biochem. 138(1):9-37
NONSTD_MAP = {
    "B": "N",   # Asx → Asn
    "Z": "E",   # Glx → Glu
    "J": "L",   # Xle → Leu
    "U": "C",   # Sec → Cys
    "O": "K",   # Pyl → Lys
    "X": "A",   # Unknown → Ala (safest fallback)
    "*": None,  # Stop → remove
    "-": None,  # Gap → remove
}

# ============================================================
# Motif-based conserved position identification
# ============================================================
# Since full IMGT numbering requires ANARCI, we use motif-based
# approximation with cascading fallbacks. Ref: Dunbar & Deane (2016)
#
# Strategy: Try strict motif first; if failed (due to mutations), use
# permissive patterns; if still failed, estimate from structural constraints.

# Primary FR2 W-motif: W-[VILFY]-[RKQ] identifies conserved W41 across species
# AND chain types. Covers: human VH WVRQ, mouse VH WVKQ, rabbit VH WIRK,
# variants WVKL, and light-chain VL WYQQ. The 2nd position is a hydrophobic
# (V/I/L/Y/F), the 3rd is basic R/K or Q — this broad pattern is the real
# conserved signature of FR2 (a v3 real-data finding; the old W-X-[RQ]-X-X-X-G
# only matched human WVRQAPG).
FR2_W_STRICT = re.compile(r"W[VILFY][RKQ]")
# Fallback: allow a hydrophobic substitution AT the W41 position itself
# (L/I/F/V/Y can occur when W41 is mutated — needed for the correction path).
FR2_W_FALLBACK = re.compile(r"[WLIFVY][VILFY][RKQ]")
# Tertiary fallback: look for [R,Q]-X-X-P-G pattern near position 35-40
FR2_RQPG = re.compile(r"[RQ].[PG]")

# Primary FR4 J-motif: [FW]-G-X-G (identifies conserved F/W118 and G119)
FR4_J_STRICT = re.compile(r"[FW]G.G")
# Fallback: allow substitutions in F118 position + more permissive
FR4_J_FALLBACK = re.compile(r"[FWLIMAVY]G.G")
# Tertiary: just find GXG near C-terminus (last 20 residues)
FR4_J_GXG = re.compile(r"G.G")


def find_conserved_positions(seq):
    """
    Identify conserved positions in an antibody V domain sequence
    using cascading motif detection with structural heuristics.

    Three-tier detection for each landmark:
      1. Strict motif (native sequence)
      2. Permissive motif (allow conservative mutations)
      3. Structural estimate (distance from other landmarks or sequence ends)

    Typical IMGT spacing in raw residues:
      C23 ← ~14 AA → W41 ← ~55 AA → C104 ← ~6 AA → F118-G119 ← ~10 AA → C-term
    """

    result = {
        "C23_pos": None, "C23_is_C": False,
        "C104_pos": None, "C104_is_C": False,
        "W41_pos": None, "W41_is_W": False,
        "F118_pos": None, "F118_is_FW": False,
        "G119_pos": None, "G119_is_G": False,
        "domain_type": "unknown",
        # Confidence flags: True only if found by motif (Tier 1/2), False if
        # estimated from structural heuristics (Tier 3 / last resort). Used to
        # reject non-antibody sequences that lack the signature motifs.
        "w41_confident": False,
        "j_confident": False,
    }
    n = len(seq)

    # Prefer exact IMGT numbering via ANARCI (motif-independent, robust to
    # G119→E / double-Cys mutations that break the motifs below).
    if _HAS_ANARCI:
        try:
            apos = _imgt_positions(seq)
        except Exception:
            apos = None
        if apos:
            c23 = apos.get(23)
            c104 = apos.get(104)
            w41 = apos.get(41)
            f118 = apos.get(118)
            g119 = apos.get(119)
            # Motif correction for FR4: ANARCI's IMGT numbering drifts on long
            # CDR3 (>13 match states), mis-placing F118/G119 into the CDR3 tail.
            # When the FR4 J-motif ([FW]-G-X-G) is intact, trust the motif for
            # F118/G119; when it is broken (a real G119 mutation), keep ANARCI's
            # detection. This is the motif↔ANARCI complement.
            # Use strict first, then fallback ([FWLIMAVY]G.G) so that an F118→A
            # mutation (which breaks [FW]G.G into AGQG) is still corrected.
            jm = FR4_J_STRICT.search(seq)
            if jm is None:
                jm = FR4_J_FALLBACK.search(seq)
            if jm and jm.start() > n - 20:
                j_f = jm.start()
                j_g = jm.start() + 1
                if seq[j_g] == "G":
                    f118 = j_f
                    g119 = j_g
            result["C23_pos"] = c23
            result["C23_is_C"] = (c23 is not None and c23 < n and seq[c23] == "C")
            result["C104_pos"] = c104
            result["C104_is_C"] = (c104 is not None and c104 < n and seq[c104] == "C")
            # C23/C104 motif correction: ANARCI drifts by ~1 position on rare
            # V genes. If ANARCI says non-C but a Cys exists at the expected
            # FR1/FR3 position, trust the motif (same complement as FR4 above).
            all_c = [i for i, aa in enumerate(seq) if aa == "C"]
            if not result["C23_is_C"]:
                cand = [p for p in all_c if p < 40]
                if cand:
                    c23 = cand[0]
                    result["C23_pos"] = c23
                    result["C23_is_C"] = True
            if not result["C104_is_C"]:
                cand = [p for p in all_c if p != c23 and p > 40]
                if cand:
                    c104 = cand[-1]
                    result["C104_pos"] = c104
                    result["C104_is_C"] = True
            result["W41_pos"] = w41
            result["W41_is_W"] = (w41 is not None and w41 < n and seq[w41] == "W")
            result["F118_pos"] = f118
            result["F118_is_FW"] = (f118 is not None and f118 < n and seq[f118] in ("F", "W"))
            result["G119_pos"] = g119
            result["G119_is_G"] = (g119 is not None and g119 < n and seq[g119] == "G")
            # W41 / F118 motif correction (same complement as C23/C104/G119):
            # ANARCI drifts on rare V genes / long CDR3; when the FR2/FR4 motif
            # is intact, trust the motif for these two positions.
            if not result["W41_is_W"]:
                m = FR2_W_STRICT.search(seq)
                if m:
                    result["W41_pos"] = m.start()
                    result["W41_is_W"] = True
            if not result["F118_is_FW"]:
                jm2 = FR4_J_STRICT.search(seq)
                if jm2 is None:
                    jm2 = FR4_J_FALLBACK.search(seq)
                if jm2 and jm2.start() > n - 20:
                    result["F118_pos"] = jm2.start()
                    result["F118_is_FW"] = True
            result["w41_confident"] = True
            result["j_confident"] = True
            return result

    # ================================================================
    # Step 0: C23 — 1st-CYS in strand B (FR1), an INDEPENDENT anchor.
    # The first Cys in the sequence (raw pos 20-26). Locate this FIRST so W41
    # can be found by distance from C23, not by an over-broad motif that
    # matched spurious "LVQ" etc. in real data.
    # ================================================================
    all_c = [i for i, aa in enumerate(seq) if aa == "C"]
    c23_candidates = [p for p in all_c if p < 40]
    if c23_candidates:
        result["C23_pos"] = c23_candidates[0]
        result["C23_is_C"] = True
    else:
        # C23 may be mutated (no C in FR1). Leave None for now; W41 will be
        # located by motif and we back-fill C23 = W41 - 14 afterwards.
        result["C23_pos"] = None
        result["C23_is_C"] = False

    # ================================================================
    # Step 1: W41 — conserved Trp in FR2, located by DISTANCE from C23
    # (W41 ≈ C23 + 14). Motifs are used only to confirm within the window,
    # not to search the whole sequence (whole-seq search matched spurious
    # "LVQ"/"LVRQ" in real data and cascaded into wrong C23).
    # ================================================================
    w41_pos = None
    w41_is_w = False
    w41_confident = False

    c23 = result["C23_pos"]
    if c23 is not None:
        w41_expected = c23 + 14
        # Search ±6 window for W
        for offset in range(0, 7):
            for sign in [1, -1]:
                cand = w41_expected + sign * offset
                if 0 <= cand < n and seq[cand] == "W":
                    w41_pos = cand
                    w41_is_w = True
                    w41_confident = True
                    break
            if w41_pos is not None:
                break
        # If W is mutated, locate via motif within the expected window only
        if w41_pos is None:
            lo = max(0, c23 + 8)
            hi = min(n, c23 + 22)
            window = seq[lo:hi]
            m = FR2_W_FALLBACK.search(window)
            if m:
                w41_pos = lo + m.start()
                w41_is_w = (seq[w41_pos] == "W")

    # Last resort: whole-sequence strict motif (low confidence)
    if w41_pos is None:
        m = FR2_W_STRICT.search(seq)
        if m:
            w41_pos = m.start()
            w41_is_w = True
            w41_confident = True
    if w41_pos is None:
        w41_pos = min(36, int(n * 0.33))
        w41_is_w = (seq[w41_pos] == "W") if w41_pos < n else False

    result["W41_pos"] = w41_pos
    result["W41_is_W"] = w41_is_w
    result["w41_confident"] = w41_confident

    # Back-fill C23 from W41 when C23 was not found (mutated Cys): W41 ≈ C23+14.
    if result["C23_pos"] is None and w41_pos is not None:
        back = w41_pos - 14
        if 0 <= back < n:
            result["C23_pos"] = back
            result["C23_is_C"] = (seq[back] == "C")

    # ================================================================
    # Step 2: J-motif detection (FR4) — 3-tier cascade
    # ================================================================
    j_pos = None
    j_is_fw = False
    j_is_g = False
    j_confident = False

    # Tier 1: Strict [FW]-G-X-G
    j_matches = list(FR4_J_STRICT.finditer(seq))
    if j_matches:
        last = j_matches[-1]
        if last.start() > n - 20:  # must be near C-terminus
            j_pos = last.start()
            j_is_fw = True
            j_is_g = True
            j_confident = True

    # Tier 2: Permissive [FWLIMAVY]-G-X-G (allow F118 substitutions)
    if j_pos is None:
        j_matches = list(FR4_J_FALLBACK.finditer(seq))
        if j_matches:
            last = j_matches[-1]
            if last.start() > n - 20:
                j_pos = last.start()
                j_is_fw = (seq[j_pos] in ("F", "W"))
                j_is_g = (seq[j_pos + 1] == "G")
                j_confident = True

    # Tier 3: G-X-G pattern near C-terminus (low confidence)
    if j_pos is None:
        tail = seq[-20:]  # last 20 residues
        gxg_matches = list(FR4_J_GXG.finditer(tail))
        if gxg_matches:
            # Take the last GXG match
            last_gxg = gxg_matches[-1]
            j_pos = n - 20 + last_gxg.start() - 1  # F/W would be one before the GXG
            if j_pos >= 0:
                j_is_fw = (seq[j_pos] in ("F", "W"))
                j_is_g = (seq[j_pos + 1] == "G")

    # Last resort: estimate from C-terminus (low confidence)
    if j_pos is None:
        j_pos = max(0, n - 12)
        j_is_fw = (seq[j_pos] in ("F", "W")) if j_pos < n else False
        j_is_g = (seq[j_pos + 1] == "G") if j_pos + 1 < n else False

    if j_pos is not None:
        result["F118_pos"] = j_pos
        result["G119_pos"] = j_pos + 1
        result["F118_is_FW"] = j_is_fw
        result["G119_is_G"] = j_is_g
        result["j_confident"] = j_confident

    # (C23 was already located independently in Step 0 — no re-detection here.)

    # ================================================================
    # Step 4: C104 detection — 2nd-CYS in strand F (FR3)
    # ================================================================
    if result["F118_pos"] is not None:
        # C104→F118 distance for VH ≈5 (project handles VH only).
        backoff = 5
        c104_expected = result["F118_pos"] - backoff
        found_c = None
        for offset in range(0, 13):
            for sign in [1, -1]:
                cand = c104_expected + sign * offset
                c23 = result.get("C23_pos") or 0
                if c23 < cand < result["F118_pos"] and seq[cand] == "C":
                    found_c = cand
                    break
            if found_c is not None:
                break
        if found_c is not None:
            result["C104_pos"] = found_c
            result["C104_is_C"] = True
        elif c23 < c104_expected < result["F118_pos"]:
            result["C104_pos"] = c104_expected
            result["C104_is_C"] = False

    # Fallback: last C before J-motif — always prefer actual C over estimate
    if result["F118_pos"] is not None:
        all_c = [i for i, aa in enumerate(seq) if aa == "C"]
        c23 = result.get("C23_pos")
        c104_candidates = [p for p in all_c if p < result["F118_pos"] and p != c23]
        if c104_candidates and not result["C104_is_C"]:
            result["C104_pos"] = c104_candidates[-1]
            result["C104_is_C"] = True

    # Last resort
    if result["C104_pos"] is None:
        result["C104_pos"] = min(n - 10, n - 1)
        result["C104_is_C"] = (seq[result["C104_pos"]] == "C") if result["C104_pos"] < n else False

    return result


# ============================================================
# Step A: Quality Check
# ============================================================

def qc_check(seq, seq_id):
    """
    Perform quality check on a single antibody sequence.

    Returns: {
        "seq_id": str,
        "nonstd_found": [(pos, char, replacement), ...],
        "conserved_issues": [(pos, expected, found, severity), ...],
        "length_issue": None or (actual, issue_type),
        "quality_warnings": [(pos, motif, warning_type), ...],
        "qc_passed": bool,
        "qc_score": float (0-1),
        "conserved_positions": dict,
    }
    """
    report = {
        "seq_id": seq_id,
        "nonstd_found": [],
        "conserved_issues": [],
        "length_issue": None,
        "quality_warnings": [],
        "qc_passed": True,
        "qc_score": 1.0,
        "conserved_positions": {},
    }

    # --- 1. Non-standard amino acid check ---
    for i, aa in enumerate(seq):
        if aa not in STANDARD_AA:
            replacement = NONSTD_MAP.get(aa)
            report["nonstd_found"].append((i, aa, replacement))
            report["qc_passed"] = False

    # --- 2. Length check ---
    length = len(seq)
    if length < 90:
        report["length_issue"] = (length, "too_short")
        report["qc_passed"] = False
    elif length > 140:
        report["length_issue"] = (length, "too_long")
        report["qc_passed"] = False

    # --- 3. Conserved position check ---
    positions = find_conserved_positions(seq)
    report["conserved_positions"] = positions

    # Tier 1: Structurally essential positions
    if positions["C23_pos"] is not None:
        aa = seq[positions["C23_pos"]]
        if aa != "C":
            report["conserved_issues"].append(
                (positions["C23_pos"], "C", aa, "critical")
            )
            report["qc_passed"] = False
    else:
        report["conserved_issues"].append((-1, "C", "?", "critical"))
        report["qc_passed"] = False

    if positions["C104_pos"] is not None:
        aa = seq[positions["C104_pos"]]
        if aa != "C":
            report["conserved_issues"].append(
                (positions["C104_pos"], "C", aa, "critical")
            )
            report["qc_passed"] = False
    else:
        report["conserved_issues"].append((-1, "C", "?", "critical"))
        report["qc_passed"] = False

    # Tier 2: Highly conserved positions
    if positions["W41_pos"] is not None:
        aa = seq[positions["W41_pos"]]
        if aa != "W":
            report["conserved_issues"].append(
                (positions["W41_pos"], "W", aa, "high")
            )
            report["qc_passed"] = False
    else:
        report["conserved_issues"].append((-1, "W", "?", "high"))

    if positions["F118_pos"] is not None:
        aa = seq[positions["F118_pos"]]
        if aa not in ("F", "W"):
            report["conserved_issues"].append(
                (positions["F118_pos"], "F/W", aa, "high")
            )
            report["qc_passed"] = False
    else:
        report["conserved_issues"].append((-1, "F/W", "?", "high"))

    if positions["G119_pos"] is not None:
        aa = seq[positions["G119_pos"]]
        if aa != "G":
            report["conserved_issues"].append(
                (positions["G119_pos"], "G", aa, "high")
            )
            report["qc_passed"] = False

    # --- 4. Quality motifs check ---
    # N-glycosylation: N-X-S/T where X≠P
    nglyc_pattern = re.compile(r"N[^P][ST]")
    for m in nglyc_pattern.finditer(seq):
        pos = m.start()
        # Check if this is likely in a CDR (approximate)
        if positions["W41_pos"] and pos > positions["W41_pos"]:
            report["quality_warnings"].append(
                (pos, m.group(), "N-glycosylation motif")
            )

    # Deamidation: NG motif
    deam_pattern = re.compile(r"NG")
    for m in deam_pattern.finditer(seq):
        pos = m.start()
        report["quality_warnings"].append(
            (pos, m.group(), "deamidation hot spot (NG)")
        )

    # Isomerization: DG, DS, DD motifs
    iso_pattern = re.compile(r"D[GSD]")
    for m in iso_pattern.finditer(seq):
        pos = m.start()
        report["quality_warnings"].append(
            (pos, m.group(), "isomerization hot spot")
        )

    # Oxidation: exposed M (all M flagged as potential risk)
    for i, aa in enumerate(seq):
        if aa == "M":
            report["quality_warnings"].append(
                (i, "M", "potential methionine oxidation site")
            )

    # Free cysteine: odd number of C residues
    c_count = seq.count("C")
    if c_count % 2 != 0:
        report["quality_warnings"].append((-1, f"C_count={c_count}", "unpaired cysteine (odd count)"))

    # --- 5. Calculate QC score ---
    score = 1.0
    # Deduct for non-standard AA
    score -= len(report["nonstd_found"]) * 0.10
    # Deduct for conserved issues
    for _, _, _, severity in report["conserved_issues"]:
        if severity == "critical":
            score -= 0.25
        elif severity == "high":
            score -= 0.12
    # Deduct for length
    if report["length_issue"]:
        score -= 0.30
    report["qc_score"] = max(0.0, score)

    return report


# ============================================================
# Step B: Correction
# ============================================================

def correct_sequence(seq, qc_report):
    """
    Apply corrections to a sequence based on QC findings.

    Corrections applied:
    1. Non-standard AA → standard AA replacement
    2. C23/C104 → C (if identified and mutated)
    3. W41 → W (if identified and mutated)
    4. F118 → F (VH) or W (VL) (if identified and mutated)

    NOT corrected (marked for audit):
    - G119 mutation (structurally unpredictable)

    Returns: {
        "corrected_seq": str,
        "corrections_made": [(pos, old, new, reason), ...],
        "uncorrectable": [(pos, old, expected, reason), ...],
    }
    """
    seq_list = list(seq)
    corrections = []
    uncorrectable = []
    positions = qc_report.get("conserved_positions", {})

    # --- 1. Fix non-standard amino acids ---
    for pos, old_aa, replacement in qc_report.get("nonstd_found", []):
        if replacement is not None:
            seq_list[pos] = replacement
            corrections.append((pos, old_aa, replacement,
                              f"Non-standard AA '{old_aa}' → '{replacement}'"))
        elif old_aa in ("*", "-"):
            seq_list[pos] = ""  # Will be removed
            corrections.append((pos, old_aa, "REMOVED",
                              f"Stop codon / gap '{old_aa}' removed"))

    # Remove empty positions (stop codons / gaps)
    corrected_seq = "".join(seq_list)

    # Recalculate positions after non-std AA removal
    seq_list = list(corrected_seq)

    # Re-identify conserved positions on cleaned sequence
    # (positions may shift if stop codons/gaps were removed)
    if any(c[1] in ("*", "-") for c in corrections):
        positions = find_conserved_positions(corrected_seq)

    # --- 2. Fix critical conserved positions ---
    # C23 → C
    if positions.get("C23_pos") is not None:
        pos = positions["C23_pos"]
        if pos < len(corrected_seq) and corrected_seq[pos] != "C":
            old = corrected_seq[pos]
            seq_list[pos] = "C"
            corrections.append((pos, old, "C",
                f"C23 (1st-CYS, strand B, intra-domain disulfide bridge) restored to Cys; "
                f"mutation '{old}' would break Ig fold (Lefranc et al., 2003; Chothia & Lesk, 1987)"))

    # C104 → C
    if positions.get("C104_pos") is not None:
        pos = positions["C104_pos"]
        if pos < len(corrected_seq) and corrected_seq[pos] != "C":
            old = corrected_seq[pos]
            seq_list[pos] = "C"
            corrections.append((pos, old, "C",
                f"C104 (2nd-CYS, strand F, intra-domain disulfide bridge) restored to Cys; "
                f"mutation '{old}' would break Ig fold (Lefranc et al., 2003)"))

    # W41 → W
    if positions.get("W41_pos") is not None:
        pos = positions["W41_pos"]
        if pos < len(corrected_seq) and corrected_seq[pos] != "W":
            old = corrected_seq[pos]
            seq_list[pos] = "W"
            corrections.append((pos, old, "W",
                f"W41 (FR2 hydrophobic core) restored to Trp; "
                f"mutation '{old}' disrupts VH/VL interface (Ewert et al., 2003; Honegger & Plückthun, 2001)"))

    # F118: VH's F118 is W (JH gene). Only correct if NOT aromatic (F or W);
    # if already aromatic, leave it alone. (Project handles VH only.)
    if positions.get("F118_pos") is not None:
        pos = positions["F118_pos"]
        if pos < len(corrected_seq):
            current = corrected_seq[pos]
            if current not in ("F", "W"):
                seq_list[pos] = "W"
                corrections.append((pos, current, "W",
                    f"F118 (FR4 J-motif aromatic) restored to Trp for VH (JH gene); "
                    f"(Lefranc, 2014)"))

    # G119 → NOT corrected (structural consequences unpredictable)
    if positions.get("G119_pos") is not None:
        pos = positions["G119_pos"]
        if pos < len(corrected_seq) and corrected_seq[pos] != "G":
            old = corrected_seq[pos]
            uncorrectable.append((pos, old, "G",
                f"G119 (FR4 J-motif Gly) cannot be reliably corrected — "
                f"tight turn requires Gly; mutation '{old}' has unpredictable structural consequences "
                f"(Lefranc et al., 2003)"))

    corrected_seq = "".join(seq_list)
    return {
        "corrected_seq": corrected_seq,
        "corrections_made": corrections,
        "uncorrectable": uncorrectable,
    }


# ============================================================
# Step C: Audit
# ============================================================

def audit_sequence(original_seq, corrected_seq, qc_report, correction_result):
    """
    Audit the corrected sequence and determine PASS/FAIL status (two outcomes).

    Scoring rubric (see docs/02_conservation_rules.md Section 2.3):
    - No non-standard AA: 15%
    - C23 intact: 25%
    - C104 intact: 25%
    - W41 intact: 10%
    - F/W118 intact: 10%
    - Length normal: 10%
    - Other quality: 5%

    Status: hard failure conditions → FAIL; otherwise PASS if score ≥ 0.94
    (data-driven cutoff from validate_v3.py), else FAIL.
    """

    seq = corrected_seq
    score = 0.0
    details = {}

    positions = find_conserved_positions(seq)
    # Re-detect on the ORIGINAL sequence to check for pre-correction structural
    # damage that should NOT be auto-repaired (e.g. both disulfide cysteines
    # missing, or no antibody signature motifs at all).
    orig_positions = find_conserved_positions(original_seq)
    double_c_missing = (not orig_positions["C23_is_C"]) and (not orig_positions["C104_is_C"])

    # Criterion 1: No non-standard AA (15%)
    nonstd_remaining = sum(1 for aa in seq if aa not in STANDARD_AA)
    if nonstd_remaining == 0:
        score += 0.15
        details["nonstd_clean"] = 1.0
    else:
        details["nonstd_clean"] = 0.0

    # Criterion 2: C23 intact (25%)
    # Full weight whether native or corrected — after correction the structure
    # is compliant, so it should score the same as a native C. The
    # "corrected" vs "intact" distinction is recorded in details only.
    has_c23 = False
    if positions["C23_pos"] is not None:
        c23_was_corrected = any(
            c[0] == positions["C23_pos"] and c[2] == "C"
            for c in correction_result["corrections_made"]
        )
        if seq[positions["C23_pos"]] == "C":
            has_c23 = True
            score += 0.25
            details["C23"] = "corrected" if c23_was_corrected else "intact"
    if not has_c23:
        details["C23"] = "missing"

    # Criterion 3: C104 intact (25%)
    has_c104 = False
    if positions["C104_pos"] is not None:
        c104_was_corrected = any(
            c[0] == positions["C104_pos"] and c[2] == "C"
            for c in correction_result["corrections_made"]
        )
        if seq[positions["C104_pos"]] == "C":
            has_c104 = True
            score += 0.25
            details["C104"] = "corrected" if c104_was_corrected else "intact"
    if not has_c104:
        details["C104"] = "missing"

    # Criterion 4: W41 intact (10%)
    if positions["W41_pos"] is not None and seq[positions["W41_pos"]] == "W":
        w41_was_corrected = any(
            c[0] == positions["W41_pos"] and c[2] == "W"
            for c in correction_result["corrections_made"]
        )
        score += 0.10
        details["W41"] = "corrected" if w41_was_corrected else "intact"
    else:
        details["W41"] = "missing_or_mutated"

    # Criterion 5: F/W118 intact (10%)
    if positions["F118_pos"] is not None:
        aa_at_118 = seq[positions["F118_pos"]]
        if aa_at_118 in ("F", "W"):
            f118_was_corrected = any(
                c[0] == positions["F118_pos"]
                for c in correction_result["corrections_made"]
            )
            score += 0.10
            details["F118"] = "corrected" if f118_was_corrected else "intact"
        else:
            details["F118"] = "mutated"
    else:
        details["F118"] = "not_found"

    # Criterion 6: Length normal (10%)
    length = len(seq)
    if 95 <= length <= 135:
        score += 0.10
        details["length"] = "normal"
    elif 90 <= length < 95 or 135 < length <= 140:
        score += 0.10 * 0.5
        details["length"] = "borderline"
    else:
        details["length"] = "abnormal"

    # Criterion 7: Other quality (5%)
    # Deduct for critical uncorrectable issues
    uncorrectable_count = len(correction_result.get("uncorrectable", []))
    g119_fail = any("G119" in u[3] for u in correction_result.get("uncorrectable", []))
    if g119_fail:
        details["quality"] = "G119_mutated_uncorrectable"
    elif uncorrectable_count > 0:
        score += 0.05 * 0.5
        details["quality"] = "warnings"
    else:
        score += 0.05
        details["quality"] = "clean"

    score = round(score, 4)

    # Determine status — two outcomes only: PASS (通过) / FAIL (不通过).
    # Hard failure conditions → FAIL; otherwise PASS if score ≥ 0.94.
    status = "PASS"
    if g119_fail:
        # G119 mutation is structurally catastrophic → FAIL
        status = "FAIL"
    elif length < 90 or length > 140:
        # Not a V domain → FAIL (length anomaly)
        status = "FAIL"
    elif not (orig_positions.get("w41_confident") or orig_positions.get("j_confident")):
        # No antibody signature motifs found → not antibody-like → FAIL
        details["quality"] = "not_antibody_like"
        status = "FAIL"
    elif double_c_missing:
        # Both intra-domain disulfide cysteines missing → FAIL
        details["quality"] = "double_disulfide_missing"
        status = "FAIL"
    elif not has_c23 or not has_c104:
        # A single conserved Cys missing → FAIL
        status = "FAIL"
    elif score >= 0.94:
        status = "PASS"
    else:
        status = "FAIL"

    return {
        "final_score": score,
        "status": status,
        "score_details": details,
        "is_pass": status == "PASS",
    }


# ============================================================
# FASTA I/O
# ============================================================

def read_fasta(filepath):
    """Read FASTA file, return list of (header, sequence) tuples."""
    entries = []
    current_header = None
    current_seq = []

    with open(filepath, "r", encoding="utf-8") as f:
        for line in f:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_header:
                    entries.append((current_header, "".join(current_seq)))
                current_header = line[1:]  # Remove '>'
                current_seq = []
            else:
                current_seq.append(line)
        if current_header:
            entries.append((current_header, "".join(current_seq)))

    return entries


def parse_fasta_header(header):
    """Parse FASTA header fields into dict."""
    info = {"raw": header}
    parts = header.split("|")
    for part in parts:
        if "=" in part:
            key, val = part.split("=", 1)
            info[key.strip()] = val.strip()
    return info


def write_fasta(filepath, header, sequence):
    """Write single sequence as FASTA."""
    os.makedirs(os.path.dirname(filepath), exist_ok=True)
    with open(filepath, "w", encoding="utf-8") as f:
        f.write(f">{header}\n")
        for i in range(0, len(sequence), 80):
            f.write(sequence[i : i + 80] + "\n")


# ============================================================
# Per-sequence processing
# ============================================================

def process_sequence(header, seq, seq_id):
    """
    Process a single sequence through QC → Correction → Audit pipeline.

    Returns: {
        "seq_id": str,
        "original_seq": str,
        "corrected_seq": str,
        "status": "PASS"|"FAIL",
        "qc_report": dict,
        "correction_result": dict,
        "audit_result": dict,
    }
    """
    result = {
        "seq_id": seq_id,
        "original_seq": seq,
        "corrected_seq": seq,
        "status": "PASS",
    }

    # Step A: Quality Check
    qc_report = qc_check(seq, seq_id)
    result["qc_report"] = qc_report

    # Step B: Correction
    correction_result = correct_sequence(seq, qc_report)
    result["correction_result"] = correction_result
    result["corrected_seq"] = correction_result["corrected_seq"]

    # Step C: Audit
    audit_result = audit_sequence(seq, correction_result["corrected_seq"],
                                   qc_report, correction_result)
    result["audit_result"] = audit_result
    result["status"] = audit_result["status"]

    return result


# ============================================================
# Batch processing
# ============================================================

def process_batch(input_fasta, output_dir, batch_id):
    """Process all sequences in a batch."""
    entries = read_fasta(input_fasta)
    print(f"Processing {len(entries)} sequences from {input_fasta}\n")

    results = []
    stats = {"PASS": 0, "FAIL": 0}

    for header, seq in entries:
        info = parse_fasta_header(header)
        # Try to extract SeqID
        seq_id = info.get("Global_Seq_ID", "")
        if not seq_id:
            for k in info:
                if "SEQ" in k or "seq" in k.lower():
                    seq_id = info[k]
                    break
        if not seq_id:
            # Derive from first token of header
            seq_id = header.split("|")[0].strip()

        print(f"  [{seq_id}] ", end="", flush=True)

        result = process_sequence(header, seq, seq_id)
        results.append(result)

        status = result["status"]
        stats[status] = stats.get(status, 0) + 1
        score = result["audit_result"]["final_score"]
        corrections = len(result["correction_result"]["corrections_made"])
        print(f"{status}  score={score:.2f}  corrections={corrections}")

        # --- Write per-sequence output files ---
        seq_dir = os.path.join(output_dir, seq_id)
        raw_fa = os.path.join(seq_dir, "raw.fa")
        processed_dir = os.path.join(seq_dir, "processed")
        logs_dir = os.path.join(seq_dir, "logs")
        reports_dir = os.path.join(seq_dir, "reports")

        os.makedirs(processed_dir, exist_ok=True)
        os.makedirs(logs_dir, exist_ok=True)
        os.makedirs(reports_dir, exist_ok=True)

        # raw.fa
        write_fasta(raw_fa, header, seq)

        # processed/v01.fa or final.fa
        if result["correction_result"]["corrections_made"]:
            # Has corrections: write v01.fa
            v01_fa = os.path.join(processed_dir, "v01.fa")
            write_fasta(v01_fa, header, result["corrected_seq"])
            # If PASS, also write as final.fa
            if status == "PASS":
                final_fa = os.path.join(processed_dir, "final.fa")
                write_fasta(final_fa, header, result["corrected_seq"])
        elif status == "PASS":
            # No corrections needed: directly final
            final_fa = os.path.join(processed_dir, "final.fa")
            write_fasta(final_fa, header, seq)

        # correction_log.json
        corr_log = {
            "seq_id": seq_id,
            "timestamp": datetime.now().isoformat(),
            "corrections": [
                {"position": pos, "original": old, "corrected": new, "reason": reason}
                for pos, old, new, reason in result["correction_result"]["corrections_made"]
            ],
            "uncorrectable": [
                {"position": pos, "original": old, "expected": exp, "reason": reason}
                for pos, old, exp, reason in result["correction_result"]["uncorrectable"]
            ],
        }
        with open(os.path.join(logs_dir, "correction_log.json"), "w", encoding="utf-8") as f:
            json.dump(corr_log, f, indent=2)

        # audit_log.json
        audit_log = {
            "seq_id": seq_id,
            "timestamp": datetime.now().isoformat(),
            "status": status,
            "final_score": result["audit_result"]["final_score"],
            "score_details": result["audit_result"]["score_details"],
            "qc_score": result["qc_report"]["qc_score"],
            "qc_passed": result["qc_report"]["qc_passed"],
            "nonstd_found": len(result["qc_report"]["nonstd_found"]),
            "conserved_issues": len(result["qc_report"]["conserved_issues"]),
            "quality_warnings": len(result["qc_report"]["quality_warnings"]),
        }
        with open(os.path.join(logs_dir, "audit_log.json"), "w", encoding="utf-8") as f:
            json.dump(audit_log, f, indent=2)

        # structure_report.md
        qr = result["qc_report"]
        ar = result["audit_result"]
        cr = result["correction_result"]

        report_lines = [
            f"# Structure Compliance Report: {seq_id}",
            f"",
            f"**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M')}",
            f"**Status:** {status}",
            f"**Final Score:** {ar['final_score']:.2f} / 1.00",
            f"",
            f"## 1. Quality Check Summary",
            f"",
            f"- QC Score: {qr['qc_score']:.2f}",
            f"- Non-standard amino acids found: {len(qr['nonstd_found'])}",
            f"- Conserved position issues: {len(qr['conserved_issues'])}",
            f"- Quality warnings: {len(qr['quality_warnings'])}",
            f"",
        ]

        if qr["nonstd_found"]:
            report_lines.append("### Non-standard Amino Acids:")
            for pos, aa, repl in qr["nonstd_found"][:10]:
                report_lines.append(f"- Position {pos}: '{aa}' → '{repl}'")
            report_lines.append("")

        if qr["conserved_issues"]:
            report_lines.append("### Conserved Position Issues:")
            for pos, expected, found, severity in qr["conserved_issues"]:
                report_lines.append(
                    f"- Position {pos}: expected '{expected}', found '{found}' [{severity}]"
                )
            report_lines.append("")

        if cr["corrections_made"]:
            report_lines.append("## 2. Corrections Applied")
            for pos, old, new, reason in cr["corrections_made"]:
                report_lines.append(f"- Pos {pos}: {old}→{new} — {reason}")
            report_lines.append("")

        if cr["uncorrectable"]:
            report_lines.append("## 3. Uncorrectable Issues")
            for pos, old, exp, reason in cr["uncorrectable"]:
                report_lines.append(f"- Pos {pos}: {old} (expected {exp}) — {reason}")
            report_lines.append("")

        report_lines.extend([
            f"## Scoring Details",
            f"",
            f"| Criterion | Result |",
            f"|-----------|--------|",
        ])
        for k, v in ar["score_details"].items():
            report_lines.append(f"| {k} | {v} |")
        report_lines.append("")

        if qr["quality_warnings"]:
            report_lines.append("## Quality Warnings")
            for pos, motif, wtype in qr["quality_warnings"]:
                report_lines.append(f"- Pos {pos}: {motif} ({wtype})")
            report_lines.append("")

        report_lines.extend([
            f"## References",
            f"- Lefranc et al. (2003) Dev. Comp. Immunol. 27(1):55-77",
            f"- Chothia & Lesk (1987) JMB 196(4):901-917",
            f"- Ewert et al. (2003) JMB 325(3):531-553",
            f"- Honegger & Plückthun (2001) JMB 309(3):657-670",
            f"- Lefranc (2014) Front. Immunol. 5:22",
        ])

        with open(os.path.join(reports_dir, "structure_report.md"), "w", encoding="utf-8") as f:
            f.write("\n".join(report_lines))

    # --- Write manifest.csv ---
    manifest_path = os.path.join(output_dir, "manifest.csv")
    with open(manifest_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["Global_Seq_ID", "Sequence", "Status"])
        for r in results:
            writer.writerow([r["seq_id"], r["corrected_seq"], r["status"]])

    # --- Write batch summary ---
    summary_dir = os.path.join(output_dir, "batch_correction_summary")
    os.makedirs(summary_dir, exist_ok=True)

    summary_lines = [
        f"# Batch Correction Summary: {batch_id}",
        f"",
        f"**Date:** {datetime.now().strftime('%Y-%m-%d %H:%M')}",
        f"**Total sequences processed:** {len(results)}",
        f"",
        f"## Results",
        f"",
        f"| Status | Count | Percentage |",
        f"|--------|-------|------------|",
    ]
    for s in ["PASS", "FAIL"]:
        pct = stats.get(s, 0) / len(results) * 100
        summary_lines.append(f"| {s} | {stats.get(s, 0)} | {pct:.1f}% |")
    summary_lines.append("")

    summary_lines.extend([
        f"## Score Distribution",
        f"",
        f"| Range | Count |",
        f"|-------|-------|",
    ])
    ranges = [("0.9-1.0", 0.9, 1.0), ("0.8-0.9", 0.8, 0.9),
              ("0.6-0.8", 0.6, 0.8), ("0.4-0.6", 0.4, 0.6),
              ("<0.4", 0.0, 0.4)]
    for label, lo, hi in ranges:
        count = sum(1 for r in results if lo <= r["audit_result"]["final_score"] < hi
                    or (label == "<0.4" and r["audit_result"]["final_score"] < 0.4))
        summary_lines.append(f"| {label} | {count} |")
    summary_lines.append("")

    summary_lines.extend([
        f"## Per-Sequence Detail",
        f"",
        f"| Seq_ID | Status | Score | Corrections | Warnings |",
        f"|--------|--------|-------|-------------|----------|",
    ])
    for r in results:
        n_corr = len(r["correction_result"]["corrections_made"])
        n_warn = len(r["qc_report"]["quality_warnings"])
        summary_lines.append(
            f"| {r['seq_id']} | {r['status']} | {r['audit_result']['final_score']:.2f} | {n_corr} | {n_warn} |"
        )

    with open(os.path.join(summary_dir, "batch_correction_summary.md"), "w", encoding="utf-8") as f:
        f.write("\n".join(summary_lines))

    # --- Print summary ---
    print(f"\n{'='*60}")
    print(f"Batch Processing Complete: {batch_id}")
    print(f"{'='*60}")
    print(f"  Total:   {len(results)}")
    print(f"  PASS:    {stats.get('PASS', 0)}")
    print(f"  FAIL:    {stats.get('FAIL', 0)}")
    print(f"\n  Manifest: {manifest_path}")
    print(f"  Summary:  {summary_dir}/batch_correction_summary.md")

    return results


# ============================================================
# Main
# ============================================================

def batch_id_from_path(filepath):
    """Derive a BatchID from an input FASTA filename.

    'Abseqs_BATCH_HER2_20260809_01.fa' -> 'BATCH_HER2_20260809_01'
    'D2_candidates_50.fasta'           -> 'D2_candidates_50'
    """
    name = os.path.basename(filepath)
    name = re.sub(r"\.(fa|fasta|fas|fna)$", "", name, flags=re.IGNORECASE)
    return re.sub(r"^Abseqs_", "", name)


def main():
    parser = argparse.ArgumentParser(
        # ASCII only: the Windows console codepage cannot render '→' in --help.
        description="Module 2: antibody VH sequence QC -> correction -> audit")
    parser.add_argument("--batch-id", "-b",
                        help="Batch ID, e.g. BATCH_HER2_20260809_01. "
                             "Inferred from the input filename if omitted.")
    parser.add_argument("--input", "-i",
                        help="Input FASTA. Defaults to "
                             "data/01_generated/<BatchID>/Abseqs_<BatchID>.fa")
    parser.add_argument("--root", default=ROOT_DIR,
                        help="Project root (default: the directory above scripts/)")
    args = parser.parse_args()

    if not args.input and not args.batch_id:
        parser.error("give at least one of -b/--batch-id or -i/--input")

    batch_id = args.batch_id or batch_id_from_path(args.input)
    input_fa = args.input or os.path.join(
        args.root, "data", "01_generated", batch_id, f"Abseqs_{batch_id}.fa")
    output_dir = os.path.join(args.root, "data", "02_correction", batch_id)

    if not os.path.exists(input_fa):
        print(f"Error: input FASTA not found: {input_fa}")
        print(f"Usage: python scripts/pipeline.py -b {batch_id} [-i <input.fa>]")
        sys.exit(1)

    os.makedirs(output_dir, exist_ok=True)
    process_batch(input_fa, output_dir, batch_id)


if __name__ == "__main__":
    main()
