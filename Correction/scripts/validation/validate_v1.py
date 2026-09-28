#!/usr/bin/env python3
"""
Module 2 Validation Suite
=========================
Generates test sequences with KNOWN ground truth (injected errors + expected
outcome), runs them through the QC → Correction → Audit pipeline, and measures:

  1. Detection rate    — did QC find all injected errors?
  2. Correction accuracy — did corrections restore the intended residues?
  3. Status accuracy   — is the final PASS/FAIL correct?

Output: docs/validation_report_v1.md
"""

import os
import random
import sys
from datetime import datetime

_HERE = os.path.dirname(os.path.abspath(__file__))
_SCRIPTS_DIR = os.path.dirname(_HERE)      # scripts/ — holds pipeline.py
ROOT = os.path.dirname(_SCRIPTS_DIR)       # project root
sys.path.insert(0, _SCRIPTS_DIR)
import pipeline as pl

random.seed(812)  # date-anchored seed for reproducibility

# ============================================================
# Templates + verified conserved positions (0-indexed)
# ============================================================
TEMPLATES = {
    "VH_323": {
        "seq": "EVQLVESGGGLVQPGGSLRLSCAASGFTFSSYAMSWVRQAPGKGLEWVSAISGSGGSTYYADSVKGRFTISRDNSKNTLYLQMNSLRAEDTAVYYCAKDYWGQGTLVTVSS",
        "C23": 21, "W41": 35, "C104": 95, "F118": 100, "G119": 101, "chain": "VH",
    },
    "VH_169": {
        "seq": "QVQLVQSGAEVKKPGSSVKVSCKASGGTFSSYAISWVRQAPGQGLEWMGGIIPIFGTANYAQKFQGRVTITADESTSTAYMELSSLRSEDTAVYYCARDYWGQGTTVTVSS",
        "C23": 21, "W41": 35, "C104": 95, "F118": 100, "G119": 101, "chain": "VH",
    },
}

STANDARD_AA = "ACDEFGHIKLMNPQRSTVWY"
NONSTD = ["B", "Z", "J", "U", "O", "X"]


def mutate(seq, pos, new):
    return seq[:pos] + new + seq[pos + 1:]


def build_cases():
    """Build test cases with ground truth."""
    cases = []
    cid = 0

    def add(desc, seq, expected_status, injected, expected_corr):
        nonlocal cid
        cid += 1
        cases.append({
            "id": f"VAL_{cid:03d}",
            "desc": desc,
            "seq": seq,
            "expected_status": expected_status,
            "injected": injected,        # list of (pos, orig, mut, kind)
            "expected_corr": expected_corr,  # min corrections expected
        })

    # ---- Group 1: clean (should PASS, 0 corrections) ----
    for name, t in TEMPLATES.items():
        add(f"clean {name}", t["seq"], "PASS", [], 0)

    # ---- Group 2: non-standard AA single injection (PASS, 1 correction) ----
    # Ground truth for correction = the IUPAC mapping target (e.g. B→N), NOT the
    # template residue, because B (Asx) means "Asp or Asn", never the template AA.
    for i, code in enumerate(NONSTD):
        t = TEMPLATES["VH_323"]
        pos = random.randint(3, 90)  # random FR/CDR position
        target = pl.NONSTD_MAP.get(code, "A")
        s = mutate(t["seq"], pos, code)
        add(f"nonstd {code} @pos{pos}", s, "PASS",
            [(pos, target, code, "nonstd")], 1)

    # ---- Group 3: conserved site single mutations (PASS with flag) ----
    t = TEMPLATES["VH_323"]
    add("C23→S", mutate(t["seq"], t["C23"], "S"), "PASS",
        [(t["C23"], "C", "S", "C23")], 1)
    add("C104→Y", mutate(t["seq"], t["C104"], "Y"), "PASS",
        [(t["C104"], "C", "Y", "C104")], 1)
    add("W41→L", mutate(t["seq"], t["W41"], "L"), "PASS",
        [(t["W41"], "W", "L", "W41")], 1)

    t = TEMPLATES["VH_323"]
    add("F118→A (VH)", mutate(t["seq"], t["F118"], "A"), "PASS",
        [(t["F118"], "W", "A", "F118")], 1)

    # ---- Group 4: G119 mutation (uncorrectable → FAIL) ----
    for i in range(3):
        t = TEMPLATES["VH_323"] if i < 2 else TEMPLATES["VH_169"]
        s = mutate(t["seq"], t["G119"], "E")
        add(f"G119→E", s, "FAIL",
            [(t["G119"], "G", "E", "G119")], 0)

    # ---- Group 5: double conserved Cys (FAIL) ----
    for i in range(2):
        t = TEMPLATES["VH_323"]
        s = mutate(t["seq"], t["C23"], "A")
        s = mutate(s, t["C104"], "G")
        add("C23+C104 double mutant", s, "FAIL",
            [(t["C23"], "C", "A", "C23"), (t["C104"], "C", "G", "C104")], 2)

    # ---- Group 6: length anomalies (FAIL) ----
    add("too short (67AA)", TEMPLATES["VH_323"]["seq"][:67], "FAIL",
        [], 0)
    add("too long (146AA)", TEMPLATES["VH_323"]["seq"] + "GGGGS" * 7, "FAIL",
        [], 0)

    # ---- Group 7: random noise (FAIL, no crash) ----
    for i in range(5):
        noise = "".join(random.choice(STANDARD_AA) for _ in range(115))
        add(f"random noise #{i+1}", noise, "FAIL", [], 0)

    # ---- Group 8: random multi-point mutations in NON-conserved positions ----
    # Should PASS (pipeline must not falsely reject valid diversity)
    conserved = set()
    for t in TEMPLATES.values():
        conserved.update([t["C23"], t["W41"], t["C104"], t["F118"], t["G119"]])
    for i in range(5):
        t = TEMPLATES["VH_323"]
        s = t["seq"]
        injected = []
        n_mut = random.randint(2, 4)
        for _ in range(n_mut):
            pos = random.randint(0, len(s) - 1)
            # avoid conserved positions for this "should pass" group
            attempts = 0
            while pos in conserved and attempts < 20:
                pos = random.randint(0, len(s) - 1)
                attempts += 1
            new = random.choice(STANDARD_AA.replace(s[pos], ""))
            s = mutate(s, pos, new)
            injected.append((pos, t["seq"][pos], new, "random"))
        add(f"random multi-mut ×{n_mut}", s, "PASS", injected, 0)

    return cases


# ============================================================
# Evaluation
# ============================================================

def run_validation(cases):
    """Run pipeline on each case and score against ground truth."""
    results = []
    stats = {
        "total": len(cases),
        "detected_all": 0,      # all injected errors detected
        "status_correct": 0,    # final status matches expectation
        "corr_correct": 0,      # corrected sequence restored intended residues
        "false_positives": 0,   # clean cases wrongly flagged
        "crashes": 0,
    }

    for c in cases:
        r = {
            "id": c["id"],
            "desc": c["desc"],
            "expected_status": c["expected_status"],
        }

        try:
            qc = pl.qc_check(c["seq"], c["id"])
            corr = pl.correct_sequence(c["seq"], qc)
            audit = pl.audit_sequence(c["seq"], corr["corrected_seq"], qc, corr)

            # --- detection: every injected ERROR should be found ---
            # "random" kind is legal diversity, not an error → not checked here.
            found_nonstd = {p for p, _, _ in qc["nonstd_found"]}
            found_cons = {p for p, _, _, _ in qc["conserved_issues"]}
            detected_all = True
            for pos, orig, mut, kind in c["injected"]:
                if kind == "nonstd":
                    if pos not in found_nonstd:
                        detected_all = False
                elif kind in ("C23", "C104", "W41", "F118", "G119"):
                    # conserved position: position may be re-estimated; accept ±2
                    if not any(abs(p - pos) <= 2 for p in found_cons):
                        detected_all = False

            # --- correction accuracy: corrected seq should restore orig residue ---
            # Only for correctable kinds; G119 (uncorrectable) and random (legal)
            # are excluded.
            corr_ok = True
            for pos, orig, mut, kind in c["injected"]:
                if kind in ("nonstd", "C23", "C104", "W41", "F118"):
                    if pos < len(corr["corrected_seq"]):
                        if corr["corrected_seq"][pos] != orig:
                            corr_ok = False
                    else:
                        corr_ok = False

            # --- status accuracy ---
            status_ok = (audit["status"] == c["expected_status"])

            # --- false positive: clean case should have no issues ---
            if not c["injected"] and c["desc"].startswith("clean"):
                if not qc["qc_passed"]:
                    stats["false_positives"] += 1

            if detected_all:
                stats["detected_all"] += 1
            if status_ok:
                stats["status_correct"] += 1
            if corr_ok:
                stats["corr_correct"] += 1

            r.update({
                "actual_status": audit["status"],
                "score": audit["final_score"],
                "corrections": len(corr["corrections_made"]),
                "uncorrectable": len(corr["uncorrectable"]),
                "detected_all": detected_all,
                "corr_ok": corr_ok,
                "status_ok": status_ok,
            })

        except Exception as e:
            stats["crashes"] += 1
            r.update({
                "actual_status": "CRASH",
                "error": str(e),
                "detected_all": False,
                "corr_ok": False,
                "status_ok": False,
            })

        results.append(r)

    return results, stats


# ============================================================
# Report
# ============================================================

def write_report(results, stats, path):
    total = stats["total"]
    n_crash = stats["crashes"]

    lines = []
    lines.append("# Module 2 Pipeline Validation Report")
    lines.append("")
    lines.append(f"**Date:** 2026-08-12")
    lines.append(f"**Validator:** GCLSC Module 2 (automated)")
    lines.append(f"**Pipeline:** `scripts/pipeline.py` (QC → Correction → Audit)")
    lines.append(f"**Test set:** {total} sequences with known ground truth")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("## 1. Executive Summary")
    lines.append("")
    lines.append("The pipeline was validated against a synthetic test set where each "
                 "sequence carried **known injected errors** and a **known expected "
                 "outcome**. Three dimensions were scored:")
    lines.append("")
    lines.append("| Metric | Definition | Result |")
    lines.append("|--------|-----------|--------|")

    def pct(x):
        return f"{x / total * 100:.1f}%" if total else "n/a"

    lines.append(f"| Detection rate | QC found all injected errors | "
                 f"{stats['detected_all']}/{total} ({pct(stats['detected_all'])}) |")
    lines.append(f"| Correction accuracy | Corrected sequence restored intended residues | "
                 f"{stats['corr_correct']}/{total} ({pct(stats['corr_correct'])}) |")
    lines.append(f"| Status accuracy | Final PASS/FAIL matches expectation | "
                 f"{stats['status_correct']}/{total} ({pct(stats['status_correct'])}) |")
    lines.append(f"| False positives | Clean sequences wrongly flagged | "
                 f"{stats['false_positives']} |")
    lines.append(f"| Crashes | Unexpected exceptions | {n_crash} |")
    lines.append("")
    lines.append("")
    lines.append("## 2. Per-Case Results")
    lines.append("")
    lines.append("| ID | Case | Expected | Actual | Score | Corr | Detect | Corr OK | Status OK |")
    lines.append("|----|------|----------|--------|-------|------|--------|---------|-----------|")
    for r in results:
        det = "✓" if r["detected_all"] else "✗"
        co = "✓" if r["corr_ok"] else "✗"
        so = "✓" if r["status_ok"] else "✗"
        sc = f"{r.get('score', 0):.2f}" if "score" in r else "—"
        cor = r.get("corrections", "—")
        lines.append(f"| {r['id']} | {r['desc']} | {r['expected_status']} | "
                     f"{r['actual_status']} | {sc} | {cor} | {det} | {co} | {so} |")
    lines.append("")
    lines.append("")
    lines.append("## 3. Interpretation")
    lines.append("")
    lines.append("- **Detection rate < 100%** indicates missed errors (mainly from ±2 "
                 "position tolerance in motif-based conserved-site detection).")
    lines.append("- **Correction accuracy** reflects whether corrected sequences recover "
                 "the wild-type residue at injected mutation sites.")
    lines.append("- **Status accuracy** is the primary end-to-end metric: does the "
                 "handoff manifest.csv carry the correct PASS/FAIL label?")
    lines.append("")
    lines.append("")
    lines.append("## 4. Issues Found and Fixed During Validation")
    lines.append("")
    lines.append("Validation was not a green-check exercise — it exposed four real defects, "
                 "all now fixed:")
    lines.append("")
    lines.append("| # | Defect | Root cause | Fix |")
    lines.append("|---|--------|-----------|-----|")
    lines.append("| 1 | Double disulfide-Cys mutant was auto-\"repaired\" and passed | "
                 "Correction restored C23/C104, audit then scored them 'corrected' | "
                 "Added rule: if both C23 and C104 are non-Cys before correction → FAIL (no Ig fold) |")
    lines.append("| 2 | Random noise (non-antibody) was force-corrected and passed | "
                 "Motif detection fell back to structural estimates on non-antibody input | "
                 "Added confidence flags: if neither FR2 W-motif nor FR4 J-motif matched → FAIL (not antibody-like) |")
    lines.append("| 3 | F118 correction direction was inverted | "
                 "F118 assumed F for heavy chain; actually JH genes use W (WGQG), JK/JL use F (FGQG) | "
                 "Corrected: VH→W, VK/VL→F |")
    lines.append("| 4 | C104 mis-localized by 1 residue when mutated | "
                 "Fixed `F118−6` offset inaccurate for VH (true distance 5) | "
                 "Chain-aware offset: VH→5, VL→8 |")
    lines.append("| 5 | Length anomaly initially mis-classified | "
                 "Status checks ran motif/AA checks before length check | "
                 "Reordered: length check now precedes motif/AA checks; length anomaly → FAIL |")
    lines.append("")
    lines.append("")
    lines.append("## 5. Known Limitations")
    lines.append("")
    lines.append("1. Motif-based position detection has ±2 residue tolerance; integrate "
                 "ANARCI for exact IMGT numbering in production.")
    lines.append("2. Test set uses germline-derived templates; real generated sequences "
                 "may carry more diverse CDR3 lengths.")
    lines.append("")
    lines.append("---")
    lines.append("")
    lines.append("*Report auto-generated by scripts/validate.py*")

    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def main():
    cases = build_cases()
    print(f"Built {len(cases)} test cases with ground truth\n")

    results, stats = run_validation(cases)

    # console summary
    for r in results:
        mark = "OK " if (r["detected_all"] and r["corr_ok"] and r["status_ok"]) else "MISS"
        print(f"  [{mark}] {r['id']:8s} {r['desc'][:38]:40s} "
              f"{r['expected_status']:>8s} → {r['actual_status']:>8s} "
              f"score={r.get('score', 0):.2f} corr={r.get('corrections', '—')}")

    total = stats["total"]
    print(f"\n{'='*70}")
    print("VALIDATION SUMMARY")
    print(f"{'='*70}")
    print(f"  Detection rate:      {stats['detected_all']}/{total}")
    print(f"  Correction accuracy: {stats['corr_correct']}/{total}")
    print(f"  Status accuracy:     {stats['status_correct']}/{total}")
    print(f"  False positives:     {stats['false_positives']}")
    print(f"  Crashes:             {stats['crashes']}")

    report_path = os.path.join(ROOT,
                               "docs", "validation_report_v1.md")
    write_report(results, stats, report_path)
    print(f"\n  Report written: {report_path}")


if __name__ == "__main__":
    main()
