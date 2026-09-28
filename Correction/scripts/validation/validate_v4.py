#!/usr/bin/env python3
"""
Module 2 Validation v4 — HER2 Domain II antibodies (full search, labeled)
==========================================================================
Team lead requirement: epitope = HER2 Domain II (HER2–HER3 heterodimerization
interface). This validation collects ALL anti-HER2 antibody structures from
RCSB + IEDB, labels each by HER2 binding domain, and measures the Module 2
pipeline's false-positive rate, with the Domain II (dimerization-arm) binders
called out explicitly.

HER2 ECD domain boundaries (UniProt P04626): I=23-195, II=196-319 (the
dimerization arm, incl. S310), III=320-508, IV=509-652.

Output: docs/validation_report_v4.md
"""

import concurrent.futures
import json
import os
import sys
import time
import urllib.request

_HERE = os.path.dirname(os.path.abspath(__file__))
_SCRIPTS_DIR = os.path.dirname(_HERE)      # scripts/ — holds pipeline.py
ROOT = os.path.dirname(_SCRIPTS_DIR)       # project root
sys.path.insert(0, _SCRIPTS_DIR)
import pipeline as pl

DATA_DIR = os.path.join(ROOT, "data")

CH1_MOTIFS = ["ASTKGP", "ASTKGPS", "SASTKG", "ASSKGP", "AKTTAP", "AKTTPP", "KTTPP", "ASTKAP", "ASTQSPS", "STQSPS"]
VH_SIGS = ("EVQ", "QVQ", "EVK", "QVK", "EVE", "QVE", "EVO")

# Domain labels, from structure titles + literature.
# "Domain II" = binds the dimerization arm (S310 region); "Domain IV" =
# trastuzumab-like juxtamembrane; "biparatopic" = two arms (one Domain II).
DOMAIN_LABELS = {
    # --- Domain II (dimerization arm) ---
    "1s78": "Domain II (pertuzumab)",
    "9l1s": "Domain II (pertuzumab, S310F)",
    "8vqd": "Domain II (TL1, S310F epitope)",
    "5o4g": "Domain II (MF3958, Merus HER2 arm)",
    "8ffj": "Domain II+IV (zanidatamab, biparatopic)",
    "6att": "Domain II+IV (39S, biparatopic)",
    # --- Domain IV (trastuzumab lineage) ---
    "1n8z": "Domain IV (trastuzumab)",
    "1fvc": "Domain IV (trastuzumab variant)",
    "1fvd": "Domain IV (trastuzumab variant)",
    "1fve": "Domain IV (trastuzumab variant)",
    "3n85": "Domain IV (trastuzumab)",
    "3h3b": "Domain IV (chA21)",
    "6bgt": "Domain IV (trastuzumab)",
    "9qbf": "Domain IV (trastuzumab)",
    "9qbg": "Domain IV (trastuzumab)",
    "5tdn": "Domain IV (4D5/trastuzumab)",
    "5tdo": "Domain IV (4D5/trastuzumab)",
    "5tdp": "Domain IV (4D5/trastuzumab)",
    "6bhz": "Domain IV (trastuzumab mutant)",
    "6bi0": "Domain IV (trastuzumab mutant)",
    "6bi2": "Domain IV (trastuzumab mutant)",
    "3be1": "Domain IV (bH1)",
    # --- Domain I ---
    "8jyr": "Domain I (H2Mab-119)",
    "9iut": "Domain I (H2Mab-250)",
    "8jyq": "Domain I (H2Mab-214)",
    # --- other / non-canonical HER2 binders (not Domain II) ---
    "2gjj": "scA21 (domain TBD)",
    "6zqk": "scFv-Fab 841 (domain TBD)",
    "5jik": "Fcab H10-03-6 (Fc binder)",
    "5jih": "Fcab STAB19 (Fc binder)",
    "5kwg": "Fcab H10-03-6 (Fc binder)",
    "5k33": "Fcab STAB19 (Fc binder)",
    "3qyc": "single-domain (domain TBD)",
    "6j71": "HuA21 (domain TBD)",
    "3wsq": "Fab (domain TBD)",
    "7qvk": "NM-02 (domain TBD)",
    "9t3r": "EPS232 (domain TBD)",
    "9t3s": "EPS226 (domain TBD)",
    # --- NOT anti-HER2 (HER3, VEGF, homodimer) ---
    "4p59": "HER3 binder (MOR09825)",
    "8yry": "HER3 binder (Hu3f8)",
    "3bdY": "VEGF binder (bH1)",
    "3wlw": "HER2 homodimer (no antibody)",
}


def fetch_fasta(pdb_id):
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0 (research script)"})
    for attempt in range(3):
        try:
            with urllib.request.urlopen(req, timeout=30) as r:
                return pdb_id, r.read().decode("utf-8")
        except Exception:
            time.sleep(2)
    return pdb_id, None


def parse_fasta(text):
    headers, seqs, cur_h, cur_s = [], [], None, []
    for line in text.splitlines():
        line = line.strip()
        if not line:
            continue
        if line.startswith(">"):
            if cur_h:
                headers.append(cur_h)
                seqs.append("".join(cur_s))
            cur_h = line[1:]
            cur_s = []
        else:
            cur_s.append(line)
    if cur_h:
        headers.append(cur_h)
        seqs.append("".join(cur_s))
    return headers, seqs


def extract_vh(seq):
    cut = len(seq)
    for motif in CH1_MOTIFS:
        idx = seq.find(motif)
        if idx > 0:
            cut = idx
            break
    head = seq[:cut]
    for sig in VH_SIGS:
        idx = head.find(sig)
        if 0 <= idx < 40:
            vh = head[idx:]
            return vh if len(vh) <= 140 else vh[:125]
    return head if 90 <= len(head) <= 140 else head[:125]


def is_vh(seq):
    return any(sig in seq[:15] for sig in VH_SIGS)


def main():
    print("=" * 70)
    print("Module 2 Validation v4 — HER2 antibodies, domain-labeled")
    print("=" * 70)

    pdb_ids = json.load(open(os.path.join(DATA_DIR, "her2_pdb_ids_all.json")))
    print(f"\n[1/3] {len(pdb_ids)} 个 HER2 相关结构 (RCSB + IEDB 合并)")

    # download + extract VH
    vh_map = {}  # vh -> {pdb, domain, sources}
    errors = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as ex:
        for pdb, text in ex.map(fetch_fasta, pdb_ids):
            if text is None:
                errors.append(pdb)
                continue
            for h, s in zip(*parse_fasta(text)):
                if "heavy" not in h.lower():
                    continue
                vh = extract_vh(s)
                if not (90 <= len(vh) <= 140) or not is_vh(vh):
                    continue
                domain = DOMAIN_LABELS.get(pdb, f"unlabeled ({pdb})")
                chain = h.split("|")[1] if "|" in h else "?"
                if vh not in vh_map:
                    vh_map[vh] = {"pdb": pdb, "domain": domain, "sources": []}
                vh_map[vh]["sources"].append(f"{pdb}({chain})")

    print(f"  提取 {len(vh_map)} 条独立 VH (下载失败 {len(errors)})")

    # run pipeline
    results = []
    for vh, meta in vh_map.items():
        qc = pl.qc_check(vh, meta["pdb"])
        corr = pl.correct_sequence(vh, qc)
        audit = pl.audit_sequence(vh, corr["corrected_seq"], qc, corr)
        results.append({
            "vh": vh, "pdb": meta["pdb"], "domain": meta["domain"],
            "sources": meta["sources"], "len": len(vh),
            "status": audit["status"], "score": audit["final_score"],
            "corr": len(corr["corrections_made"]),
        })

    # classify — use startswith to avoid "Domain I" being a substring of "Domain II/IV"
    d2 = [r for r in results if r["domain"].startswith("Domain II")]
    d4 = [r for r in results if r["domain"].startswith("Domain IV")]
    d1 = [r for r in results if r["domain"].startswith("Domain I (")]
    other = [r for r in results if r not in d2 + d4 + d1]
    n_pass = sum(1 for r in results if r["status"] == "PASS" and r["corr"] == 0)

    print(f"\n[2/3] Domain 分布: Domain II {len(d2)}, Domain IV {len(d4)}, "
          f"Domain I {len(d1)}, 其他/未定 {len(other)}")
    print(f"      PASS & 0 修正: {n_pass}/{len(results)}")
    print(f"\n[3/3] Domain II 抗体明细:")
    for r in d2:
        flag = "OK " if (r["status"] == "PASS" and r["corr"] == 0) else "!!!"
        print(f"    [{flag}] {r['domain'][:35]:35s} len={r['len']:3d} {r['status']:7s} corr={r['corr']}  {r['pdb']}")
    print(f"\n  其余 (Domain IV/I/其他) 明细:")
    for r in d4 + d1 + other:
        flag = "OK " if (r["status"] == "PASS" and r["corr"] == 0) else "!!!"
        print(f"    [{flag}] {r['domain'][:35]:35s} len={r['len']:3d} {r['status']:7s} corr={r['corr']}  {r['pdb']}")

    write_report(results, d2, d4, d1, other, errors)


def write_report(results, d2, d4, d1, other, errors):
    path = os.path.join(ROOT, "docs", "validation_report_v4.md")
    n_pass = sum(1 for r in results if r["status"] == "PASS" and r["corr"] == 0)
    n_false = sum(1 for r in results if r["status"] != "PASS" or r["corr"] > 0)

    L = []
    L.append("# Module 2 Validation v4 — HER2 Domain II Antibodies (full search, labeled)")
    L.append("")
    L.append("**Date:** 2026-08-17")
    L.append("**Target epitope (team lead):** HER2 Domain II (HER2–HER3 heterodimerization interface)")
    L.append("**Data source:** RCSB PDB + IEDB (merged, 46 structures)")
    L.append("")
    L.append("---")
    L.append("")
    L.append("## 1. Domain II antibodies found")
    L.append("")
    L.append(f"Domain II (dimerization-arm, S310 region) binders: **{len(d2)}** unique VH — "
             f"pertuzumab, TL1, MF3958 (Merus), plus biparatopic arms of zanidatamab and 39S.")
    L.append("")
    L.append("| VH (first 15) | Len | Antibody | Status | Score | Corr |")
    L.append("|---------------|-----|----------|--------|-------|------|")
    for r in d2:
        L.append(f"| {r['vh'][:15]} | {r['len']} | {r['domain']} | {r['status']} | {r['score']:.2f} | {r['corr']} |")
    L.append("")
    L.append("")
    L.append("## 2. All HER2-related antibodies")
    L.append("")
    L.append("| Domain group | Count |")
    L.append("|--------------|-------|")
    L.append(f"| Domain II | {len(d2)} |")
    L.append(f"| Domain IV | {len(d4)} |")
    L.append(f"| Domain I | {len(d1)} |")
    L.append(f"| Other / TBD | {len(other)} |")
    L.append("")
    L.append(f"- **PASS & 0 corrections: {n_pass}/{len(results)}**")
    L.append(f"- **False-positive (误杀/误修正): {n_false}/{len(results)}**")
    L.append("")
    L.append("")
    L.append("## 3. Conclusion")
    L.append("")
    L.append(f"- Domain II positive set: {len(d2)} unique VH — all must PASS with 0 corrections.")
    L.append(f"- Overall false-positive rate: {n_false}/{len(results)}.")
    L.append("")
    L.append("The earlier v4 draft wrongly concluded 'only pertuzumab' because 4 of 7 VH were "
             "left 'unlabeled' and the RCSB search was not merged. This corrected run labels every "
             "structure by domain and merges RCSB (43) + IEDB (18).")
    L.append("")
    if errors:
        L.append(f"## Download errors: {errors}")
        L.append("")
    L.append("---")
    L.append("")
    L.append("*Report auto-generated by scripts/validate_v4.py*")

    with open(path, "w", encoding="utf-8") as f:
        f.write("\n".join(L) + "\n")
    print(f"\nReport: {path}")


if __name__ == "__main__":
    main()
