#!/usr/bin/env python3
"""
Download RCSB antibody structures, extract unique VH sequences, and ARCHIVE
them to data/positive_vh.json so validation runs read locally instead of
re-downloading (which caused sample-count fluctuations).

The 300 PDB IDs are already cached in data/pdb_antibody_ids.json.
"""

import concurrent.futures
import json
import os
import sys
import time
import urllib.request

_HERE = os.path.dirname(os.path.abspath(__file__))
_SCRIPTS_DIR = os.path.dirname(_HERE)      # scripts/
ROOT = os.path.dirname(_SCRIPTS_DIR)       # project root
DATA = os.path.join(ROOT, "data")

UA = {"User-Agent": "Mozilla/5.0 (research script)"}
CH1 = ["ASTKGP", "ASTKGPS", "SASTKG", "ASSKGP", "AKTTAP", "AKTTPP", "KTTPP",
       "ASVAAP", "ASPTSP", "ASTKAP", "ASTQSPS", "STQSPS", "APTKAP", "ASATL"]
VH_SIGS = ("EVQ", "QVQ", "EVK", "QVK", "EVE", "QVE", "EVO")


def fetch_fasta(pdb_id, retries=6):
    url = f"https://www.rcsb.org/fasta/entry/{pdb_id}"
    for attempt in range(retries):
        try:
            req = urllib.request.Request(url, headers=UA)
            with urllib.request.urlopen(req, timeout=40) as r:
                return pdb_id, r.read().decode("utf-8")
        except Exception:
            time.sleep(2 + attempt * 2)
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
    for m in CH1:
        i = seq.find(m)
        if i > 0:
            cut = i
            break
    head = seq[:cut]
    for sig in VH_SIGS:
        i = head.find(sig)
        if 0 <= i < 40:
            vh = head[i:]
            return vh if len(vh) <= 140 else vh[:125]
    return head if 90 <= len(head) <= 140 else head[:125]


def is_vh(seq):
    return any(sig in seq[:15] for sig in VH_SIGS)


def main():
    pdb_ids = json.load(open(os.path.join(DATA, "pdb_antibody_ids.json")))
    print(f"开始下载 {len(pdb_ids)} 个结构 ...", flush=True)

    seen = {}  # vh -> [pdb, ...]
    with concurrent.futures.ThreadPoolExecutor(max_workers=10) as ex:
        for pdb, text in ex.map(fetch_fasta, pdb_ids):
            if text is None:
                continue
            for h, s in zip(*parse_fasta(text)):
                if "heavy" not in h.lower():
                    continue
                vh = extract_vh(s)
                if 90 <= len(vh) <= 140 and is_vh(vh):
                    seen.setdefault(vh, []).append(pdb)

    # 存档
    out = {"vh_sequences": list(seen.keys()),
           "sources": {vh: src for vh, src in seen.items()},
           "count": len(seen)}
    with open(os.path.join(DATA, "positive_vh.json"), "w", encoding="utf-8") as f:
        json.dump(out, f, ensure_ascii=False)

    print(f"提取独立 VH {len(seen)} 条，已存档到 data/positive_vh.json", flush=True)


if __name__ == "__main__":
    main()
