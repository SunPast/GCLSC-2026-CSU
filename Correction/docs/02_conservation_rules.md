# Module 2: Antibody Sequence Correction — Conservative Site Rules

## 1. Literature Basis

### 1.1 IMGT Numbering Scheme

The IMGT (ImMunoGeneTics) unique numbering system for immunoglobulin variable domains is the international standard for antibody sequence annotation (Lefranc et al., 2003). Key boundaries:

| Region | IMGT Positions |
|--------|---------------|
| FR1    | 1–26          |
| CDR1   | 27–38         |
| FR2    | 39–55         |
| CDR2   | 56–65         |
| FR3    | 66–104        |
| CDR3   | 105–117       |
| FR4    | 118–128       |

**Reference:** Lefranc, M. P., et al. (2003). "IMGT unique numbering for immunoglobulin and T cell receptor variable domains and Ig superfamily V-like domains." *Developmental & Comparative Immunology*, 27(1), 55–77. DOI: 10.1016/S0145-305X(02)00039-3

### 1.2 Critical Conserved Framework Positions

The following positions are structurally essential for the immunoglobulin fold and are conserved at >98% across all functional antibody V domains.

#### Tier 1: Structurally Essential (nearly 100% conserved — mutation → FAIL)

| IMGT Pos | Region | Required AA | Structural Role | Reference |
|----------|--------|-------------|-----------------|-----------|
| C23      | FR1 (strand B) | C (Cys)     | Intra-domain disulfide bridge (S–S bond with C104). 1st-CYS. Essential for Ig fold stability. | Lefranc et al., 2003; Chothia & Lesk, 1987 |
| C104     | FR3 (strand F) | C (Cys)     | Intra-domain disulfide bridge (S–S bond with C23). 2nd-CYS. Breakage causes domain unfolding. | Lefranc et al., 2003; Chothia & Lesk, 1987 |

**Reference:** Chothia, C., & Lesk, A. M. (1987). "Canonical structures for the hypervariable regions of immunoglobulins." *Journal of Molecular Biology*, 196(4), 901–917. DOI: 10.1016/0022-2836(87)90412-8

#### Tier 2: Highly Conserved (>95% — mutation → CORRECT)

| IMGT Pos | Region | Required AA | Structural Role | Reference |
|----------|--------|-------------|-----------------|-----------|
| W41      | FR2    | W (Trp)     | Core hydrophobic packing; VH/VL interface. Mutation disrupts domain folding. | Ewert et al., 2003; Honegger & Plückthun, 2001 |
| F118/W118| FR4    | F or W      | J-region aromatic residue (IMGT "J-PHE or J-TRP"). **W in JH (heavy chain), F in JK/JL (light chain).** Critical for Ig fold C-terminal closure. | Lefranc, 2014 |
| G119     | FR4    | G (Gly)     | Part of J-motif "F/W-G-X-G". Tight turn required; any side chain would clash. | Lefranc et al., 2003 |

#### Tier 3: Moderately Conserved (>85% — mutation → WARN)

| IMGT Pos | Region | Commonly Found | Structural Role | Reference |
|----------|--------|---------------|-----------------|-----------|
| L89      | FR3    | L, I, V       | Hydrophobic core residue; stabilizes VH/VL packing. | Ewert et al., 2003 |
| D86/Y90  | FR3    | Variable      | VH/VL interface salt bridge participants (varies by germline). | Vargas-Madrazo & Paz-García, 2002 |
| T122     | FR4    | T, S          | J-region conserved residue. | Lefranc et al., 2003 |

**References:**
- Ewert, S., Huber, T., Honegger, A., & Plückthun, A. (2003). "Biophysical properties of human antibody variable domains." *Journal of Molecular Biology*, 325(3), 531–553. DOI: 10.1016/S0022-2836(02)01237-8
- Honegger, A., & Plückthun, A. (2001). "Yet another numbering scheme for immunoglobulin variable domains: an automatic modeling and analysis tool." *Journal of Molecular Biology*, 309(3), 657–670. DOI: 10.1006/jmbi.2001.4662
- Lefranc, M. P. (2014). "Immunoglobulin and T cell receptor genes: IMGT and the birth of immunoinformatics." *Frontiers in Immunology*, 5, 22. DOI: 10.3389/fimmu.2014.00022
- Vargas-Madrazo, E., & Paz-García, E. (2002). "Modifications in the VH/VL interface of antibodies." *Journal of Molecular Recognition*, 15(5), 291–299.

### 1.3 Non-Standard Amino Acids

Standard proteinogenic amino acids (20):
```
A C D E F G H I K L M N P Q R S T V W Y
```

Non-standard codes encountered in sequence data and their handling:

| Code | Name | Replacement | Rationale |
|------|------|-------------|-----------|
| B    | Asx (Asp or Asn) | N | Asn more common at antibody surface positions; conservative choice |
| Z    | Glx (Glu or Gln) | E | Glu more common at antibody surface positions |
| J    | Xle (Leu or Ile) | L | Leu more common in antibody framework hydrophobic core |
| U    | Sec (Selenocysteine) | C | Chemically closest; Cys forms native disulfide bonds |
| O    | Pyl (Pyrrolysine) | K | Closest standard analog; both basic |
| X    | Unknown/Any | A | Ala is smallest, sterically least disruptive; "safe fallback" |
| *    | Stop codon | — | Remove; sequence ends at stop |
| -    | Gap | — | Remove; represents deletion |

**Reference:** IUPAC-IUB Joint Commission on Biochemical Nomenclature (1984). "Nomenclature and Symbolism for Amino Acids and Peptides." *European Journal of Biochemistry*, 138(1), 9–37.

### 1.4 Additional Sequence Quality Rules

| Rule | Condition | Action | Reference |
|------|-----------|--------|-----------|
| Length check | V domain length 95–135 AA | Outside range → WARN; <90 or >140 → FAIL | Lefranc et al., 2003 |
| N-glycosylation in CDR | N-X-S/T motif (X≠P) in CDR1/2/3 | WARN — may affect antigen binding | Jefferis, R. (2009). *Nature Reviews Drug Discovery*, 8(3), 226–234. |
| Unpaired Cys | Odd number of C residues | WARN — free cysteine may cause aggregation | Buchanan, A., et al. (2013). *mAbs*, 5(2), 255–262. |
| Proline in CDR | Multiple P in CDR loop | WARN — may constrain loop flexibility | Collis, A. V., et al. (2003). *Journal of Molecular Biology*, 325(2), 337–354. |
| Deamidation motif | N-G in CDR | WARN — potential deamidation hot spot | Sydow, J. F., et al. (2014). *PLoS ONE*, 9(6), e100736. |
| Isomerization motif | D-G, D-S, D-D in CDR | WARN — potential isomerization hot spot | Wakankar, A. A., & Borchardt, R. T. (2006). *Journal of Pharmaceutical Sciences*, 95(11), 2321–2336. |
| Oxidation motif | M in CDR | WARN — potential methionine oxidation | Gao, J., et al. (2015). *Molecular Pharmaceutics*, 12(7), 2382–2393. |

## 2. Correction Strategy (Evidence-Based)

### 2.1 Auto-Correction Rules (applied by pipeline)

| Issue | Action | Confidence |
|-------|--------|------------|
| Non-standard AA (B, Z, J, U, O, X) | Replace with mapped standard AA | High |
| Stop codon (*) / Gap (-) in V domain | Remove residue | High |
| C23 not Cys → mutated | Correct to C → PASS with flag | High — C23 is invariant |
| C104 not Cys → mutated | Correct to C → PASS with flag | High — C104 is invariant |
| W41 not Trp → mutated | Correct to W → PASS with flag | High — >99% conserved |
| F118 not F/W | Correct to W (VH/JH) or F (VL/JK-JL) → PASS with flag | Medium — depends on J-gene |
| G119 not Gly → mutated | FAIL — cannot reliably correct | N/A — tight turn, any change unpredictable |

### 2.2 Sequences Marked FAIL

Hard failure conditions — any one of these is FAIL regardless of score:

- G119 mutation — cannot predict structural consequences, not correctable
- Length <90 or >140 — not a valid V domain
- Not antibody-like — neither the FR2 W-motif nor the FR4 J-motif is found
- Missing both conserved Cys (C23 and C104) — no Ig fold possible
- Missing either C23 or C104 alone — a disulfide bond is broken

### 2.3 Audit Scoring Rubric

| Criterion | Weight | Scoring |
|-----------|--------|---------|
| No non-standard AA | 15% | Pass: 1.0, Fail: 0.0 |
| C23 intact | 25% | Pass/CORRECTED: 1.0, Fail: 0.0 |
| C104 intact | 25% | Pass/CORRECTED: 1.0, Fail: 0.0 |
| W41 intact | 10% | Pass/CORRECTED: 1.0, Fail: 0.0 |
| F/W118 intact | 10% | Pass/CORRECTED: 1.0, Fail: 0.0 |
| Length normal | 10% | 95–135: 1.0, 90–94 or 136–140: 0.5, outside: 0.0 |
| Other quality rules | 5% | No uncorrectable issue: 1.0, warnings: 0.5, G119 mutated: 0.0 |

A position that is correct **after** correction scores full weight; whether it was
native or corrected is recorded in `audit_log.json` only, not scored.

| Final Score | Status |
|-------------|--------|
| ≥ 0.94 | PASS |
| < 0.94 | FAIL |

The 0.94 cutoff is the highest threshold with zero false rejections on the
calibration set (positive minimum score 0.950); see
`docs/02_规则说明_中文.md` Section 5. A sequence that trips a hard failure
condition (Section 2.2) is FAIL regardless of score. There is no intermediate
WARN status — `manifest.csv` carries PASS or FAIL only.

## 3. Position Identification Approach

Conserved positions are located by **exact IMGT numbering via ANARCI** (Dunbar & Deane, 2016), run through a pyhmmer backend (no hmmscan binary required). ANARCI's HMM-based numbering is motif-independent, so it correctly locates C23/C104/W41/F118/G119 even when a mutation destroys the local motif (e.g. G119→E, double-Cys breaks).

**ANARCI ↔ motif complement:**

- ANARCI drifts by ~1 position on rare V genes, and by more on long CDR3 (>13 match states). It can also slip by one match state on an N-terminally truncated sequence: deleting exactly 15 residues from a full-length VH stacks the IMGT gap at position 10 with the alignment's own gap, after which IMGT 23 reads `S` and IMGT 104 reads `Y` although both cysteines are present (see `docs/02_IMGT注释说明.md` Section 4 for the controlled experiment). In those cases the intact motifs below are used to correct it.
- The motif fallback (used only when ANARCI is unavailable or its numbering is implausible) is:
  1. **C23**: the first Cys in FR1 (raw pos 20–26)
  2. **C104**: the Cys pairing with C23 (raw pos ~92–96)
  3. **W41**: conserved W within the cross-species FR2 motif `W[VILFY][RKQ]` (human WVRQ, mouse WVKQ, rabbit WIRK, light-chain WYQQ)
  4. **FR4 J-motif** `[FW]-G-X-G`: the F/W at position 118 and G at position 119

**Reference:** Dunbar, J., & Deane, C. M. (2016). "ANARCI: antigen receptor numbering and receptor classification." *Bioinformatics*, 32(2), 298–300. DOI: 10.1093/bioinformatics/btv552

---

*Document Version: 1.1 | Last Updated: 2026-08-17*
*Prepared for: GCLSC Project — Module 2 (Structure Correction & Multi-Agent QC)*
