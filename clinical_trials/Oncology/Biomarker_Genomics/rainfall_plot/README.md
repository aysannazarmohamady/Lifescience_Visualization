# Rainfall Plot — Inter-Mutation Distance Across the Sequenced Gene Panel

A publication-ready rainfall plot displaying the log₁₀-transformed distance between consecutive somatic mutations (ordered by chromosomal cytoband position within each chromosome) across 138 somatic variants in 66 patients. Points positioned below a cluster threshold line indicate locally hypermutated regions or mutational hotspots.

**Dataset:** ONCVIZ-001 · N = 138 somatic point mutations / indels · 66 patients · 15 genes | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

Rainfall plots were developed by Alexandrov et al. (2013) to visualize kataegis — localized hypermutation events characterized by clusters of C→T and C→G mutations in TpC context at short inter-mutation distances. They are generated routinely by tools such as maftools and SigProfiler but require whole-genome sequencing data for full kataegis detection. Applying the rainfall format to targeted panel data (where inter-mutation distance is computed between variants in the sequenced gene panel, not genome-wide) provides a panel-level view of mutation clustering that is informative even without WGS.

| | WGS kataegis rainfall | This work (panel-level) |
|---|---|---|
| Data type | Whole-genome variants | Targeted panel variants (n=138) |
| x-axis | Genome-wide position (Mb) | Ordinal variant index within chromosome |
| y-axis | log₁₀(genomic distance in bp) | log₁₀(protein position distance in AA) |
| Clustering definition | Distance < 1 kb | Illustrative threshold at distance < 300 AA |
| Color | SBS mutation context (96-type) | Variant classification (6 types) |
| Kataegis detection | Formal algorithm | Visual annotation only |
| Reproducibility | Tool-dependent | `set.seed(10)` · bitwise-identical |

---

## Visual anatomy

```
  log₁₀(distance)
  3.0  ┤  ○  ○      ○   ○    ○       ○   ○    ○   ○
       │    ○    ○      ○  ○   ○  ○     ○
  2.0  ┤       ○        ○           ○       ○
       │  ○         ○     ○    ○     ○
  1.0  ┤    ○                         ● ● (cluster)
  ─────┼──────────────────────────── threshold log₁₀(300) ≈ 2.48 (red dashed)
  0.5  ┤  ● (cluster)
       ├──────────────────────────────────────────────
       chr1 chr3  chr7 chr9 chr10 chr12 chr13 chr17 chr18 chr19
```

| Element | Description |
|---|---|
| Each point | One somatic variant · y = log₁₀(distance to previous variant in same chr) |
| Red dashed line | Illustrative cluster threshold at distance = 300 AA |
| Point color | Variant classification (6-class palette) |
| X-axis | Chromosome labels at group midpoints (ordinal index, not bp) |
| Missing first point | The first variant on each chromosome has no predecessor → `NA` dropped |

---

## Output figure

| File | Dimensions | DPI |
|---|---|---|
| `plots/10_rainfall.png` | 15 × 6 in | 180 |

---

## Synthetic dataset — ONCVIZ-001 ADaM v1 (ADMUT)

```
n mutations    138 somatic point mutations / indels
n patients     66 (patients with ≥1 mutation)
n genes        15
Chromosomes    10 (same as Manhattan/Miami plots)
Seed           10 · fully reproducible
```

### Variables used

| Variable | Source | Description |
|---|---|---|
| `HUGO_SYMBOL` | ADMUT | Gene name (used for CHR assignment) |
| `PROTEIN_POS` | ADMUT | Amino acid position (surrogate genomic coordinate) |
| `VARIANT_CLASS` | ADMUT | Mutation type (color encoding) |
| `CHR` | Derived | Chromosome from `CHR_MAP[HUGO_SYMBOL]` |

### Variant class colors

| Variant class | Color |
|---|---|
| Missense Mutation | `#2c5f8a` (Blue) |
| Nonsense Mutation | `#c0392b` (Red) |
| Frame Shift (Del) | `#e67e22` (Orange) |
| Frame Shift (Ins) | `#f1c40f` (Yellow) |
| Splice Site | `#8e44ad` (Purple) |
| In-Frame Del | `#27ae60` (Green) |
| Other Mutation | `#7f8c8d` (Gray) |

---

## Distance computation

```r
DIST <- abs(PROTEIN_POS − lag(PROTEIN_POS))   # within chromosome
DIST <- pmax(DIST, 1)                          # prevent log(0)
LOG_DIST <- log10(DIST)
```

Variants are first sorted by `CHR` (chromosomal order) then by `PROTEIN_POS` within chromosome. The first variant on each chromosome has `DIST = NA` and is excluded from the plot.

**Important caveat:** `PROTEIN_POS` is the amino acid position within the gene, not the chromosomal base-pair coordinate. For genes on the same chromosome (e.g., EGFR, MET, BRAF on chr7), the distance between the last variant in EGFR and the first in MET reflects protein-coordinate discontinuity, not genomic distance. This is explicitly noted in the figure caption as "genomic-position proxy."

---

## Key parameters

| Parameter | Value | Description |
|---|---|---|
| Cluster threshold | 300 AA | `log₁₀(300) ≈ 2.48` · illustrative |
| Distance floor | 1 | `pmax(DIST, 1)` prevents `log10(0)` |
| Sort order | CHR → PROTEIN_POS | Within-chromosome sort |
| Point size | 2.0 | Moderate; 138 points total |
| Point alpha | 0.85 | Slight transparency |
| Figure dimensions | 15 × 6 in | Wide landscape (standard for rainfall plots) |
| DPI | 180 | Publication quality |

---

## Limitations

- **Protein-position proxy:** True kataegis analysis requires genomic base-pair coordinates from VCF files. The protein position proxy used here introduces two artifacts: (1) inter-gene distances within the same chromosome are not genomically meaningful, and (2) non-coding variants (splice sites, UTR mutations) may be mispositioned.
- **Targeted panel ≠ WGS:** With only 15 genes sequenced, inter-mutation distances are inherently sparse. Kataegis events in non-panel regions are invisible. This plot is best interpreted as a panel-level mutation clustering overview, not a kataegis detector.
- **Low point density:** 138 variants across 10 chromosome groups produce a sparse plot. Rainfall plots are most informative with ≥1,000 variants (i.e., WGS or large panel data).

---

## Requirements

```r
library(ggplot2)   # >= 3.4
library(dplyr)     # >= 1.1
```

---

## Files

| File | Description |
|---|---|
| `10_rainfall.R` | Self-contained script · synthetic data generated internally |
| `plots/10_rainfall.png` | Output figure · 15 × 6 in · 180 DPI |

---

## References

Alexandrov LB, et al. Signatures of mutational processes in human cancer. *Nature.* 2013;500(7463):415–421.

Roberts SA, et al. Clustered mutations in cancer: APOBEC's role in breast tumor evolution. *Nat Genet.* 2013;45(9):970–976.

Mayakonda A, et al. Maftools: efficient and comprehensive analysis of somatic variants in cancer. *Genome Res.* 2018;28(11):1747–1756.
