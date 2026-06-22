# Volcano Plot — On-Treatment Biomarker Change vs. Baseline

A publication-ready volcano plot mapping mean percent change from baseline (x-axis) against statistical significance (−log₁₀ p-value, y-axis) for 11 circulating and tissue biomarkers measured across multiple on-treatment cycles in the ONCVIZ-001 treatment arm, with biomarkers color-coded by biological family and labeled at significance thresholds.

**Dataset:** ONCVIZ-001 · Treatment Arm · N = 62 patients · ADBM · 11 biomarkers × 6 cycles | **Language:** R · ggplot2 · ggrepel | **License:** CC BY 4.0

---

## The gap this fills

Volcano plots are ubiquitous in transcriptomics and proteomics but are underused in clinical biomarker reporting, where the dominant visualization is a waterfall or boxplot per marker. A multi-biomarker volcano allows simultaneous comparison of effect size and significance across biologically heterogeneous markers within the same trial arm — a format particularly relevant for basket trials with broad biomarker panels.

| | Prior art | This work |
|---|---|---|
| Scope | Single omic layer (RNA-seq, proteomics) | Multi-family panel: genomic · serum · immune · cytokine |
| x-axis metric | log₂ fold-change | Mean % change from baseline (clinically interpretable) |
| Statistical test | DESeq2/limma model | One-sample t-test vs. zero change (one-sample, no-change null) |
| Color grouping | Significance threshold | Biological family (4 color groups) |
| Point size | Uniform | Significant points enlarged (3.5 vs 2.2) |
| Shaded region | Not shown | Shaded band above p = 0.05 threshold |
| Labels | All or by threshold | `ggrepel` labels for −log₁₀(p) > 0.5 |
| Visit granularity | Summarized | Each biomarker × cycle is a separate point |

---

## Visual anatomy

```
  −log₁₀(p)
  4.0  ┤                    ●CEA (C4)
       │  shaded: p < 0.05  ◆PDL1 (C2)
  3.0  ┤                    ●
       │                        ● ctDNA (C3)
  2.0  ┤   ○        ○
       │
  1.0  ┤         ○     ○        ○
  ─────┼─────────────────────────────────── p = 0.05 (dashed)
  0.0  ┤   ○    ○    ○    ○
       ├────────────────────────────────────
      -30   0   +30  +60  +90  +110
              Mean % Change from Baseline
```

| Element | Description |
|---|---|
| Orange shaded band | Region above p = 0.05 significance threshold |
| Dashed horizontal line | p = 0.05 reference |
| Vertical line at 0 | No-change reference |
| Colored points | Biomarker family (Tumor genomics · Serum marker · Immune cell · Cytokine) |
| Large points | Statistically significant (p < 0.05) |
| `ggrepel` labels | Biomarker name + cycle, for −log₁₀(p) > 0.5 |

---

## Output figure

| File | Dimensions | DPI |
|---|---|---|
| `plots/03_volcano.png` | 13 × 9 in | 180 |

---

## Synthetic dataset — ONCVIZ-001 ADaM v1 (ADBM)

```
Arm            Treatment (n = 62)
Biomarkers     11: TMB · ctDNA · PDL1 · MSI · CEA · CA125 ·
                   CD8 · CD4 · NK · TREG · IFNγ
Cycles         C1 · C2 · C3 · C4 · C6 · C8
Points total   66 (11 × 6 cycles)
Seed           3 · fully reproducible
```

### Biomarker families

| Family | Biomarkers | Color |
|---|---|---|
| Tumor genomics | TMB · ctDNA · PDL1 · MSI | `#2c5f8a` (Blue) |
| Serum tumor marker | CEA · CA125 | `#c0392b` (Red) |
| Immune cell subset | CD8 · CD4 · NK · TREG | `#27ae60` (Green) |
| Cytokine | IFNγ | `#8e44ad` (Purple) |

### Variables used

| Variable | Description |
|---|---|
| `PARAMCD` | Biomarker code |
| `AVISIT` | On-treatment cycle (CYCLE 1–8 DAY 1) |
| `PCHG` | Percent change from baseline |
| `ARM` | Treatment · Control (filtered to Treatment only) |
| `ANL01FL` | Analysis flag (filtered to "Y") |

---

## Statistical method

Each point represents one (biomarker, cycle) combination. The test statistic is a **one-sample t-test** of individual patient PCHG values against a null hypothesis of zero change:

```r
pval <- t.test(PCHG_vector)$p.value
```

The mean PCHG across patients at that visit is the x-coordinate. This approach asks: *"Did this biomarker change significantly from baseline at this visit?"* It does not compare arms.

**Appropriate when:** Characterizing on-treatment pharmacodynamic activity in a single arm.
**Not appropriate when:** Comparing arms (use a two-sample test per biomarker); controlling family-wise error rate across 66 tests (apply Bonferroni or BH correction and annotate the adjusted threshold).

---

## Key parameters

| Parameter | Value | Description |
|---|---|---|
| Statistical test | One-sample t-test | vs. null: mean PCHG = 0 |
| Significance threshold | 0.05 | Nominal; no multiplicity correction |
| Label threshold | −log₁₀(p) > 0.5 | ~p < 0.32 |
| Point size: significant | 3.5 | `scale_size_manual` |
| Point size: non-significant | 2.2 | `scale_size_manual` |
| Point alpha | 0.85 | Slight transparency |
| X-axis limits | c(−32, 115) | Accommodates CEA/CA125 increases |
| Figure dimensions | 13 × 9 in | Near-square for balanced scatter |
| DPI | 180 | Publication quality |

---

## Limitations

- **No multiplicity correction:** With 66 simultaneous tests (11 × 6 cycles), ~3 false positives are expected by chance at α = 0.05. A Benjamini–Hochberg FDR line at q = 0.05 is recommended for formal reporting.
- **Cycle independence assumption:** The same patient contributes to multiple cycle points; the one-sample t-test treats each cycle independently and does not account for within-patient correlation. A mixed-effects model is more appropriate for longitudinal data.
- **Aggregate summary:** The x-axis shows mean PCHG; highly skewed distributions (common for serum markers) may make the median a more robust summary.

---

## Requirements

```r
library(ggplot2)   # >= 3.4
library(ggrepel)   # >= 0.9
library(dplyr)     # >= 1.1
```

---

## Files

| File | Description |
|---|---|
| `03_volcano.R` | Self-contained script · synthetic data generated internally |
| `plots/03_volcano.png` | Output figure · 13 × 9 in · 180 DPI |

---

## References

Gao J, et al. Integrative analysis of complex cancer genomics and clinical profiles using the cBioPortal. *Sci Signal.* 2013;6(269):pl1.

Ritchie ME, et al. limma powers differential expression analyses for RNA-sequencing and microarray studies. *Nucleic Acids Res.* 2015;43(7):e47.

Benjamini Y, Hochberg Y. Controlling the False Discovery Rate: A Practical and Powerful Approach to Multiple Testing. *J R Stat Soc Series B.* 1995;57(1):289–300.
