# VAF Distribution Plot — Clonal Architecture

A three-panel publication-ready figure characterizing the distribution of variant allele frequencies (VAF) across 138 somatic mutations in 66 patients: (top) a cohort-wide stacked histogram by clonality status; (bottom-left) per-tumor-type violin + boxplot + jitter with clonality overlay; (bottom-right) per-gene violin + boxplot + jitter for the 8 most frequently mutated genes.

**Dataset:** ONCVIZ-001 · N = 138 somatic point mutations / indels · 66 patients · ADMUT | **Language:** R · ggplot2 · patchwork | **License:** CC BY 4.0

---

## The gap this fills

VAF distribution plots are a standard quality-control and biological interpretation tool in cancer genomics: they reveal tumor purity, ploidy shifts, the presence of subclonal populations, and potential contamination. Most publications show only a single histogram or boxplot; combining histology-level and gene-level breakdowns in a synchronized three-panel layout provides a comprehensive clonal architecture overview.

| | Prior art | This work |
|---|---|---|
| Panel count | Single histogram or single boxplot | Three synchronized panels |
| Clonality encoding | Not shown | Stacked histogram + jitter color by clonal flag |
| Clonality threshold | Hard VAF = 0.5 line | Dashed red line at VAF = 0.4 with annotation |
| Tumor-type breakdown | Separate figure or facet | Bottom-left panel · violin + box + jitter |
| Gene-level breakdown | Not shown | Bottom-right panel · 8 genes · violin + box + jitter |
| Patchwork layout | N/A | `/ (| )` operator · heights `c(1, 1.2)` |
| Reproducibility | No fixed seed | `set.seed(12)` · bitwise-identical |

---

## Visual anatomy

```
  [TOP PANEL — Cohort-wide histogram]
  15 ┤         ████
  10 ┤      ████████                    (blue = clonal)
   5 ┤  ████████████████████            (light blue = subclonal)
   0 └─────────────────────────────── VAF
     0.0  0.2  0.4  0.6  0.8   ← VAF = 0.4 (red dashed)

  [BOTTOM-LEFT]          [BOTTOM-RIGHT]
  Violin + box:          Violin + box:
  NSCLC BRCA HCC CRC    PTEN KEAP1 PIK3CA SMAD4 CDKN2A KRAS ARID1A TP53
  PDAC                   Points colored by clonality
```

| Element | Description |
|---|---|
| Top histogram | `geom_histogram(binwidth=0.05, position="stack")` · color by `CLONAL` |
| Dashed red line | VAF = 0.4 illustrative clonal threshold in all three panels |
| Violin | `geom_violin(alpha=0.38, color=NA)` · fill by tumor type or gray |
| Boxplot | `geom_boxplot(width=0.12, outlier.shape=NA)` · IQR and median |
| Jitter | `geom_jitter(width=0.08, size=1.8, alpha=0.72)` · colored by `CLONAL` |
| Clonal color | `#1a3f8f` (dark blue) = clonal · `#aec6ef` (light blue) = subclonal |

---

## Output figure

| File | Dimensions | DPI |
|---|---|---|
| `plots/12_vaf.png` | 17 × 12 in | 180 |

---

## Synthetic dataset — ONCVIZ-001 ADaM v1 (ADMUT)

```
n mutations    138 somatic point mutations / indels
n patients     66
VAF range      0.05 – 0.80
Variant type   Point mutations + indels (no CNV)
Seed           12 · fully reproducible
```

### Variables used

| Variable | Source | Description |
|---|---|---|
| `USUBJID` | ADMUT | Subject identifier |
| `HUGO_SYMBOL` | ADMUT | Gene name |
| `VAF` | ADMUT | Variant allele frequency (0–1) |
| `CLONAL` | ADMUT | "Y" = clonal · "N" = subclonal |
| `TUMORTYPE` | ADSL (joined) | Histology |
| `VARIANT_CLASS` | ADMUT | Mutation type (not used directly in this figure) |

### Clonality colors

| Status | Color | Description |
|---|---|---|
| Clonal (Y) | `#1a3f8f` | Dark blue — early, high-CCF mutations |
| Subclonal (N) | `#aec6ef` | Light blue — late, low-CCF mutations |

### Gene-level panel: 8 genes shown

`PTEN · KEAP1 · PIK3CA · SMAD4 · CDKN2A · KRAS · ARID1A · TP53`

---

## Key parameters

| Parameter | Value | Description |
|---|---|---|
| Histogram binwidth | 0.05 | 20 bins across [0, 1] |
| Histogram position | "stack" | Clonal/subclonal stacked within each bin |
| Clonal threshold line | VAF = 0.4 | `geom_hline(yintercept=0.4)` · dashed red |
| Violin alpha | 0.38 | Translucent fill · underlying points visible |
| Boxplot width | 0.12 | Narrow overlay on violin |
| Outlier shape | `NA` | Outliers shown as jitter points instead |
| Jitter width | 0.08 | Tight column jitter |
| Y-axis limits | `c(0, 0.88)` | Uniform across violin panels |
| Patchwork layout | `/ (| )` · heights `c(1, 1.2)` | Histogram over two violin panels |
| Figure dimensions | 17 × 12 in | Portrait-landscape hybrid |
| DPI | 180 | Publication quality |

---

## Statistical considerations

- **VAF as a CCF proxy:** In a diploid tumor at 100% purity, VAF ≈ 0.5 for a clonal heterozygous mutation and VAF ≈ 1.0 for homozygous. With tumor purity < 100%, VAF is scaled by purity. A hard threshold at VAF = 0.4 is illustrative; true clonal vs. subclonal classification requires purity and ploidy correction (e.g., via FACETS or PyClone).
- **Bimodal distribution:** In many real tumor datasets, the VAF distribution is bimodal with peaks near 0.5 (clonal) and 0.15–0.25 (subclonal). The synthetic data generates uniform VAF; real data should produce a more structured histogram.
- **Per-gene VAF patterns:** Genes with consistently high VAF (e.g., TP53 near 0.5) are typically early clonal drivers; genes with lower VAF (e.g., PIK3CA near 0.2) may represent late subclonal events or co-driver expansions.

---

## Limitations

- The `CLONAL` flag in the synthetic dataset is randomly assigned independently of VAF value; in real data, clonality should be inferred from VAF adjusted for copy number and purity, not randomly assigned.
- Violins with very few points (small gene/histology groups) can render as thin lines that misrepresent the distribution; a minimum group size check is recommended before rendering violin plots.
- The three-panel layout becomes crowded if more than 6 tumor types or more than 10 genes are added. Consider `facet_wrap` with free y-scales for larger gene panels.

---

## Requirements

```r
library(ggplot2)    # >= 3.4
library(dplyr)      # >= 1.1
library(patchwork)  # >= 1.2
```

---

## Files

| File | Description |
|---|---|
| `12_vaf.R` | Self-contained script · synthetic data generated internally |
| `plots/12_vaf.png` | Output figure · 17 × 12 in · 180 DPI |

---

## References

McGranahan N, Swanton C. Clonal Heterogeneity and Tumor Evolution: Past, Present, and the Future. *Cell.* 2017;168(4):613–628.

Bolli N, et al. Heterogeneity of genomic evolution and mutational profiles in multiple myeloma. *Nat Commun.* 2014;5:2997.

Dentro SC, et al. Characterizing genetic intra-tumor heterogeneity across 2,658 human cancer genomes. *Cell.* 2021;184(8):2239–2254.
