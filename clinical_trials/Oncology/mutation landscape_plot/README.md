# Mutation Landscape Plot — Somatic Alteration Frequencies Across Histologies

A dual-panel publication-ready figure displaying somatic point mutation prevalence for the 15 most frequently altered genes across five tumor histologies (left panel, grouped bar chart) alongside pan-cohort alteration frequency (right panel), with genes ranked by pan-cohort frequency. Produced from the ONCVIZ-001 synthetic basket trial dataset.

**Dataset:** ONCVIZ-001 · N = 80 patients · ADMUT / ADSL · 15 genes | **Language:** R · ggplot2 · patchwork | **License:** CC BY 4.0

---

## The gap this fills

Mutation landscape bar charts (also called "lollipop frequency charts" or "alteration frequency plots") summarize how often each driver gene is mutated across a cohort and its histological subgroups. Published examples typically show either a single stacked bar per gene (collapsing histology) or a faceted per-histology panel (losing the pan-cohort summary). Combining both in one synchronized dual-panel layout with shared y-axis (gene) and a right-aligned pan-cohort bar requires the `patchwork` compositor and careful factor-level alignment.

| | Prior art | This work |
|---|---|---|
| Histology breakdown | Separate facet or collapsed | Side-by-side per-histology bars within each gene row |
| Pan-cohort summary | Separate figure | Right-aligned panel via `patchwork`, same gene order |
| Gene ranking | Alphabetical | Ascending pan-cohort frequency (most frequent at top) |
| Denominator | Absolute count | % of patients within each histology's cohort size |
| Color palette | Default ggplot2 | ONCVIZ-001 tumor-type palette (shared across all 14 figures) |
| Reproducibility | No fixed seed | `set.seed(5)` · bitwise-identical |

---

## Visual anatomy

```
  Gene (ranked by pan-cohort %) ← less common              Pan-cohort
  STK11  █  ░░                                             │  8% ██
  RET    ██ ░░ ░                                           │ 10% ██
  EGFR   ████ ░ ░░                                         │ 14% ███
  ...                                                      │
  TP53   ████████  ████████  ██████                        │ 48% ██████
         NSCLC BRCA HCC CRC PDAC     ← dodge bars         │
                Mutation Prevalence Within Tumor Type (%)  │ Overall n=80
```

| Element | Description |
|---|---|
| Grouped bar (left) | Per-histology mutation prevalence · `position_dodge(width=0.8)` |
| Bar colors | 5-histology ONCVIZ-001 palette |
| Right panel bars | Pan-cohort % · gray fill (`#444444`) |
| Right panel text | Percentage label · `hjust = -0.1` |
| Gene order | Ascending pan-cohort frequency (least common at bottom) |
| Patchwork layout | `widths = c(4, 1)` |

---

## Output figure

| File | Dimensions | DPI |
|---|---|---|
| `plots/05_mutation_landscape.png` | 15 × 10 in | 180 |

---

## Synthetic dataset — ONCVIZ-001 ADaM v1 (ADMUT / ADSL)

```
n patients     80 (NSCLC=25 · BRCA=13 · HCC=12 · CRC=16 · PDAC=14)
Genes shown    15: TP53 · KRAS · PTEN · ARID1A · CDKN2A · PIK3CA ·
                   KEAP1 · SMAD4 · RB1 · NF1 · MET · EGFR · BRAF · RET · STK11
Variant type   Somatic point mutations and indels (VARIANT_TYPE ≠ "CNV")
Seed           5 · fully reproducible
```

### Variables used

| Variable | Source | Description |
|---|---|---|
| `USUBJID` | ADMUT / ADSL | Subject identifier |
| `HUGO_SYMBOL` | ADMUT | Gene name |
| `TUMORTYPE` | ADSL | Histology |
| `ANL01FL` | ADMUT | Analysis flag (filtered to "Y") |

### Denominator per histology

Mutation prevalence is calculated as `n_distinct(patients with ≥1 mutation in gene) / N_histology`:

| Histology | N |
|---|---|
| NSCLC | 25 |
| BRCA | 13 |
| HCC | 12 |
| CRC | 16 |
| PDAC | 14 |

---

## Key parameters

| Parameter | Value | Description |
|---|---|---|
| Genes displayed | 15 | Fixed list from source ADMUT |
| Gene ordering | Ascending pan-cohort frequency | `arrange(pan) %>% pull(HUGO_SYMBOL)` |
| Bar position | `position_dodge(0.8)` | Side-by-side per histology |
| Bar width | 0.7 | Within dodge slot |
| Bar alpha | 0.88 | Slight transparency |
| X-axis left | `limits = c(0, 110)` | Headroom for bars near 100% |
| X-axis right | `limits = c(0, 65)` | Pan-cohort panel |
| Layout widths | `c(4, 1)` | 4:1 main vs. summary panel |
| Color: NSCLC | `#2c5f8a` | Blue |
| Color: BRCA | `#c0392b` | Red |
| Color: HCC | `#27ae60` | Green |
| Color: CRC | `#e67e22` | Orange |
| Color: PDAC | `#8e44ad` | Purple |
| Figure dimensions | 15 × 10 in | Portrait-landscape hybrid |
| DPI | 180 | Publication quality |

---

## Statistical considerations

- **Prevalence vs. frequency:** Bars show the proportion of *patients* with ≥1 mutation in a gene within a histology, not mutation count. A patient with 3 TP53 mutations is counted once.
- **Histology-specific denominator:** Each bar is normalized to its own histology's cohort size. Absolute count bars (using the full n = 80 as denominator for all histologies) would visually overweight larger histologies (NSCLC).
- **No significance test:** The chart is descriptive. To test whether mutation frequency differs between histologies for a given gene, a Fisher's exact test or chi-squared test across the 5 × 2 contingency table is appropriate.

---

## Limitations

- Variant type filtering (point mutations / indels only, no CNV) means amplifications and deletions in these genes are not reflected.
- With histology-specific denominators as small as n = 12 (HCC), a single patient's mutation changes the bar by ~8.3%; interpret small-histology bars with caution.
- Gene list is fixed at 15 and would need updating if a different frequency ranking is desired from real data.

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
| `05_mutation_landscape.R` | Self-contained script · synthetic data generated internally |
| `plots/05_mutation_landscape.png` | Output figure · 15 × 10 in · 180 DPI |

---

## References

Bailey MH, et al. Comprehensive Characterization of Cancer Driver Genes and Mutations. *Cell.* 2018;173(2):371–385.

Cancer Genome Atlas Research Network. Comprehensive genomic characterization of squamous cell lung cancers. *Nature.* 2012;489(7417):519–525.

Sanchez-Vega F, et al. Oncogenic Signaling Pathways in The Cancer Genome Atlas. *Cell.* 2018;173(2):321–337.
