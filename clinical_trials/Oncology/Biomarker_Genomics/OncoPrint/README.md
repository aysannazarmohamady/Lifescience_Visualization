# OncoPrint — Somatic Mutation Landscape, Top 12 Genes

A publication-ready OncoPrint heatmap displaying the somatic mutation landscape across 80 patients and 12 frequently altered genes in the ONCVIZ-001 synthetic basket trial, co-displayed with per-patient tumor mutational burden (TMB), best overall response (BOR), and tumor type tracks, plus a per-gene alteration frequency bar panel.

**Dataset:** ONCVIZ-001 · N = 80 patients · ADMUT / ADSL | **Language:** R · ggplot2 · patchwork | **License:** CC BY 4.0

---

## The gap this fills

OncoPrint is the genomic landscape visualization format popularized by cBioPortal and the TCGA consortium. Reproducing it natively in ggplot2 with clinical annotation tracks and a synchronized frequency panel is non-trivial: most R implementations either depend on the `ComplexHeatmap` Bioconductor package (which requires Bioconductor installation and a distinct grammar) or produce figures that lack synchronized clinical tracks.

| | Prior art | This work |
|---|---|---|
| Dependency | ComplexHeatmap (Bioconductor) | Base ggplot2 + patchwork only |
| Clinical annotation tracks | Limited | TMB bar · BOR tile · Tumor type tile |
| Gene frequency panel | Separate call | Right-aligned panel via `patchwork` layout |
| Patient ordering | Alphabetical or arbitrary | Ranked by total alteration count (most-altered first) |
| Multi-hit handling | Split tiles | Single dominant variant per patient per gene (priority hierarchy) |
| Color palette | Default | ONCVIZ-001 palette shared across all 14 figures |

---

## Visual anatomy

```
  TMB ████▌██▌▌▌████▌▌████                        (top bar track)
  BOR ■■■■■■■■■■■■■■■■■■                          (tile: CR/PR/SD/PD)
  TT  ■■■■■■■■■■■■■■■■■■                          (tile: tumor type)
  ─────────────────────────────────────────────────────────────────
  TP53  ■ ■   ■ ■ ■   ■ ■   ■ ■ ■     ■ ■   ■   │ 48% ████████
  KRAS    ■ ■   ■   ■   ■   ■     ■ ■   ■   ■   │ 30% ████
  ...                                             │
  KEAP1 ■     ■   ■                               │ 10% █
  ─────────────────────────────────────────────────────────────────
       Patient columns (n=80, ordered by alteration count) →
```

| Element | Description |
|---|---|
| TMB bar track | `geom_col` per patient · gray fill |
| BOR tile track | `geom_tile` · CR/PR/SD/PD color |
| Tumor type track | `geom_tile` · 5-histology ONCVIZ-001 palette |
| Main heatmap | `geom_tile` · colored by dominant variant class · gray95 = wildtype |
| Right bar panel | `geom_col` per gene · % patients altered · text labels |

---

## Output figure

| File | Dimensions | DPI |
|---|---|---|
| `plots/02_oncoprint.png` | 20 × 10 in | 180 |

---

## Synthetic dataset — ONCVIZ-001 ADaM v1 (ADMUT / ADSL)

```
n patients     80 (TRT = 62 · CTL = 18)
Genes shown    12 (TP53 · KRAS · CDKN2A · PTEN · SMAD4 · PIK3CA ·
                   ARID1A · RB1 · NF1 · MET · BRAF · KEAP1)
Histologies    NSCLC · BRCA · HCC · CRC · PDAC
Seed           2 · fully reproducible
```

### Variables used

| Variable | Source | Description |
|---|---|---|
| `USUBJID` | ADMUT / ADSL | Subject identifier |
| `HUGO_SYMBOL` | ADMUT | Gene name |
| `VARIANT_CLASS` | ADMUT | Mutation type (see palette below) |
| `TMB` | ADSL | Tumor mutational burden (mut/Mb) |
| `BESTRSPC` | ADSL | Best overall response (CR · PR · SD · PD) |
| `TUMORTYPE` | ADSL | Histology (NSCLC · BRCA · HCC · CRC · PDAC) |

### Variant priority hierarchy (multi-hit tie-breaking)

When a patient carries more than one variant in the same gene, one is selected for display using the following priority (lower = higher priority):

| Priority | Variant class |
|---|---|
| 0 | Amplification |
| 1 | Nonsense_Mutation |
| 2 | Frame_Shift_Del / Frame_Shift_Ins |
| 3 | Splice_Site |
| 4 | In_Frame_Del |
| 5 | Missense_Mutation |
| 6 | Other_Mutation |

---

## Key parameters

| Parameter | Value | Description |
|---|---|---|
| Genes displayed | 12 | Fixed list; ranked by pan-cohort frequency in source script |
| Patient ordering | By alteration count (desc) | Most-altered patients at left |
| Layout widths | `c(19, 1)` | 19 parts heatmap · 1 part frequency bar |
| Track heights | `c(TMB, BOR, TT, heat)` | Assembled with `patchwork` `/` operator |
| Color: Missense | `#2c5f8a` | Blue |
| Color: Nonsense | `#c0392b` | Red |
| Color: Frame_Shift_Del | `#e67e22` | Orange |
| Color: Frame_Shift_Ins | `#f1c40f` | Yellow |
| Color: Splice_Site | `#8e44ad` | Purple |
| Color: In_Frame_Del | `#27ae60` | Green |
| Color: Amplification | `#7b0000` | Dark red |
| Color: CR | `#1a3f8f` | Dark blue |
| Color: PR | `#4f8fd4` | Light blue |
| Color: SD | `#e8a020` | Amber |
| Color: PD | `#c0392b` | Red |
| Figure dimensions | 20 × 10 in | Wide landscape for 80-patient columns |
| DPI | 180 | Publication quality |

---

## Statistical considerations

- **Patient ordering:** Ordering by total alteration count (rather than by gene or hierarchical clustering) maximizes the visual impression of mutual exclusivity and co-occurrence patterns at the left margin. Hierarchical clustering of the binary alteration matrix is an alternative that reveals biological co-occurrence structure more rigorously.
- **Frequency calculation:** Gene alteration frequency = `n_distinct(USUBJID with ≥1 mutation in gene) / 80`. This counts each patient once regardless of the number of mutations in that gene.
- **Wildtype cells:** Patients with no detected alteration in a gene receive `na.value = "grey95"`. This represents the true wildtype state only if sequencing coverage is adequate; in practice, missing data should be distinguished from confirmed wildtype.

---

## Limitations

- Multi-hit display (split tiles for patients with two variant classes in one gene) is not implemented; the dominant variant by priority hierarchy is shown.
- The `patchwork` layout does not enforce pixel-aligned column widths between the TMB/BOR/TT tracks and the main heatmap; minor misalignment may occur at very large n.
- Mutual exclusivity and co-occurrence statistics (e.g., Fisher's exact test per gene pair) are not annotated; these can be computed via `maftools::somaticInteractions()` or a custom pairwise Fisher's test loop.

---

## Requirements

```r
library(ggplot2)    # >= 3.4
library(dplyr)      # >= 1.1
library(tidyr)      # >= 1.3
library(patchwork)  # >= 1.2
```

---

## Files

| File | Description |
|---|---|
| `02_oncoprint.R` | Self-contained script · synthetic data generated internally |
| `plots/02_oncoprint.png` | Output figure · 20 × 10 in · 180 DPI |

---

## References

Cerami E, et al. The cBio Cancer Genomics Portal: An Open Platform for Exploring Multidimensional Cancer Genomics Data. *Cancer Discov.* 2012;2(5):401–404.

Gu Z, Eils R, Schlesner M. Complex heatmaps reveal patterns and correlations in multidimensional genomic data. *Bioinformatics.* 2016;32(18):2847–2849.

Alexandrov LB, et al. Signatures of mutational processes in human cancer. *Nature.* 2013;500(7463):415–421.

Pedersen TL. patchwork: The Composer of Plots. R package version 1.2.0. https://CRAN.R-project.org/package=patchwork
