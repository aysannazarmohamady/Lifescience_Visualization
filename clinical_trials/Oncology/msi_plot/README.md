# MSI Plot — Microsatellite Instability Score Across Cohort

A three-panel publication-ready figure characterizing microsatellite instability (MSI) status across the ONCVIZ-001 cohort: (top) MSI score ranked bar chart with MSI-High threshold annotation; (bottom-left) per-tumor-type violin + boxplot + jitter; (bottom-right) stacked horizontal bar showing MSI-H vs. MSS proportion by histology.

**Dataset:** ONCVIZ-001 · N = 80 patients · ADSL | **Language:** R · ggplot2 · patchwork | **License:** CC BY 4.0

---

## The gap this fills

MSI status is a predictive biomarker for immunotherapy response (pembrolizumab pan-tumor approval, 2017) and is routinely reported in basket trial publications. Standard reporting shows only a pie chart or a table. A ranked bar chart combined with tumor-type distributions and proportion panels provides a more complete characterization of MSI heterogeneity across histologies.

| | Prior art | This work |
|---|---|---|
| MSI visualization | Pie chart or binary table | Ranked continuous score bar chart |
| Threshold annotation | Binary cutoff mentioned in text | Annotated dashed line with MSI-H patient count |
| Histology breakdown | Table footnote | Bottom-left violin + boxplot + jitter |
| Proportion by histology | Not shown | Bottom-right stacked bar with percentage labels |
| Score type | Binary MSI-H / MSS | Continuous MSI score + binary classification |
| Reproducibility | No fixed seed | `set.seed(14)` · bitwise-identical |

---

## Visual anatomy

```
  [TOP PANEL — Ranked bar chart]
  MSI  4.4 ─────────────────────────────────────────
  Score     ██ ██ (4 MSI-H patients above threshold)
       3.63 ──────── MSI-H threshold ─────── (red dashed)
            ████████████████████████████████ (76 MSS patients)
       0.0  ─────────────────────────────────────────
            Patients ranked by MSI score (n=80) →

  [BOTTOM-LEFT]              [BOTTOM-RIGHT]
  Violin by tumor type:      Stacked proportion by histology:
  NSCLC BRCA HCC CRC PDAC   NSCLC ████████████████████ 0% MSI-H
  Points colored by status   BRCA  ████████████████████
                             HCC   ████████████████████
                             CRC   ██████████████████░░ ~8% MSI-H
                             PDAC  ████████████████████
```

| Element | Description |
|---|---|
| Top ranked bars | `geom_col` · fill by MSI status · ranked descending |
| Dashed threshold line | `geom_hline(yintercept = 3.63)` · red dashed |
| Threshold annotation | MSI-H count and threshold value via `annotate("text")` |
| Bottom-left violin | `geom_violin + geom_boxplot + geom_jitter` · fill by tumor type |
| Bottom-left jitter | Points colored by MSI status (red = MSI-H · blue = MSS) |
| Bottom-right bars | `geom_col(position="stack")` · MSS (blue) + MSI-H (red) |
| Bar labels | White bold percentage via `position_stack(vjust=0.5)` |
| Patchwork layout | `/ (| )` · heights `c(1.2, 1)` |

---

## Output figure

| File | Dimensions | DPI |
|---|---|---|
| `plots/14_msi.png` | 16 × 10 in | 180 |

---

## Synthetic dataset — ONCVIZ-001 ADaM v1 (ADSL)

```
n patients     80 (NSCLC=25 · BRCA=13 · HCC=12 · CRC=16 · PDAC=14)
MSI score      Continuous · MSS range: [0.5, 3.5] · MSI-H range: [3.7, 4.4]
MSI-H          4 patients (5.0%) · threshold = 3.63
Seed           14 · fully reproducible
```

### MSI status classification

| Status | Abbreviation | Color | Definition |
|---|---|---|---|
| MSI-High | MSI-H | `#c0392b` (Red) | MSI score ≥ 3.63 |
| Microsatellite Stable | MSS | `#5b9bd5` (Blue) | MSI score < 3.63 |

### MSI score generation (synthetic)

```r
MSI_SCORE <- c(runif(76, 0.5, 3.5),   # MSS patients
               runif(4, 3.7, 4.4))     # MSI-H patients
MSISTS    <- ifelse(MSI_SCORE >= 3.63, "MSI-H", "MSS")
```

### Variables used

| Variable | Source | Description |
|---|---|---|
| `USUBJID` | ADSL | Subject identifier |
| `MSI_SCORE` | ADSL | Continuous MSI score |
| `MSISTS` | ADSL | Binary MSI status ("MSI-H" · "MSS") |
| `TUMORTYPE` | ADSL | Histology |

---

## Key parameters

| Parameter | Value | Description |
|---|---|---|
| MSI-H threshold | 3.63 | `geom_hline(yintercept = 3.63)` · from ADSL |
| Top bar width | 0.92 | Dense fill · minimal white gaps between bars |
| Top bar alpha | 0.88 | Slight transparency |
| Violin alpha | 0.38 | Translucent · reveals jitter points |
| Boxplot width | 0.12 | Narrow overlay |
| Jitter width | 0.10 | Tight column |
| Stacked bar labels | White bold · `position_stack(vjust=0.5)` | Centered in each segment |
| Y-axis limits (violin) | `c(0, 4.5)` | Uniform across tumor types |
| Patchwork layout | `/ (| )` · heights `c(1.2, 1)` | Ranked bar taller than bottom panels |
| Figure dimensions | 16 × 10 in | Portrait-landscape hybrid |
| DPI | 180 | Publication quality |

---

## Statistical considerations

- **MSI prevalence (5%):** With only 4 MSI-H patients (5.0% of 80), statistical comparisons between MSI-H and MSS groups (e.g., OS by MSI status) have very low power. This is representative of pan-tumor basket trials where MSI-H prevalence varies from ~1% (ovarian) to ~20% (endometrial).
- **Threshold selection:** The MSI-H threshold of 3.63 is derived from the ADSL dataset's data cutoff distribution and does not correspond to a validated assay-specific threshold. Real MSI scoring uses fragment length analysis (MSI PCR), next-generation sequencing (MSISensor, MANTIS), or immunohistochemistry (MMR protein loss) with assay-specific cutoffs.
- **Continuous score reporting:** Reporting the continuous MSI score alongside binary classification is recommended by CAP (College of American Pathologists) to allow laboratories to adjust thresholds as assay-specific data accumulates.
- **Histology-specific rates:** In the right panel, tumor types with 0 MSI-H patients will show 100% blue bars. With n = 12–25 per histology and 4 MSI-H events distributed, the observed rates have very wide confidence intervals.

---

## Limitations

- The threshold at 3.63 is arbitrary in the synthetic dataset. For real data, the threshold must be clinically validated and will vary by assay (MSISensor, MANTIS, MSI PCR) and tumor type.
- With 4 MSI-H patients, the violin and jitter in the bottom-left panel for any tumor type with ≤1 MSI-H patient will be visually dominated by the MSS distribution; the MSI-H points are visible as isolated red dots.
- The stacked proportion bar for any histology with 0 MSI-H will be 100% blue, which is visually uninformative. A minimum MSI-H n label or "n = 0" annotation is recommended.

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
| `14_msi.R` | Self-contained script · synthetic data generated internally |
| `plots/14_msi.png` | Output figure · 16 × 10 in · 180 DPI |

---

## References

Le DT, et al. PD-1 Blockade in Tumors with Mismatch-Repair Deficiency. *N Engl J Med.* 2015;372(26):2509–2520.

Marabelle A, et al. Efficacy of Pembrolizumab in Patients with Noncolorectal High Microsatellite Instability/Mismatch Repair–Deficient Cancer: Results From the Phase II KEYNOTE-158 Study. *J Clin Oncol.* 2020;38(1):1–10.

Bonneville R, et al. Landscape of Microsatellite Instability Across 39 Cancer Types. *JCO Precis Oncol.* 2017;1:PO.17.00073.

Salipante SJ, et al. Microsatellite instability detection by next generation sequencing. *Clin Chem.* 2014;60(9):1192–1199.

Kautto EA, et al. Performance evaluation for rapid detection of pan-cancer microsatellite instability with MANTIS. *Oncotarget.* 2017;8(5):7452–7463.
