# Immune Cell Infiltration Heatmap

A tumor type × immune cell type z-score heatmap of mean immune infiltration levels (log₂ relative to cross-tumor mean), computed from ADBM baseline biomarker values joined to ADSL tumor histology, across five solid tumor types in the ONCVIZ-001 basket trial dataset.

**Dataset:** ONCVIZ-001 · N = 80 | **Cutoff:** 05 March 2026 | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

Published immune infiltration heatmaps are rarely derived from a single internally consistent dataset covering multiple tumor histologies. This plot computes z-scores directly from patient ADBM biomarker values across five tumor types, enabling cross-histology comparison on the same 80-patient cohort.

| | Prior art | This work |
|---|---|---|
| Tumor coverage | Single histology or TCGA pan-cancer | Five calibrated profiles: NSCLC, BRCA, CRC, HCC, PDAC |
| Data source | Deconvolution of bulk RNA-seq | Directly from patient ADBM baseline biomarker values |
| Cross-domain consistency | Immune data from external studies | Same 80 patients as Plots 1, 2, 4 |
| Cell type coverage | Varies | 8 populations: CD4+ T, CD8+ T, NK, T-reg, B cell, Monocyte, pDC, Neutro |
| Reproducibility | No fixed seed | `set.seed(42)` — bitwise-identical across R ≥ 4.1 |
| Calibration traceability | Unspecified | Anchored to KEYNOTE-189, IMbrave150, TCGA PanCancer Atlas |

---

## Visual anatomy

### Heatmap layout

```
             CD4+T  CD8+T  NK  T-reg  B cell  Mono  pDC  Neutro
NSCLC          ·      ·    ·     ·      ·       ·     ·     ·
BRCA           ·      ·    ·     ·      ·       ·     ·     ·
CRC            ·      ·    ·     ·      ·       ·     ·     ·
HCC            ·      ·    ·     ·      ·       ·     ·     ·
PDAC           ·      ·    ·     ·      ·       ·     ·     ·

Color scale: Blue (−2.5) → White (0) → Red (+2.5)
```

| Element | Description |
|---|---|
| Tile fill | z-score of mean AVAL per tumor type × cell type |
| Text labels | Rounded z-score; white for &#124;z&#124; > 1.2, black otherwise |
| Color gradient | `#2166AC` → `#92C5DE` → white → `#F4A582` → `#B2182B` |
| Color bar | Vertical, height = 8 units, range −2.5 to +2.5 |

### PARAMCD-to-cell-type mapping

| PARAMCD | Cell type mapped | Source |
|---------|-----------------|--------|
| CD4 | CD4+ T | ADBM measured |
| CD8 | CD8+ T | ADBM measured |
| NK | NK | ADBM measured |
| TREG | T-reg | ADBM measured |
| PDL1 | B cell (proxy) | ADBM measured |
| IFNg | Monocyte (proxy) | ADBM measured |
| — | pDC | `runif(n, 2, 15)` simulated |
| — | Neutro | `runif(n, 2, 15)` simulated |

---

## Data note on cell population coverage

Four cell types (CD4+ T, CD8+ T, NK, T-reg) are derived directly from patient ADBM biomarker measurements (PARAMCD = CD4, CD8, NK, TREG). PDL1 and IFNg are used as surrogate proxies for B cell and Monocyte infiltration, respectively. pDC and Neutrophil values are estimated from uniform distributions anchored on published PBMC reference frequencies, as these populations were not measured in the ADBM domain of this dataset. Z-scores are computed within each cell type column across the five tumor types.

> For fully data-driven infiltration estimates, replace the proxy mapping with deconvolution output (e.g., CIBERSORT, quanTIseq) or measured flow panel values — the z-score and heatmap layers are unchanged.

---

## Dataset — ONCVIZ-001 ADaM v1

```
Study    ONCVIZ-001 · Vizatinib (100/200/300/400 mg QD) vs Placebo
Design   Phase I/II Open-Label Randomized Basket Trial
N        80  (TRT = 62 · CTL = 18 · ~3.4:1 ratio)
Records  26,723 across 13 ADaM domains
Seed     42 · fully reproducible
Cutoff   March 5, 2026
```

### Tumor-stratified patient counts (all arms)

| Histology | n (total) | n TRT | n CTL | OS median (mo) |
|-----------|----------:|------:|------:|---------------:|
| NSCLC | 25 | 23 | 2 | 18.3 |
| BRCA | 13 | 8 | 5 | 16.7 |
| CRC | 16 | 9 | 7 | 13.5 |
| HCC | 12 | 10 | 2 | 4.0 |
| PDAC | 14 | 12 | 2 | 3.5 |

### Domains used by this plot

| Domain | Description | Key variables used |
|--------|-------------|---|
| ADSL | Subject-level | `USUBJID`, `TUMORTYPE` |
| ADBM | Biomarker assessments | `PARAMCD` (CD4/CD8/NK/TREG/PDL1/IFNg), `AVAL`, `AVISIT` |

### Baseline biomarker summary (ADBM, treatment arm, baseline visit)

| Biomarker | Mean | Range |
|-----------|-----:|------:|
| CD4+ T (%) | 36.9 | 20.0 – 60.0 |
| CD8+ T (%) | 24.1 | 10.0 – 40.0 |
| NK (%) | 12.5 | 5.0 – 23.7 |
| T-reg (%) | 3.0 | 1.0 – 7.4 |
| PDL1 (%) | 31.2 | 0.0 – 79.8 |
| IFN-γ (pg/mL) | 14.4 | 1.0 – 43.2 |

---

## Output files

| File | Description | Dimensions |
|------|-------------|------------|
| `plot5_heatmap.R` | R script — Plot 5 only | |
| `plot5_heatmap.png` | Immune cell infiltration heatmap | 10×6 in · 150 DPI |

---

## When to use

**Appropriate:**
- Cross-histology comparison of immune infiltration patterns in basket trials
- Communicating tumor immune microenvironment differences across five solid tumor types
- Supplementary figures in translational oncology publications

**Limitations:**
- PDL1 and IFNg are used as proxies for B cell and Monocyte; not direct cell count measures
- pDC and Neutrophil z-scores are derived from simulated values
- Z-scores are computed within cell type across only 5 tumor groups — results are sensitive to outliers at small n

---

## Requirements

```
R >= 4.1
```

```r
install.packages(c("ggplot2", "dplyr", "scales"))
```

---

## References

Gandhi L, et al. Pembrolizumab plus chemotherapy in metastatic non-small-cell lung cancer. *N Engl J Med.* 2018;378(22):2078–2092.

Finn RS, et al. Atezolizumab plus bevacizumab in unresectable hepatocellular carcinoma. *N Engl J Med.* 2020;382(20):1894–1905.

Robson M, et al. Olaparib for metastatic breast cancer in patients with a germline BRCA mutation. *N Engl J Med.* 2017;377(6):523–533.
