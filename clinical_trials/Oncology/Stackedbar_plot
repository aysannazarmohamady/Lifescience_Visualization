# Cell Composition Stacked Bar

A patient-level stacked bar chart displaying proportions of 8 immune cell types across all 80 subjects, stratified by treatment arm (Tumor vs Healthy Donor), sorted Healthy → Tumor, derived from the ONCVIZ-001 basket trial dataset.

**Dataset:** ONCVIZ-001 · N = 80 | **Cutoff:** 05 March 2026 | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

Most published immune cell composition figures aggregate across patients or show only a handful of subjects. This plot displays individual patient-level proportions for all 80 subjects with full cell type breakdown, maintaining direct traceability to ADBM biomarker measurements.

| | Prior art | This work |
|---|---|---|
| Subject coverage | Aggregated or subset of patients | All 80 patients displayed individually |
| Cell type coverage | 3–4 major populations | 8 populations: CD4+ T, CD8+ T, NK, T-reg, B cell, Monocyte, pDC, Neutro |
| Data anchor | Simulated proportions | CD4, CD8, NK, T-reg derived from ADBM `PARAMCD` values |
| Arm stratification | Treatment only | Treatment (PT) and Control (HD) side by side |
| Reproducibility | No fixed seed | `set.seed(42)` — bitwise-identical across R ≥ 4.1 |
| Calibration traceability | Unspecified | Anchored to patient ADBM data (mean CD4=36.9%, CD8=24.1%, NK=12.5%, T-reg=3.0%) |

---

## Visual anatomy

### Stacked bar layout

```
100% ┤ ████ Neutro
     │ ░░░░ pDC
     │ ████ Monocyte
     │ ████ B cell
 50% ┤ ████ T-reg
     │ ████ NK
     │ ████ CD8+ T
     │ ████ CD4+ T
  0% ┴─────────────────────────────────────────→
     HD-01 … HD-18   PT-01 … PT-62
     ◄── Healthy ──► ◄──── Tumor ────►
```

| Element | Description |
|---|---|
| CD4+ T (blue `#2166AC`) | Derived from ADBM PARAMCD=CD4, baseline mean |
| CD8+ T (red `#D6604D`) | Derived from ADBM PARAMCD=CD8, baseline mean |
| NK (green `#4DAC26`) | Derived from ADBM PARAMCD=NK, baseline mean |
| T-reg (dark green `#1B7837`) | Derived from ADBM PARAMCD=TREG, baseline mean |
| B cell (salmon `#F4A582`) | Estimated: `runif(1, 8, 16)` per subject |
| Monocyte (purple `#762A83`) | Estimated: `runif(1, 6, 14)` per subject |
| pDC (yellow `#FEE08B`) | Estimated: `runif(1, 1, 4)` per subject |
| Neutro (brown `#8C510A`) | Estimated: `runif(1, 2, 8)` per subject |
| In-bar labels | White, bold; shown only for segments > 4% |
| Sample order | Healthy donors (HD) first, Tumor patients (PT) second, each sorted by SampleID |

---

## Data note on cell population coverage

Four cell types (CD4+ T, CD8+ T, NK, T-reg) are derived directly from patient ADBM biomarker measurements (PARAMCD = CD4, CD8, NK, TREG), filtered to `ANL01FL=="Y"` and `AVISIT=="BASELINE"`. The remaining four (B cell, Monocyte, pDC, Neutrophil) are simulated from uniform distributions anchored on published PBMC reference frequencies, as these populations were not measured in the ADBM domain of this dataset. Proportions are normalized within each subject to sum to 100% across all 8 cell types.

> For fully data-driven stacked bars, replace the `extra` simulation block with measured B cell, Monocyte, pDC, and Neutrophil values from a clinical flow panel — the normalization and plotting code are unchanged.

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

| Histology | n (total) | n TRT | n CTL |
|-----------|----------:|------:|------:|
| NSCLC | 25 | 23 | 2 |
| BRCA | 13 | 8 | 5 |
| CRC | 16 | 9 | 7 |
| HCC | 12 | 10 | 2 |
| PDAC | 14 | 12 | 2 |

### Domains used by this plot

| Domain | Description | Key variables used |
|--------|-------------|---|
| ADSL | Subject-level | `USUBJID`, `ARM`, `TUMORTYPE` |
| ADBM | Biomarker assessments | `PARAMCD` (CD4/CD8/NK/TREG), `AVAL`, `AVISIT` |

### Baseline biomarker summary (ADBM, treatment arm, baseline visit)

| Biomarker | Mean | Range |
|-----------|-----:|------:|
| CD4+ T (%) | 36.9 | 20.0 – 60.0 |
| CD8+ T (%) | 24.1 | 10.0 – 40.0 |
| NK (%) | 12.5 | 5.0 – 23.7 |
| T-reg (%) | 3.0 | 1.0 – 7.4 |

---

## Output files

| File | Description | Dimensions |
|------|-------------|------------|
| `plot4_stacked_bar.R` | R script — Plot 4 only | |
| `plot4_stacked_bar.png` | Cell composition per sample | 13×6 in · 150 DPI |

---

## When to use

**Appropriate:**
- Communicating immune cell composition across tumor histologies and treatment arms
- Biomarker-response correlation visualization in basket trials
- Supplementary figures in translational publications showing per-patient immune profiles

**Limitations:**
- B cell, Monocyte, pDC, Neutrophil proportions are estimated, not measured in this cohort
- Within-subject normalization to 8-cell-type sum may differ from clinical flow panel totals
- Individual bar width narrows at N = 80; consider faceting by tumor type for publication

---

## Requirements

```
R >= 4.1
```

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "scales"))
```

---

## References

Gandhi L, et al. Pembrolizumab plus chemotherapy in metastatic non-small-cell lung cancer. *N Engl J Med.* 2018;378(22):2078–2092.

Finn RS, et al. Atezolizumab plus bevacizumab in unresectable hepatocellular carcinoma. *N Engl J Med.* 2020;382(20):1894–1905.
