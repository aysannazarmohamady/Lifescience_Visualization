# Plot 1 — Flow Cytometry Hierarchical Gating

A three-panel hierarchical gating strategy visualizing immune cell subsets from FSC/SSC scatter through lineage markers to CD4/CD8 T-cell resolution, derived from patient ADBM baseline biomarker data in the ONCVIZ-001 basket trial dataset.

**Dataset:** ONCVIZ-001 · N = 80 | **Cutoff:** 05 March 2026 | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

Most published flow cytometry figures either show raw FCS data from a single experiment or simulate scatter plots with no biological anchor. This plot bridges the two by deriving gate annotations directly from patient-level ADBM biomarker measurements.

| | Prior art | This work |
|---|---|---|
| Gate annotations | Arbitrary thresholds | CD4⁺/CD8⁺/T-reg labels derived from ADBM baseline means |
| Scatter coordinates | Raw FCS or pure simulation | Simulated from distributions anchored on ADBM proportions |
| Cell type coverage | Typically 1–2 populations | CD4+ T, CD8+ T, B cell, Monocyte, T-reg, Debris |
| Reproducibility | No fixed seed | `set.seed(42)` — bitwise-identical across R ≥ 4.1 |
| Cross-domain traceability | None | FSC/SSC proportions consistent with ADBM CD4=36.9%, CD8=24.1% |

---

## Visual anatomy

### Panel A — FSC vs SSC (All Events)

```
SSC ▲
    |          ┌─────────────┐
    |          │  Lymphocyte │
    |          │    Gate     │
    |__________|_____________|___→ FSC
```

| Element | Description |
|---|---|
| Red rectangle | Lymphocyte gate: FSC 3–7.5, SSC 3.5–7.5 |
| Points | 18 × n_subjects events, colored by simulated cell type |

### Panel B — CD3 vs CD19 (Lymphocyte Gate)

```
CD19 ▲
     |  CD3⁻CD19⁺  │  CD3⁺CD19⁻
     |  (B cells)  │  (T cells)
  ───┼─────────────┼───  CD19 = 3.5
     |     DN      │
     └─────────────┴──→ CD3
                CD3 = 4.2
```

| Element | Description |
|---|---|
| Dashed vertical | CD3 threshold at 4.2 |
| Dashed horizontal | CD19 threshold at 3.5 |
| CD3⁻CD19⁺ quadrant | B cells (salmon) |
| CD3⁺CD19⁻ quadrant | T cells (blue/red) |
| DN label | Double-negative population (grey) |

### Panel C — CD4 vs CD8 (T-cell Gate)

```
CD8 ▲
    │ T-reg (low CD4, low CD8)   │ CD8⁺ T
────┼─────────────────────────────┼──  CD8 = 3.5
    │                             │ CD4⁺ T
    └─────────────────────────────┴──→ CD4
                              CD4 = 3.5
```

| Element | Description |
|---|---|
| CD4⁺ label | Mean % from patient ADBM data (36.9%) |
| CD8⁺ label | Mean % from patient ADBM data (24.1%) |
| T-reg label | Mean % from patient ADBM data (3.0%) |
| Dashed crosshairs | Quadrant boundaries at CD4=3.5, CD8=3.5 |

---

## Data note on cell population coverage

CD4⁺ T, CD8⁺ T, and T-reg annotations are derived directly from patient ADBM biomarker measurements (PARAMCD = CD4, CD8, TREG). FSC/SSC scatter coordinates and CD3/CD19/CD4/CD8 fluorescence values are simulated from Gaussian distributions anchored on these proportions. B cell and Monocyte events are estimated from published PBMC reference frequencies, as direct counts were not measured in the ADBM domain of this dataset.

> For a primary gating figure, replace simulated scatter coordinates with raw FCS files processed via `flowCore::read.FCS()` — the annotation and theme layers are unchanged.

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

### Domains used by this plot

| Domain | Description | Key variables used |
|--------|-------------|---|
| ADSL | Subject-level | `USUBJID`, `ARM`, `TUMORTYPE` |
| ADBM | Biomarker assessments | `PARAMCD` (CD4/CD8/TREG), `AVAL`, `AVISIT` |

### Baseline biomarker summary (ADBM, treatment arm, baseline visit)

| Biomarker | Mean | Range |
|-----------|-----:|------:|
| CD4+ T (%) | 36.9 | 20.0 – 60.0 |
| CD8+ T (%) | 24.1 | 10.0 – 40.0 |
| T-reg (%) | 3.0 | 1.0 – 7.4 |

---

## Output files

| File | Description | Dimensions |
|------|-------------|------------|
| `plot1_flow_gating.R` | R script — Plot 1 only | |
| `plot1_flow_gating.png` | Flow cytometry hierarchical gating | 14×5 in · 150 DPI |

---

## When to use

**Appropriate:**
- Illustrating a hierarchical gating strategy in early-phase oncology translational reports
- Communicating CD4/CD8/T-reg distributions derived from patient biomarker data
- Supplementary figures in publications where full FCS data are unavailable

**Limitations:**
- FSC/SSC coordinates are simulated; not derived from raw cytometer FCS output
- Unsuitable as a primary gating figure without real FCS files
- B cell and Monocyte events are estimated from reference frequencies, not measured

---

## Requirements

```
R >= 4.1
```

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "gridExtra", "ggpubr"))
```

---

## References

Gandhi L, et al. Pembrolizumab plus chemotherapy in metastatic non-small-cell lung cancer. *N Engl J Med.* 2018;378(22):2078–2092.

Weber LM, et al. HDCytoData: Collection of high-dimensional cytometry benchmark datasets in Bioconductor object formats. *F1000Research.* 2019;8:1459.
