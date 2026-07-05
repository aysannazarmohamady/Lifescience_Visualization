# CyTOF Dot Plot

A marker expression dot plot displaying 18 CyTOF markers across 8 immune cell populations, where dot size encodes percentage of cells expressing each marker and dot color encodes mean normalized expression (0–1), derived from canonical literature profiles anchored on the ONCVIZ-001 basket trial immune landscape.

**Dataset:** ONCVIZ-001 · N = 80 | **Cutoff:** 05 March 2026 | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

CyTOF dot plots in the literature are typically shown for a single experiment with narrow marker panels. This plot provides a comprehensive 18-marker × 8-cell-type reference panel covering lineage, T-reg identity, cytotoxicity, exhaustion/checkpoint, cytokine, and proliferation markers in a single figure consistent with the ONCVIZ-001 immune landscape.

| | Prior art | This work |
|---|---|---|
| Marker coverage | 5–10 markers per panel | 18 markers across 6 functional groups |
| Cell type coverage | 3–5 populations | 8 populations: CD4+ T, CD8+ T, NK, T-reg, B cell, Monocyte, pDC, Neutro |
| Functional groups | Lineage only | Lineage, T-reg identity, cytotoxicity, exhaustion/checkpoint, cytokines, proliferation |
| Cross-domain consistency | Standalone figure | Consistent cell type palette with Plots 1–5 |
| Reproducibility | No fixed seed | `set.seed(42)` — bitwise-identical across R ≥ 4.1 |
| Calibration traceability | Unspecified | Expression profiles anchored to Bodenmiller et al. 2012 and published flow/CyTOF literature |

---

## Visual anatomy

### Dot plot layout

```
Marker     CD4+T  CD8+T  NK  T-reg  B cell  Mono  pDC  Neutro
CD3          ●      ●     ○     ●      ○       ○     ○     ○
CD4          ●      ○     ○     ●      ○       ○     ○     ○
CD8a         ○      ●     ○     ○      ○       ○     ○     ○
CD56         ○      ○     ●     ○      ○       ○     ○     ○
CD19         ○      ○     ○     ○      ●       ○     ○     ○
...
Ki-67        ○      ○     ○     ○      ○       ○     ○     ○

● large/dark = high expression · ○ small/light = low expression
```

| Element | Description |
|---|---|
| Dot size | % cells expressing (range 1–12 pt, scale 0–100%) |
| Dot fill | Mean normalized expression 0–1 (white → yellow → orange → red → dark red) |
| Dot outline | `grey60`, stroke = 0.3, shape = 21 |
| Size breaks | 25%, 50%, 75%, 100% |
| Color gradient | white → `#FEE08B` → `#FC8D59` → `#D73027` → `#67001F` |

### Marker groups

| Group | Markers |
|---|---|
| Lineage | CD3, CD4, CD8a, CD56, CD19, CD14, CD16 |
| T-reg identity | CD25, CD127, FoxP3 |
| Cytotoxicity | GZMB |
| Exhaustion / checkpoint | PD-1, TIM-3, LAG-3 |
| Cytokines | IFN-γ, TNF-α, IL-2 |
| Proliferation | Ki-67 |

### Key expression highlights (mean expression values)

| Marker | Highest-expressing cell type | Value |
|--------|------------------------------|-------|
| CD3 | CD4+ T, CD8+ T, T-reg | 0.95 |
| CD4 | CD4+ T | 0.85 |
| CD8a | CD8+ T | 0.90 |
| CD56 | NK | 0.90 |
| CD19 | B cell | 0.95 |
| CD14 | Monocyte | 0.90 |
| FoxP3 | T-reg | 0.92 |
| GZMB | CD8+ T | 0.65 |
| PD-1 | CD8+ T | 0.45 |
| IFN-γ | CD8+ T | 0.60 |

---

## Data note on cell population coverage

All 18 marker expression profiles are canonical literature-derived values anchored on Bodenmiller et al. 2012 and standard flow/CyTOF reference panels. No patient-derived single-cell data were used. Percentage-expressing values are computed as a deterministic transform of mean expression (`MeanExpr × 110 + runif(0, 8)`, capped at 100%) to simulate realistic dot size variation. This approach is consistent with methods used in biomarker-exploratory analyses when patient-level CyTOF panels are unavailable.

> For fully data-driven CyTOF dot plots, replace `expr_profiles` with per-cell-type expression means computed from `HDCytoData::Bodenmiller_BCR_XL_SE()` — the `ggplot2` dot plot layer is unchanged.

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
| — | No ADaM domain read directly | Cell type palette (`CELL_COLS`) shared with all plots |

### Baseline biomarker summary (ADBM, treatment arm — for context)

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
| `plot6_cytof_dotplot.R` | R script — Plot 6 only | |
| `plot6_cytof_dotplot.png` | CyTOF marker expression dot plot | 11×9 in · 150 DPI |

---

## When to use

**Appropriate:**
- Illustrating canonical CyTOF marker profiles across immune populations in translational reports
- Supplementary reference panels showing lineage and functional marker co-expression
- Communicating checkpoint and exhaustion marker distribution alongside clinical biomarker data

**Limitations:**
- All expression profiles are canonical literature-derived; not from patient-derived single-cell data
- % expressing is a deterministic transform of mean expression with small random jitter — not empirical
- Does not capture patient-to-patient variability in marker expression

---

## Requirements

```
R >= 4.1
```

```r
install.packages(c("ggplot2", "scales"))
```

---

## References

Bodenmiller B, et al. Multiplexed mass cytometry profiling of cellular states perturbed by small-molecule regulators. *Nat Biotechnol.* 2012;30(9):858–867.

Weber LM, et al. HDCytoData: Collection of high-dimensional cytometry benchmark datasets in Bioconductor object formats. *F1000Research.* 2019;8:1459.
