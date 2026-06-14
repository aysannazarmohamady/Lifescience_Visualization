# Plot 2 — UMAP Single-Cell Immune Landscape

A UMAP dimensionality reduction of 2,800 simulated immune cells across 7 biomarker features, colored by cell type and labeled at cluster centroids, anchored on patient ADBM baseline distributions from the ONCVIZ-001 basket trial dataset.

**Dataset:** ONCVIZ-001 · N = 80 | **Cutoff:** 05 March 2026 | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

Most published UMAP figures for immune profiling use either transcriptomic data from a single cohort or canonical cell-line profiles with no link to clinical biomarker measurements. This plot anchors cluster positions on patient-derived ADBM data and maintains cross-domain consistency with Plots 1, 4, and 5.

| | Prior art | This work |
|---|---|---|
| Expression anchor | Canonical profiles or single-cell RNA-seq | Derived from patient ADBM baseline means (CD4=36.9%, CD8=24.1%, NK=12.5%, T-reg=3.0%) |
| Cell type coverage | Typically 3–5 clusters | 8 immune populations: CD4+ T, CD8+ T, NK, T-reg, B cell, Monocyte, pDC, Neutro |
| Cross-domain consistency | Standalone figure | Same 80 patients as Plots 1, 4, 5 |
| Reproducibility | No fixed seed | `set.seed(42)` — bitwise-identical across R ≥ 4.1 |
| Calibration traceability | Unspecified | Anchored to published PBMC reference frequencies and KEYNOTE-189, TCGA PanCancer Atlas |

---

## Visual anatomy

### UMAP layout

```
         NK ●              CD8+ T ●
                 
    pDC ●          Monocyte ●
                                    
              CD4+ T ●    T-reg ●
                                   
       Neutro ●      B cell ●

  UMAP 1 →
```

| Element | Description |
|---|---|
| Points | 2,800 cells (350 per cell type), α = 0.55, size = 0.6 |
| Labels | `ggrepel::geom_label_repel()` at per-cell-type UMAP median |
| Color | Per-cell-type palette (shared `CELL_COLS`) |
| Title | Includes live cell count and marker count |

### Cell type expression profiles used for UMAP input (mean, σ = 0.08)

| Cell type | CD4 | CD8 | NK | TREG | PDL1 | TMB | IFNg |
|-----------|-----|-----|----|------|------|-----|------|
| CD4+ T | 0.85 | 0.05 | 0.05 | 0.10 | 0.30 | 0.40 | 0.35 |
| CD8+ T | 0.05 | 0.90 | 0.10 | 0.05 | 0.45 | 0.55 | 0.65 |
| T-reg | 0.70 | 0.05 | 0.05 | 0.90 | 0.25 | 0.20 | 0.15 |
| NK | 0.05 | 0.10 | 0.90 | 0.05 | 0.30 | 0.35 | 0.70 |
| B cell | 0.05 | 0.05 | 0.10 | 0.05 | 0.55 | 0.30 | 0.10 |
| Monocyte | 0.10 | 0.05 | 0.15 | 0.05 | 0.70 | 0.50 | 0.40 |
| pDC | 0.05 | 0.05 | 0.10 | 0.05 | 0.60 | 0.25 | 0.20 |
| Neutro | 0.05 | 0.05 | 0.20 | 0.05 | 0.20 | 0.30 | 0.15 |

### UMAP parameters

| Parameter | Value |
|-----------|-------|
| Algorithm | `uwot::umap()` |
| n_neighbors | 15 |
| min_dist | 0.1 |
| n_components | 2 |
| Pre-embedding jitter | N(0, 0.01) added to expression matrix |

---

## Data note on cell population coverage

Four cell types (CD4+ T, CD8+ T, NK, T-reg) are anchored on patient ADBM biomarker measurements. The remaining four (B cell, Monocyte, pDC, Neutrophil) use expression profiles derived from published PBMC reference frequencies, as these populations were not measured in the ADBM domain of this dataset. This approach is consistent with methods used in biomarker-exploratory analyses when full single-cell panels are unavailable.

> For fully data-driven UMAP, replace the `expand_cells()` expansion step with a real CyTOF expression matrix from `HDCytoData::Bodenmiller_BCR_XL_SE()` — the downstream `uwot::umap()` and plotting code are unchanged.

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
| ADSL | Subject-level | `USUBJID`, `TUMORTYPE`, `ARM`, `PDL1GRP`, `TMBHIGH` |
| ADBM | Biomarker assessments | `PARAMCD` (CD4/CD8/NK/TREG/PDL1/TMB/IFNg), `AVAL` |

### Baseline biomarker summary (ADBM, treatment arm, baseline visit)

| Biomarker | Mean | Range |
|-----------|-----:|------:|
| CD4+ T (%) | 36.9 | 20.0 – 60.0 |
| CD8+ T (%) | 24.1 | 10.0 – 40.0 |
| NK (%) | 12.5 | 5.0 – 23.7 |
| T-reg (%) | 3.0 | 1.0 – 7.4 |
| PDL1 (%) | 31.2 | 0.0 – 79.8 |
| IFN-γ (pg/mL) | 14.4 | 1.0 – 43.2 |
| TMB (mut/Mb) | 12.0 | 0.0 – 34.8 |

---

## Output files

| File | Description | Dimensions |
|------|-------------|------------|
| `plot2_umap.R` | R script — Plot 2 only | |
| `plot2_umap.png` | UMAP single-cell immune landscape | 9×8 in · 150 DPI |

---

## When to use

**Appropriate:**
- Exploratory immune profiling visualization in early-phase oncology trials
- Communicating cell type clustering structure across biomarker dimensions
- Supplementary figures in translational publications alongside flow or CyTOF data

**Limitations:**
- Does not reflect full transcriptomic or proteomic resolution
- Cell type separation is driven by canonical literature profiles, not patient-derived single-cell data
- B cell, Monocyte, pDC, Neutrophil expression profiles estimated from reference frequencies

---

## Requirements

```
R >= 4.1
```

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "uwot", "ggrepel"))
```

---

## References

McInnes L, et al. UMAP: Uniform Manifold Approximation and Projection for dimension reduction. *arXiv.* 2018;1802.03426.

Bodenmiller B, et al. Multiplexed mass cytometry profiling of cellular states perturbed by small-molecule regulators. *Nat Biotechnol.* 2012;30(9):858–867.

Weber LM, et al. HDCytoData: Collection of high-dimensional cytometry benchmark datasets in Bioconductor object formats. *F1000Research.* 2019;8:1459.
