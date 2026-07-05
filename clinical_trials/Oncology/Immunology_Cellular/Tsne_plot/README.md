# t-SNE Single-Cell Immune Landscape

A t-SNE dimensionality reduction of the same 2,800-cell × 7-marker matrix used in Plot 2, providing an orthogonal single-cell view for cross-validation of cluster structure against UMAP, anchored on patient ADBM baseline distributions from the ONCVIZ-001 basket trial dataset.

**Dataset:** ONCVIZ-001 · N = 80 | **Cutoff:** 05 March 2026 | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

UMAP and t-SNE preserve different aspects of high-dimensional structure. Most immune profiling publications present one or the other but not both on the same dataset. By running t-SNE on the identical cell × marker matrix as Plot 2, this plot enables direct visual comparison of local and global neighborhood preservation.

| | Prior art | This work |
|---|---|---|
| Orthogonal DR view | UMAP or t-SNE separately, different datasets | t-SNE on exact same matrix as UMAP (Plot 2) |
| Cell type coverage | Typically 3–5 clusters | 8 immune populations: CD4+ T, CD8+ T, NK, T-reg, B cell, Monocyte, pDC, Neutro |
| Cross-domain consistency | Standalone figure | Same 80 patients, same expression matrix as Plots 2, 4, 5 |
| Reproducibility | No fixed seed | `set.seed(42)` — bitwise-identical across R ≥ 4.1 |
| Calibration traceability | Unspecified | Anchored to published PBMC reference frequencies and KEYNOTE-189, TCGA PanCancer Atlas |

---

## Visual anatomy

### t-SNE layout

```
    CD4+ T ●        T-reg ●
                  
       B cell ●            NK ●
                  
    Monocyte ●      CD8+ T ●
                  
        pDC ●       Neutro ●

  t-SNE 1 →
```

| Element | Description |
|---|---|
| Points | 2,800 cells (350 per cell type), α = 0.55, size = 0.6 |
| Labels | `ggrepel::geom_label_repel()` at per-cell-type t-SNE median |
| Color | Per-cell-type palette (shared `CELL_COLS`) |
| Title | Includes live cell count and perplexity value |

### t-SNE parameters

| Parameter | Value |
|-----------|-------|
| Algorithm | `Rtsne::Rtsne()` |
| dims | 2 |
| perplexity | min(40, ⌊(n − 1) / 3⌋) |
| max_iter | 1,000 |
| check_duplicates | FALSE |
| Input matrix | Same `mat_sc` as Plot 2 (2,800 × 7) |

### Comparison with UMAP (Plot 2)

| Property | UMAP | t-SNE |
|----------|------|-------|
| Global structure | Better preserved | Local neighborhoods prioritized |
| Inter-cluster distances | Meaningful | Not interpretable |
| Speed at n = 2,800 | Faster | Slower |
| Key parameter | n_neighbors = 15 | perplexity = 40 (capped) |

---

## Data note on cell population coverage

Four cell types (CD4+ T, CD8+ T, NK, T-reg) are anchored on patient ADBM biomarker measurements. The remaining four (B cell, Monocyte, pDC, Neutrophil) use expression profiles derived from published PBMC reference frequencies, as these populations were not measured in the ADBM domain of this dataset. The expression matrix is identical to Plot 2 — only the embedding algorithm differs.

> For fully data-driven t-SNE, replace the `expand_cells()` expansion step with a real CyTOF expression matrix from `HDCytoData::Bodenmiller_BCR_XL_SE()` — the downstream `Rtsne()` and plotting code are unchanged.

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
| `plot3_tsne.R` | R script — Plot 3 only | |
| `plot3_tsne.png` | t-SNE single-cell immune landscape | 9×8 in · 150 DPI |

---

## When to use

**Appropriate:**
- Cross-validating UMAP cluster structure with an orthogonal algorithm
- Exploratory immune profiling in early-phase oncology translational reports
- Supplementary figures paired with Plot 2 to demonstrate local neighborhood robustness

**Limitations:**
- t-SNE inter-cluster distances are not biologically interpretable
- Does not reflect full transcriptomic or proteomic resolution
- B cell, Monocyte, pDC, Neutrophil expression profiles estimated from reference frequencies

---

## Requirements

```
R >= 4.1
```

```r
install.packages(c("ggplot2", "dplyr", "tidyr", "Rtsne", "uwot", "ggrepel"))
```

---

## References

van der Maaten L, Hinton G. Visualizing data using t-SNE. *J Mach Learn Res.* 2008;9:2579–2605.

Bodenmiller B, et al. Multiplexed mass cytometry profiling of cellular states perturbed by small-molecule regulators. *Nat Biotechnol.* 2012;30(9):858–867.

Weber LM, et al. HDCytoData: Collection of high-dimensional cytometry benchmark datasets in Bioconductor object formats. *F1000Research.* 2019;8:1459.
