# Circos Plot — Genomic Alteration Overview, Top 8 Mutated Genes

A polar bar chart approximating a circos-style genomic overview, displaying the alteration frequency (% of 80 patients) for the 8 most commonly mutated genes in the ONCVIZ-001 synthetic basket trial. Arc height encodes alteration prevalence; each gene is assigned a distinct color.

**Dataset:** ONCVIZ-001 · N = 80 patients · ADMUT | **Language:** R · ggplot2 (`coord_polar`) | **License:** CC BY 4.0

---

## The gap this fills

True circos plots — with inter-gene chord arcs representing co-occurrence or translocation partners — require the `circlize` Bioconductor package, which introduces a non-ggplot2 grammar and significant setup overhead. For a high-level summary of alteration frequencies in a visually distinctive circular format (common in publications for visual appeal), a `coord_polar` bar chart in ggplot2 achieves the same aesthetic with zero additional dependencies.

| | Full circlize circos | This work (ggplot2 polar bar) |
|---|---|---|
| Package dependency | `circlize` (Bioconductor) | Base `ggplot2` only |
| Co-occurrence chords | ✓ | ✗ (omitted by design; noted in caption) |
| Alteration frequency arcs | ✓ | ✓ via `geom_col + coord_polar` |
| Color per gene | ✓ | ✓ |
| Labels with frequency | Optional | Inline · `geom_text` above each arc |
| Reproducibility | Varies | `set.seed(7)` · bitwise-identical |

---

## Visual anatomy

```
              TP53 (48%)
          ████████████████
         ██              ██
  KEAP1 ██                ██ KRAS (30%)
  (10%) ██                ██
         ██              ██
         ██   (center)   ██
          ██            ██  SMAD4 (9%)
  PTEN   ██            ██
  (12%)   ████████████████
          CDKN2A  PIK3CA
          (11%)   (11%)  ARID1A (12%)
```

| Element | Description |
|---|---|
| Arc (bar) | Alteration frequency · `geom_col` · height ∝ % altered |
| Color | Per gene · 8 distinct colors |
| Label | Gene name + % · positioned above arc tip with `geom_text` |
| Coordinate system | `coord_polar()` · bars converted to radial arcs |
| Note | Co-occurrence chord arcs are not implemented (requires `circlize`) |

---

## Output figure

| File | Dimensions | DPI |
|---|---|---|
| `plots/07_circos.png` | 13 × 13 in | 180 |

---

## Synthetic dataset — ONCVIZ-001 ADaM v1 (ADMUT)

```
n patients     80
Genes shown    8: ARID1A · PIK3CA · CDKN2A · PTEN · KRAS · SMAD4 · TP53 · KEAP1
Frequencies    Hard-coded from pan-cohort prevalence (pre-computed from ADMUT)
Seed           7
```

### Alteration frequencies (pan-cohort %)

| Gene | % Altered | Color |
|---|---|---|
| TP53 | 48% | `#6c3483` (Dark purple) |
| KRAS | 30% | `#c0392b` (Red) |
| PTEN | 12% | `#8e44ad` (Purple) |
| ARID1A | 12% | `#2980b9` (Blue) |
| CDKN2A | 11% | `#27ae60` (Green) |
| PIK3CA | 11% | `#e67e22` (Orange) |
| KEAP1 | 10% | `#1abc9c` (Teal) |
| SMAD4 | 9% | `#16a085` (Dark teal) |

---

## Key parameters

| Parameter | Value | Description |
|---|---|---|
| Y-axis limits | `c(0, 65)` | Provides headroom above highest arc (48%) |
| Bar width | 0.7 | Within polar coordinate |
| Bar alpha | 0.85 | Slight transparency |
| Label y-position | `pct + 4` | 4 percentage-points above arc tip |
| Coordinate | `coord_polar()` | Converts bar chart to radial layout |
| Figure dimensions | 13 × 13 in | Square — optimal for circular layout |
| DPI | 180 | Publication quality |

---

## Statistical considerations

- **Frequency, not statistical test:** This is a descriptive frequency plot. No significance test is applied; the arc height represents simple prevalence.
- **No chord arcs:** Co-occurrence relationships (e.g., TP53 + KRAS co-mutation) require computing a pairwise co-occurrence matrix and rendering chord segments via `circlize::chordDiagram()`. This is the most informative additional layer for a true circos figure.
- **Gene selection:** The 8 genes shown are pre-selected by pan-cohort frequency rank. For a different gene panel or ranking criterion (e.g., clinical actionability), the frequency vector should be recomputed from the source ADMUT.

---

## Limitations

- `coord_polar()` in ggplot2 does not support chord arcs, segment curvature, or multi-track rings. For publication-grade circular figures with co-occurrence chords, use `circlize` or `ggforce`.
- The arc for SMAD4 (9%) is visually short; in polar coordinates, small differences in low-frequency genes are harder to discern than in a standard bar chart. A minimum arc height can be added artificially but would distort the frequency encoding.
- Labels may overlap when multiple genes have similar frequencies (PTEN 12% vs. ARID1A 12%); manual `nudge_x` or `nudge_y` adjustments may be needed for real data.

---

## Requirements

```r
library(ggplot2)   # >= 3.4
```

---

## Files

| File | Description |
|---|---|
| `07_circos.R` | Self-contained script · synthetic data generated internally |
| `plots/07_circos.png` | Output figure · 13 × 13 in · 180 DPI |

---

## References

Krzywinski MI, et al. Circos: An information aesthetic for comparative genomics. *Genome Res.* 2009;19(9):1639–1645.

Gu Z, et al. circlize implements and enhances circular visualization in R. *Bioinformatics.* 2014;30(19):2811–2812.

Cancer Genome Atlas Research Network. Comprehensive molecular profiling of lung adenocarcinoma. *Nature.* 2014;511(7511):543–550.
