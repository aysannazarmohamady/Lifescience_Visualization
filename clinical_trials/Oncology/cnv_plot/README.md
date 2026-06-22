# CNV Plot — Gene-Level Copy Number Alterations

A three-panel publication-ready figure displaying copy number variation (CNV) scores across 10 key oncogenes and tumor suppressors in 49 patients: a tumor-type track (top), a diverging-color heatmap (center) ordered by tumor type, and a signed summary bar chart (right) showing the number of patients with amplifications vs. deletions per gene.

**Dataset:** ONCVIZ-001 · N = 49 patients with ≥1 CNV event · ADMUT | **Language:** R · ggplot2 · patchwork | **License:** CC BY 4.0

---

## The gap this fills

Copy number visualization in oncology is dominated by genome-wide circular or linear segment plots (e.g., CNVkit, GISTIC2 circus plots). These require whole-genome or whole-exome data and specialized tools. For targeted panel data — where copy number is reported per gene rather than per segment — a gene × patient heatmap with diverging color scale and marginal summary bars communicates the landscape more efficiently.

| | Prior art | This work |
|---|---|---|
| Data type assumed | Genome-wide segments | Gene-level score (targeted panel) |
| Color scale | Binary (amp/del) | Seven-level diverging scale (Deep Del → Amp) |
| Patient ordering | Random or alphabetical | By tumor type (clinical grouping) |
| Amplitude summary | Not shown | Signed bar: amplification (red) vs. deletion (blue) per gene |
| Cytogenetic annotation | Not shown | Gene label includes cytoband (e.g., "TP53\n(17p13.1)") |
| Tumor type track | Not shown | Tile track above heatmap |
| Reproducibility | No fixed seed | `set.seed(6)` · bitwise-identical |

---

## Visual anatomy

```
  Tumor type  ■NSCLC■■NSCLC■■BRCA■■HCC■■CRC■■CRC■■PDAC■  (tile track)
  ─────────────────────────────────────────────────────────────────────
  PIK3CA      ░░░░  ████  ░░░   ██  ░░░░      ████    │  amp: ██ (n=8)
  (3q26.3)    ░░░░  ████  ░░░             ████        │  del: █  (n=4)
  EGFR        ████  ░░░░        ████  ░░░░            │  amp: ████
  (7p11.2)                                            │  del: ██
  ...                                                 │
  SMAD4       ░░░░  ████  ░░░░  ████                  │
  (18q21.2)                                           │
  ─────────────────────────────────────────────────────────────────────
  Patient columns (n=49, ordered by tumor type) →      ← Del  Amp →
```

| Element | Description |
|---|---|
| Tumor type track | `geom_tile` · ONCVIZ-001 5-histology palette |
| Heatmap | `geom_tile` · diverging color: deep blue (−4) → white (0) → dark red (+4) |
| Gray cells | `NA` score = no CNV detected at this gene in this patient |
| Right bar chart | Signed `geom_col`: positive (red) = amplification, negative (blue) = deletion |

---

## Output figure

| File | Dimensions | DPI |
|---|---|---|
| `plots/06_cnv.png` | 17 × 9 in | 180 |

---

## Synthetic dataset — ONCVIZ-001 ADaM v1 (ADMUT)

```
n patients (CNV)  49
Genes             10: PIK3CA · EGFR · MET · BRAF · CDKN2A ·
                      PTEN · KRAS · RB1 · NF1 · SMAD4
Variant type      CNV only (VARIANT_TYPE == "CNV")
Seed              6 · fully reproducible
```

### CNV score encoding

| Score | Variant class | Clinical interpretation |
|---|---|---|
| +4 | Amplification | High-level amplification (typically ≥6 copies) |
| +3 | High_Gain | Moderate gain |
| +1 | Gain | Low-level gain |
| 0 | Copy_Neutral | Diploid |
| −1 | Shallow_Deletion | Heterozygous loss |
| −2 | Deletion | Homozygous or deep loss |
| −4 | Deep_Deletion | Homozygous deletion (CDKN2A, PTEN, RB1 typical) |

### Cytogenetic locations

| Gene | Cytoband | Alteration type (typical) |
|---|---|---|
| PIK3CA | 3q26.3 | Amplification (gain-of-function) |
| EGFR | 7p11.2 | Amplification (gain-of-function) |
| MET | 7q31.2 | Amplification / high gain |
| BRAF | 7q34 | Amplification |
| CDKN2A | 9p21.3 | Deep deletion (loss-of-function) |
| PTEN | 10q23.3 | Deletion (loss-of-function) |
| KRAS | 12p12.1 | Amplification (gain-of-function) |
| RB1 | 13q14.2 | Deep deletion (loss-of-function) |
| NF1 | 17q11.2 | Deletion |
| SMAD4 | 18q21.2 | Deletion (loss-of-function) |

### Variables used

| Variable | Source | Description |
|---|---|---|
| `USUBJID` | ADMUT | Subject identifier |
| `HUGO_SYMBOL` | ADMUT | Gene name |
| `VARIANT_TYPE` | ADMUT | Filtered to "CNV" |
| `VARIANT_CLASS` | ADMUT | CNV category (Amplification · Deletion · etc.) |
| `TUMORTYPE` | ADSL (joined) | Histology for patient ordering and track color |

---

## Key parameters

| Parameter | Value | Description |
|---|---|---|
| CNV score range | [−4, 4] | Fixed for consistent diverging scale |
| Color: deep blue | `#1a5276` | Score = −4 (Deep_Deletion) |
| Color: white | `#f8f9fa` | Score = 0 (Copy_Neutral) |
| Color: dark red | `#922b21` | Score = +4 (Amplification) |
| Patient ordering | By tumor type | Groups histologies visually |
| Layout heights | `c(1, 12)` | Track : heatmap ratio |
| Layout widths | `c(13, 2)` | Heatmap : summary bar ratio |
| Bar fill: amp | `#c0392b` | Red (positive scores) |
| Bar fill: del | `#1a5276` | Blue (negative scores) |
| Figure dimensions | 17 × 9 in | Wide landscape |
| DPI | 180 | Publication quality |

---

## Statistical considerations

- **Score as ordinal proxy:** The integer score is a discretized representation of copy number ratio; it is not a continuous log₂ ratio as produced by CNVkit or GATK CNV. Treat it as categorical for visualization purposes.
- **Absence of data vs. Copy_Neutral:** Cells with `NA` score (gray) represent genes not detected in the panel for that patient — they should not be interpreted as copy-neutral without reviewing the sequencing report.
- **Frequency summary:** The right bar panel counts patients with *any* positive score (amp) or *any* negative score (del). A patient with a score of +1 (Gain) and a score of −2 (Deletion) at the same gene in different clones cannot be represented in this summary; in practice, the dominant clone score should be used.

---

## Limitations

- This plot requires gene-level copy number calls (integer or discrete score). It is not suitable for genome-wide raw copy number ratios (use a circos or genome-wide segment plot instead).
- With n = 49 patients, patient column widths may become narrow; increasing `width` in `ggsave()` or subsetting to fewer patients is recommended for larger cohorts.
- The summary bar does not distinguish between patients with different score magnitudes (a patient with +1 Gain counts the same as +4 Amplification). A score-weighted bar or frequency breakdown by score level is a refinement.

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
| `06_cnv.R` | Self-contained script · synthetic data generated internally |
| `plots/06_cnv.png` | Output figure · 17 × 9 in · 180 DPI |

---

## References

Talevich E, et al. CNVkit: Genome-Wide Copy Number Detection and Visualization from Targeted DNA Sequencing. *PLOS Comput Biol.* 2016;12(4):e1004873.

Mermel CH, et al. GISTIC2.0 facilitates sensitive and confident localization of the targets of focal somatic copy-number alteration in human cancers. *Genome Biol.* 2011;12(4):R41.

Cancer Genome Atlas Network. Comprehensive molecular portraits of human breast tumours. *Nature.* 2012;490(7418):61–70.
