# Manhattan Plot — Gene-Level Association with TMB-High Status

A publication-ready Manhattan-style plot displaying gene-level association p-values (Fisher's exact test) between somatic mutation status and TMB-High (tumor mutational burden ≥ threshold) classification, across 15 genes positioned by chromosomal cytoband location, with alternating chromosome shading and a nominal significance threshold line.

**Dataset:** ONCVIZ-001 · N = 80 patients · ADMUT / ADSL · 15 genes · 138 variants | **Language:** R · ggplot2 · ggrepel | **License:** CC BY 4.0

---

## The gap this fills

Manhattan plots are the canonical visualization for genome-wide association studies (GWAS) but are rarely used in clinical genomics for gene-panel association analyses. Applying the Manhattan layout to targeted panel data — with genes positioned at their cytogenetic loci rather than at SNP base-pair coordinates — provides a spatially-organized alternative to ranked bar charts or volcano plots for biomarker discovery in trial datasets.

| | GWAS Manhattan | This work (panel-level) |
|---|---|---|
| x-axis unit | Base-pair position | Cytoband position (gene-level, jittered) |
| Markers | SNPs (millions) | Somatic variants (138) |
| Multiple testing threshold | Bonferroni 5×10⁻⁸ | Nominal p = 0.05 (exploratory) |
| Test | Additive genotype model | Fisher's exact test (2×2 contingency: mutant vs. TMB-H/L) |
| Labels | Top SNPs | Gene-level peak variant |
| Alternating shading | By chromosome | By chromosome (10 chromosomes shown) |
| Reproducibility | No fixed seed | `set.seed(8)` · bitwise-identical |

---

## Visual anatomy

```
  −log₁₀(p)
  3.0  ┤             ●TP53          ●KRAS
       │
  2.0  ┤   ●ARID1A        ●PTEN
       │
  1.0  ┤   ○  ○    ○   ○   ○  ○   ○   ○  ○   ○
  ─────┼─────────────────────────────────── p = 0.05 (dashed, amber)
  0.0  ┤   ○  ○    ○
       ├──────────────────────────────────────────────
       chr1 chr3  chr7 chr9 chr10 chr12 chr13 chr17 chr18 chr19
             (alternating gray/white chromosome bands)
```

| Element | Description |
|---|---|
| Alternating gray bands | Even-numbered chromosome columns · `geom_rect + alpha=0.7` |
| Dashed amber line | Nominal p = 0.05 threshold · `−log₁₀(0.05) ≈ 1.30` |
| Dark blue points | Notable associations (p < 0.15) |
| Light blue points | Non-notable associations |
| Larger points | Notable associations · `size = 2.5 vs 1.8` |
| `ggrepel` labels | Gene name at peak variant position (p < 0.5) |
| Jitter | ±0.28 horizontal jitter per chromosome (within chromosome band) |

---

## Output figure

| File | Dimensions | DPI |
|---|---|---|
| `plots/08_manhattan.png` | 15 × 6.5 in | 180 |

---

## Synthetic dataset — ONCVIZ-001 ADaM v1 (ADMUT / ADSL)

```
n patients     80
n variants     138 somatic point mutations / indels
n genes        15
Chromosomes    10 (chr1 · chr3 · chr7 · chr9 · chr10 · chr12 · chr13 · chr17 · chr18 · chr19)
Comparison     TMB-High (n=32) vs TMB-Low (n=48)
Seed           8 · fully reproducible
```

### Gene-chromosome mapping (cytoband positions)

| Gene | Chromosome | Cytoband |
|---|---|---|
| ARID1A | chr1 | 1p36.11 |
| PIK3CA | chr3 | 3q26.3 |
| EGFR | chr7 | 7p11.2 |
| MET | chr7 | 7q31.2 |
| BRAF | chr7 | 7q34 |
| CDKN2A | chr9 | 9p21.3 |
| RET | chr10 | 10q11.21 |
| PTEN | chr10 | 10q23.3 |
| KRAS | chr12 | 12p12.1 |
| RB1 | chr13 | 13q14.2 |
| TP53 | chr17 | 17p13.1 |
| NF1 | chr17 | 17q11.2 |
| SMAD4 | chr18 | 18q21.2 |
| KEAP1 | chr19 | 19p13.2 |
| STK11 | chr19 | 19p13.3 |

### Variables used

| Variable | Source | Description |
|---|---|---|
| `USUBJID` | ADMUT | Subject identifier |
| `HUGO_SYMBOL` | ADMUT | Gene name |
| `PROTEIN_POS` | ADMUT | Amino acid position (used for horizontal jitter) |
| `TMBHIGH` | ADSL | TMB-High flag ("Y"/"N") |

---

## Statistical method

For each gene, a **2×2 Fisher's exact test** is computed:

```
              TMB-High  TMB-Low
  Mutant        a          c
  Wildtype      b          d

Fisher's exact p-value = fisher.test(matrix(c(a,b,c,d), 2))$p.value
```

`−log₁₀(p)` is plotted on the y-axis. The nominal threshold (p = 0.05) is rendered as a dashed amber line at y = 1.30. Each individual somatic variant inherits the p-value of its gene (gene-level association), so multiple points per gene appear at the same y-level with horizontal jitter.

**Appropriate when:** Exploring which genes in a targeted panel are enriched for mutations in TMB-High patients (hypothesis-generating).
**Not appropriate when:** Drawing causal conclusions; controlling for tumor type confounding; analyzing whole-genome associations.

---

## Key parameters

| Parameter | Value | Description |
|---|---|---|
| Statistical test | Fisher's exact (2×2) | Per gene: mutant vs. TMB-H/L |
| Significance threshold | p = 0.05 (nominal) | No Bonferroni correction |
| Notable threshold | p < 0.15 | For point size/color distinction |
| Label threshold | p < 0.5 | `ggrepel` gene labels |
| Horizontal jitter | ±0.28 | `runif(n, −0.28, 0.28)` within chromosome |
| Even-chr band alpha | 0.7 | Alternating gray background |
| Color: notable | `#2c5f8a` | Blue |
| Color: non-notable | `#7fb3d3` | Light blue |
| Figure dimensions | 15 × 6.5 in | Wide landscape (GWAS aspect ratio) |
| DPI | 180 | Publication quality |

---

## Limitations

- **Gene-level resolution only:** Unlike GWAS, this does not resolve which specific variant within a gene drives the association. Individual variant-level tests would require larger cohort sizes.
- **Confounding by tumor type:** Some genes are preferentially mutated in specific histologies (e.g., KRAS in CRC/PDAC), and TMB-High prevalence also varies by histology. A stratified or adjusted test is needed to separate gene–TMB association from the histology confounder.
- **Protein position as genomic coordinate proxy:** The x-axis jitter uses `PROTEIN_POS` (amino acid position), not true genomic base-pair coordinates. Multiple genes on the same chromosome are therefore not ordered by genomic position within the chromosome band.
- **No multiple testing correction:** With 15 genes, 1 false positive is expected at α = 0.05. Apply Bonferroni (threshold: p < 0.0033) or Benjamini–Hochberg FDR for formal reporting.

---

## Requirements

```r
library(ggplot2)   # >= 3.4
library(ggrepel)   # >= 0.9
library(dplyr)     # >= 1.1
```

---

## Files

| File | Description |
|---|---|
| `08_manhattan.R` | Self-contained script · synthetic data generated internally |
| `plots/08_manhattan.png` | Output figure · 15 × 6.5 in · 180 DPI |

---

## References

Turner S. qqman: an R package for visualizing GWAS results using Q-Q and manhattan plots. *J Open Source Softw.* 2018;3(25):731.

Purcell S, et al. PLINK: A Tool Set for Whole-Genome Association and Population-Based Linkage Analyses. *Am J Hum Genet.* 2007;81(3):559–575.

Alexandrov LB, et al. Mutational signatures associated with tobacco smoking in human cancer. *Science.* 2016;354(6312):618–622.
