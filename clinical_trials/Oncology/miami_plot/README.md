# Miami Plot — Gene-Level Association with TMB-High (top) vs. MSI-High (bottom)

A publication-ready Miami plot (mirrored double Manhattan plot) simultaneously displaying gene-level Fisher's exact test associations with two binary biomarker endpoints — TMB-High (top panel, y ≥ 0) and MSI-High (bottom panel, y ≤ 0) — across 15 genes positioned by chromosomal cytoband. The mirrored layout enables direct visual comparison of which genes associate with each endpoint.

**Dataset:** ONCVIZ-001 · N = 80 patients · ADMUT / ADSL · 15 genes · 138 variants | **Language:** R · ggplot2 · ggrepel | **License:** CC BY 4.0

---

## The gap this fills

Miami plots (also called mirrored or back-to-back Manhattan plots) originate in epigenomics, where they are used to contrast differentially methylated positions between two conditions. Their application to somatic mutation association studies is uncommon but highly informative in oncology basket trials, where multiple biomarker endpoints (TMB, MSI, PDL1) need to be compared side-by-side for the same gene panel.

| | Standard EWAS/GWAS Miami | This work |
|---|---|---|
| Data type | Methylation probes / SNPs | Somatic mutations (targeted panel) |
| Top panel comparison | One condition | TMB-High enrichment |
| Bottom panel comparison | Other condition | MSI-High enrichment |
| x-axis unit | Genomic coordinate | Cytoband position (gene-level, jittered) |
| Statistical test | Linear/logistic model | Fisher's exact test (2×2 per gene) |
| Threshold lines | Bonferroni | Nominal p = 0.05 (amber dashed, per panel) |
| Color coding | By chromosome | By endpoint (blue = TMB-H · red = MSI-H) |
| Reproducibility | No fixed seed | `set.seed(9)` · bitwise-identical |

---

## Visual anatomy

```
  −log₁₀(p)  [TOP PANEL — TMB-High association]
  3.0  ┤           ●TP53         ●KRAS
  2.0  ┤  ●ARID1A       ●PTEN
  1.0  ┤  ○   ○   ○  ○  ○  ○   ○  ○   ○
  ─────┼────────────────────────── p = 0.05 (amber dashed)
  0.0  ════════════════════════════════════ (zero axis)
  ─────┼────────────────────────── p = 0.05 (amber dashed)
  −1.0 ┤  ○   ○   ○  ○  ○  ○   ○  ○   ○
  −2.0 ┤     ●PTEN         ●ARID1A
  −3.0 ┤                            (none at this level)
       [BOTTOM PANEL — MSI-High association]
       chr1 chr3  chr7 chr9 chr10 chr12 chr13 chr17 chr18 chr19
```

| Element | Description |
|---|---|
| Top panel (y ≥ 0) | TMB-High association · blue points |
| Bottom panel (y ≤ 0) | MSI-High association · red points · y = −log₁₀(p) (negated) |
| Amber dashed lines | Nominal p = 0.05 threshold in each panel |
| `ggrepel` labels | Gene name at peak absolute y per gene per panel |
| `facet_wrap` | Two panels · `scales = "free_y"` · `strip.placement = "outside"` |
| Y-axis labels | `labels = abs` — shows absolute values (p-value magnitude) |

---

## Output figure

| File | Dimensions | DPI |
|---|---|---|
| `plots/09_miami.png` | 16 × 10 in | 180 |

---

## Synthetic dataset — ONCVIZ-001 ADaM v1 (ADMUT / ADSL)

```
n patients       80
n variants       138 somatic point mutations / indels
n genes          15
Chromosomes      10 (chr1 · chr3 · chr7 · chr9 · chr10 · chr12 · chr13 · chr17 · chr18 · chr19)
Top comparison   TMB-High (Y/N) — n(Y) ≈ 32 · n(N) ≈ 48
Bottom comparsn  MSI-High (MSI-H / MSS) — n(MSI-H) ≈ 4 · n(MSS) ≈ 76
Seed             9 · fully reproducible
```

### Variables used

| Variable | Source | Description |
|---|---|---|
| `USUBJID` | ADMUT | Subject identifier |
| `HUGO_SYMBOL` | ADMUT | Gene name |
| `TMBHIGH` | ADSL | TMB-High flag ("Y"/"N") |
| `MSISTS` | ADSL | MSI status ("MSI-H" · "MSS") |

---

## Statistical method

For each gene and each endpoint, a **2×2 Fisher's exact test** is computed:

```
Top panel (TMB-High):           Bottom panel (MSI-High):
              TMB-H  TMB-L                   MSI-H  MSS
  Mutant       a      c        Mutant          a     c
  Wildtype     b      d        Wildtype        b     d
```

Top panel y-coordinate: `+log₁₀(p_tmb)`
Bottom panel y-coordinate: `−log₁₀(p_msi)` (negated for mirror)

`scale_y_continuous(labels = abs)` ensures both panels display positive tick labels despite negative y-values in the bottom panel.

---

## Key parameters

| Parameter | Value | Description |
|---|---|---|
| Statistical test | Fisher's exact (2×2) | Per gene per endpoint |
| Top panel color | `#2c5f8a` | Blue (TMB-High) |
| Bottom panel color | `#c0392b` | Red (MSI-High) |
| Significance threshold | p = 0.05 | Amber dashed line in each panel |
| Label threshold | `abs(log10p_dir) > 0.5` | ~p < 0.32 |
| Horizontal jitter | ±0.28 | Independent seed per panel |
| Facet layout | `ncol=1 · scales="free_y"` | Stacked panels |
| Strip position | "outside" | Panel labels on left side |
| Figure dimensions | 16 × 10 in | Tall to accommodate both panels |
| DPI | 180 | Publication quality |

---

## Statistical considerations

- **MSI-High arm is very small (n ≈ 4):** Fisher's exact test is the correct choice (not chi-squared) for small expected cell counts, but with n = 4 MSI-H patients, statistical power is extremely low. The bottom panel of this Miami plot is primarily illustrative for a cohort this size.
- **Independence of tests:** TMB and MSI are biologically correlated (MSI-H tumors typically have high TMB). Associations in both panels for the same gene (e.g., PTEN) may reflect the same underlying biology rather than independent signals.
- **No multiple testing correction:** With 15 × 2 = 30 simultaneous tests, ~1.5 false positives are expected at α = 0.05.

---

## Limitations

- True Miami plots from EWAS/GWAS have millions of probes/SNPs. With 138 variants across 15 genes, this figure is sparsely populated — the "Manhattan" visual is more schematic than epidemiological.
- The bottom panel's statistical power is severely limited by the small MSI-H population (n ≈ 4). Results should be treated as hypothesis-generating only.
- Genes with the same chromosome (e.g., EGFR, MET, BRAF on chr7) are plotted at the same chromosome x-position with jitter; their relative ordering within the chromosome is not genomically accurate.

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
| `09_miami.R` | Self-contained script · synthetic data generated internally |
| `plots/09_miami.png` | Output figure · 16 × 10 in · 180 DPI |

---

## References

Jaffe AE, et al. Bump hunting to identify differentially methylated regions in epigenetic epidemiology studies. *Int J Epidemiol.* 2012;41(1):200–209.

Le DT, et al. PD-1 Blockade in Tumors with Mismatch-Repair Deficiency. *N Engl J Med.* 2015;372(26):2509–2520.

Goodman AM, et al. Tumor Mutational Burden as an Independent Predictor of Response to Immunotherapy in Diverse Cancers. *Mol Cancer Ther.* 2017;16(11):2598–2608.
