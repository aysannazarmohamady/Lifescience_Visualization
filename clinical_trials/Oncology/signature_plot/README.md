# Mutational Signature Plot — COSMIC SBS Signature Contributions

A dual-panel publication-ready figure showing (left) per-patient mutational signature decomposition as stacked 100% bar charts ordered by tumor type and dominant signature, and (right) mean signature weight per tumor type as a stacked horizontal bar. Twelve COSMIC v3.3 SBS signatures are displayed, covering known aetiologies from tobacco smoking and APOBEC activity to homologous recombination deficiency and mismatch repair.

**Dataset:** ONCVIZ-001 · N = 80 patients · ADSIG | **Signatures:** 12 COSMIC v3.3 SBS | **Language:** R · ggplot2 · patchwork | **License:** CC BY 4.0

---

## The gap this fills

Mutational signature decomposition is routinely performed with tools such as SigProfiler, deconstructSigs, or MutationalPatterns. The resulting weight matrices are typically presented as either (1) a per-signature bar chart collapsed across all patients, or (2) a pie chart per patient — both of which lose the cross-patient heterogeneity that is clinically meaningful in basket trials enrolling multiple histologies.

| | Prior art | This work |
|---|---|---|
| Per-patient resolution | Collapsed or pie chart | Full stacked 100% bar per patient |
| Cross-histology comparison | Separate figures or facets | Synchronized left + right panels via `patchwork` |
| Patient ordering | Arbitrary | By tumor type → then by dominant signature weight |
| Signature labeling | Short SBS code | Extended legend: SBS code + full aetiology description |
| Tumor-type summary | Not shown | Right panel: mean weight per tumor type |
| Reproducibility | Tool-dependent | `set.seed(11)` · bitwise-identical |

---

## Visual anatomy

```
  [LEFT PANEL — Per-patient stacked 100% bar]
  1.0 ├─────────────────────────────────────────────────
      │ SBS3 SBS13 SBS40 SBS4 SBS1...  (colored stacks)
  0.5 ┤                                                 [RIGHT PANEL]
      │ ← NSCLC (25 pts) → ← BRCA → ← HCC → ← CRC →  Mean weight
  0.0 └─────────────────────────────────────────────────  by tumor type
       Patient (ordered by tumor type, then dominant sig) → ████ NSCLC
                                                            ███  BRCA
  Legend: SBS3 – Homologous recombination deficiency        ██   HCC
          SBS13 – APOBEC activity                           ████ CRC
          SBS4 – Tobacco smoking                            ███  PDAC
          ...
```

| Element | Description |
|---|---|
| Left stacked bars | `geom_col(position="stack")` · width=1 (no gap between patients) |
| Patient ordering | Tumor type → dominant signature weight descending |
| Legend | 3-column layout · SBS code + aetiology · `guide_legend(ncol=4)` |
| Right panel | Mean signature weight per tumor type per signature · `position="stack"` |
| Patchwork layout | `widths = c(13, 4)` |

---

## Output figure

| File | Dimensions | DPI |
|---|---|---|
| `plots/11_signature.png` | 18 × 10 in | 180 |

---

## Synthetic dataset — ONCVIZ-001 ADaM v1 (ADSIG)

```
n patients     80
n signatures   12 (COSMIC SBS v3.3 subset)
Weight sum     1.0 per patient (Dirichlet-like decomposition)
Seed           11 · fully reproducible
```

### Signatures displayed

| SBS | Aetiology | Color |
|---|---|---|
| SBS1 | Spontaneous deamination of 5-methylcytosine (clock-like) | `#f1948a` |
| SBS2 | APOBEC cytidine deaminase activity | `#85c1e9` |
| SBS3 | Homologous recombination deficiency (HRD) | `#3498db` |
| SBS4 | Tobacco smoking | `#e74c3c` |
| SBS5 | Unknown clock-like | `#1abc9c` |
| SBS6 | Mismatch repair deficiency (MMRd) | `#d7bde2` |
| SBS13 | APOBEC cytidine deaminase activity | `#aed6f1` |
| SBS17 | Unknown aetiology (gastric/esophageal) | `#922b21` |
| SBS18 | Oxidative damage (ROS) | `#d5dbdb` |
| SBS22 | Aristolochic acid exposure | `#f9e79f` |
| SBS40 | Unknown clock-like | `#f39c12` |
| SBS44 | Mismatch repair deficiency (treatment-related) | `#8e44ad` |

### Variables used

| Variable | Source | Description |
|---|---|---|
| `USUBJID` | ADSIG / ADSL | Subject identifier |
| `SIG_NAME` | ADSIG | COSMIC SBS signature code |
| `SIG_WEIGHT` | ADSIG | Signature contribution weight (0–1, sum to 1 per patient) |
| `TUMORTYPE` | ADSL (joined) | Histology (for ordering and right panel) |
| `ANL01FL` | ADSIG | Analysis flag (filtered to "Y") |

---

## Key parameters

| Parameter | Value | Description |
|---|---|---|
| Bar width | 1.0 | Contiguous (no gap between patients) |
| Patient ordering | Tumor type → `desc(dominant SIG_WEIGHT)` | Groups histologies, then separates high-burden patients |
| Y-axis limits | `c(0, 1.02)` | 0.02 headroom for visual clarity |
| Right panel bar | `position="stack"` | Stacked mean weights per tumor type |
| Legend columns | 4 | `guide_legend(ncol=4)` |
| Patchwork widths | `c(13, 4)` | Main panel : summary panel |
| Color palette | COSMIC v3.3 standard | See table above |
| Figure dimensions | 18 × 10 in | Wide landscape for 80 patient columns |
| DPI | 180 | Publication quality |

---

## Statistical considerations

- **Weights sum to 1 per patient:** The signature decomposition assumes that all mutations in a patient's tumor are attributable to the 12 displayed signatures. In reality, COSMIC v3.3 contains 79 SBS signatures; limiting to 12 forces the remaining burden to be distributed among the displayed set, potentially inflating weights for biologically absent signatures.
- **Dominant signature:** Defined here as `slice_max(SIG_WEIGHT, n=1)` per patient. In patients with highly mixed signatures (e.g., SBS3 = 0.35, SBS1 = 0.34), the dominant signature is arbitrary and ordering by it produces visual noise.
- **Clinical relevance:** SBS3 (HRD) predicts sensitivity to PARP inhibitors; SBS4 (tobacco) is enriched in NSCLC; SBS6/SBS44 (MMRd) predict immunotherapy response. The right panel's tumor-type breakdown makes these clinical associations visually apparent.

---

## Limitations

- Synthetic weights are generated by normalizing absolute-value normal draws (`abs(rnorm(k))`), which does not perfectly simulate the Dirichlet distribution used in real signature fitting. Weights should be replaced with SigProfiler or deconstructSigs output for real data.
- The 12 signature subset is arbitrary; real decomposition should use all cosmically supported signatures for the relevant sequencing assay (WGS vs. WES vs. targeted panel).
- With 80 patient columns at 18 in width, each column is ~0.225 in wide; individual patient bars may not be distinguishable when rendered at screen resolution. This is acceptable for publication (the individual patient signal is in the right-panel summary).

---

## Requirements

```r
library(ggplot2)    # >= 3.4
library(dplyr)      # >= 1.1
library(patchwork)  # >= 1.2
```

---

## Files

| File | Description |
|---|---|
| `11_signature.R` | Self-contained script · synthetic data generated internally |
| `plots/11_signature.png` | Output figure · 18 × 10 in · 180 DPI |

---

## References

Alexandrov LB, et al. The repertoire of mutational signatures in human cancer. *Nature.* 2020;578(7793):94–101.

COSMIC Mutational Signatures v3.3. https://cancer.sanger.ac.uk/signatures/

Maura F, et al. A practical guide for mutational signature analysis in hematological malignancies. *Nat Commun.* 2019;10(1):2969.

Tate JG, et al. COSMIC: the Catalogue Of Somatic Mutations In Cancer. *Nucleic Acids Res.* 2019;47(D1):D941–D947.
