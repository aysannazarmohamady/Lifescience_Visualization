# Lollipop Plot — TP53 Somatic Mutations

A publication-ready lollipop plot displaying the position, variant allele frequency (VAF), variant classification, clonality status, and tumor-type origin of somatic point mutations in the **TP53** gene, rendered across all five histologies of the ONCVIZ-001 synthetic basket trial dataset.

**Dataset:** ONCVIZ-001 · N = 38 TP53 mutations across 38 patients | **Gene:** TP53 (17p13.1) | **Language:** R · ggplot2 · ggrepel | **License:** CC BY 4.0

---

## The gap this fills

Lollipop plots are the canonical representation for visualizing the distribution of somatic mutations along a protein sequence, yet most implementations in the oncology literature either (1) rely on external tools such as MutationMapper or cBioPortal that require internet connectivity and produce fixed-style outputs, or (2) produce minimal R plots that omit protein domain annotations, clonality encoding, tumor-type stratification, and high-VAF labels simultaneously.

| | Prior art | This work |
|---|---|---|
| Domain annotation | Separate track or omitted | Inline colored bar below x-axis with domain label |
| Clonality encoding | Not shown | Filled circle = clonal; open circle = subclonal |
| Tumor-type origin | Omitted | Colored tick marks below the domain bar |
| VAF-based labeling | Not shown | `ggrepel` labels at VAF ≥ 0.55 (HGVSp short form) |
| Color scheme | Default ggplot2 | Consistent ONCVIZ-001 palette shared across all 14 figures |
| Reproducibility | No fixed seed | `set.seed(1)` · bitwise-identical output |

---

## Visual anatomy

```
  VAF
  0.80 ──────────────────────────────────────────────────
       │         ●  R175H                 ○
  0.60 │    ●         ○          ●
       │  ●                ○  ●
  0.40 │       ○    ●         ○     ●    ○
       │  |    |    |    |    |    |    |    |   (stems)
  0.00 ├────────────────────────────────────────────────
       ████████████████████████ DNA-Binding Domain ████
       | | |  | |  | tumor-type ticks
        0   200  400  600  800 1000
             Amino Acid Position (codon)
```

| Element | Description |
|---|---|
| Vertical stem | `geom_segment` from y=0 to y=VAF, colored by variant class |
| Filled circle (●) | Clonal mutation (`CLONAL == "Y"`) |
| Open circle (○) | Subclonal mutation (`CLONAL == "N"`, `shape=1, stroke=1.2`) |
| Blue rectangle | DNA-Binding Domain annotation (codons 100–1100) |
| Colored ticks below domain | Tumor-type origin (same palette as all ONCVIZ-001 figures) |
| `ggrepel` labels | HGVSp short form for VAF ≥ 0.55 mutations |

---

## Output figure

| File | Dimensions | DPI |
|---|---|---|
| `plots/01_lollipop.png` | 16 × 7.5 in | 180 |

---

## Synthetic dataset — ONCVIZ-001 ADaM v1 (ADMUT subset)

```
n mutations    38
Gene           TP53 (17p13.1)
Histologies    NSCLC · BRCA · HCC · CRC · PDAC
Seed           1 · fully reproducible
VAF range      0.08 – 0.78
```

### Variables used

| Variable | Type | Description |
|---|---|---|
| `USUBJID` | character | Subject identifier |
| `PROTEIN_POS` | integer | Amino acid position (1–1100) |
| `VAF` | numeric | Variant allele frequency (0–1) |
| `VARIANT_CLASS` | character | Missense_Mutation · Nonsense_Mutation · Frame_Shift_Del · Splice_Site |
| `CLONAL` | character | "Y" = clonal · "N" = subclonal |
| `TUMORTYPE` | character | NSCLC · BRCA · HCC · CRC · PDAC |
| `HGVSP` | character | Protein-level change (HGVS notation, e.g. p.R175H) |

---

## Key parameters

| Parameter | Value | Description |
|---|---|---|
| VAF label threshold | 0.55 | Mutations with VAF ≥ 0.55 receive `ggrepel` labels |
| Domain bar | codons 100–1100 | Represents the TP53 DNA-Binding Domain (exons 4–8) |
| Clonal shape | 16 (filled) | `geom_point` for `CLONAL == "Y"` |
| Subclonal shape | 1 (open) | `geom_point` for `CLONAL == "N"`, `stroke = 1.2` |
| Stem alpha | 0.7 | Reduces visual clutter at dense regions |
| Color: Missense | `#2c5f8a` | Blue |
| Color: Nonsense | `#c0392b` | Red |
| Color: Frame Shift Del | `#e67e22` | Orange |
| Color: Splice Site | `#8e44ad` | Purple |
| TT Color: NSCLC | `#2c5f8a` | Blue (tick marks) |
| TT Color: BRCA | `#c0392b` | Red (tick marks) |
| TT Color: HCC | `#27ae60` | Green (tick marks) |
| TT Color: CRC | `#e67e22` | Orange (tick marks) |
| TT Color: PDAC | `#8e44ad` | Purple (tick marks) |
| Figure dimensions | 16 × 7.5 in | Landscape, optimized for journal two-column |
| DPI | 180 | Publication quality |

---

## Statistical considerations

- **VAF as a clonality proxy:** A VAF ≥ 0.40 is commonly used as an empirical threshold for clonal mutations in diploid tumor regions assuming ~100% tumor purity; subclonal mutations (VAF < 0.40) typically arise later in tumor evolution. This plot uses the `CLONAL` flag pre-annotated in ADMUT rather than a hard VAF cutoff.
- **Multiple mutations per patient:** The data permits at most one row per patient in this figure because the synthetic dataset generates one TP53 event per subject; in real ADMUT data with multi-hit genes, a priority hierarchy (e.g., nonsense > missense) should be applied before plotting.
- **Protein position as a surrogate coordinate:** Lollipop plots use amino acid position rather than genomic coordinate, which conflates exon boundaries but is the standard for protein-centric mutation visualization.

---

## Limitations

- Domain annotation is hard-coded to TP53; adapting to another gene requires updating `HUGO_SYMBOL`, domain coordinates, and domain label.
- Clonality is pre-annotated in the synthetic dataset; in real data, clonality inference requires computational tools (e.g., FACETS, PureCN, or PyClone).
- The `ggrepel` label placement is non-deterministic across R sessions unless a seed is set *inside* `geom_text_repel()` via `seed=`; set `seed=42` in `geom_text_repel()` for fully reproducible label positions.

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
| `01_lollipop.R` | Self-contained script · synthetic data generated internally |
| `plots/01_lollipop.png` | Output figure · 16 × 7.5 in · 180 DPI |

---

## References

Freed-Pastor WA, Prives C. Mutant p53: one name, many proteins. *Genes Dev.* 2012;26(12):1268–1286.

Joerger AC, Fersht AR. Structural Biology of the Tumor Suppressor p53. *Annu Rev Biochem.* 2008;77:557–579.

Alexandrov LB, et al. Signatures of mutational processes in human cancer. *Nature.* 2013;500(7463):415–421.

Wickham H. *ggplot2: Elegant Graphics for Data Analysis.* Springer-Verlag New York, 2016.

Slowikowski K. ggrepel: Automatically Position Non-Overlapping Text Labels. R package version 0.9.5. https://CRAN.R-project.org/package=ggrepel
