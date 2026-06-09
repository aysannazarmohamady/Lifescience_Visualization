# Response Heatmap

A dose-cohort × tumor-type heatmap visualizing Objective Response Rate (ORR) across two categorical dimensions simultaneously, allowing rapid identification of which patient subgroups respond best to treatment, with per-cell response breakdowns (CR/PR/SD/PD) and patient counts, across a fully synthetic phase II/III basket trial dataset.

**Dataset:** ONCVIZ-001 · N = 263 (Treatment Arm) | **Cutoff:** 05 Mar 2026 | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

ORR heatmaps are common in oncology dose-escalation publications, yet most open-source implementations are built on toy datasets with a single tumor type and no calibration to real trial benchmarks. Multi-histology examples are rare, and those that exist rarely encode response breakdowns at the cell level.

| | Prior art | This work |
|---|---|---|
| Dataset size | 15–50 synthetic patients, single histology | 263 treatment-arm patients from ONCVIZ-001 |
| Tumor histologies | Single histology or uncalibrated mix | Five independently calibrated profiles |
| Cross-domain consistency | Heatmap data unrelated to other plot types | Same 400 patients, 13 ADaM domains, 131,690 records |
| Cell content | ORR % only | ORR % + CR/PR/SD/PD counts + n= per cell |
| Empty cells | Dropped or zeroed | Explicitly rendered as grey tiles |
| Text contrast | Fixed color | Automatic white/dark switching at ORR ≥ 65% or = 0% |
| Axis ordering | Alphabetical or arbitrary | Dose-ascending y-axis · fixed tumor-type x-order |
| Reproducibility | No fixed seed | `set.seed(42)` — bitwise-identical across R ≥ 4.1 |
| Calibration traceability | Unspecified | Anchored to KEYNOTE-189, IMbrave150, OlympiAD, KEYNOTE-177, NAPOLI-1 |

---

## Why the heatmap complements the waterfall and swimmer plots

The waterfall plot answers *how much* — depth of tumor response at best assessment. The swimmer plot answers *how long* — duration and outcome per patient. The heatmap answers *where* — which subgroup combination yields the highest response rate. Used together from the same dataset they provide the complete RECIST efficacy picture.

| Question | Waterfall | Swimmer | Heatmap |
|---|---|---|---|
| Depth of response (% SLD change) | ✓ | — | — |
| Duration of response | — | ✓ | — |
| ORR by cohort and tumor type | — | — | ✓ |
| Response category breakdown per subgroup | — | — | ✓ |
| Patient-level outcome (ongoing / death) | — | ✓ | — |
| Dose-response relationship across histologies | — | — | ✓ |
| Responders vs. non-responders at a glance | ✓ | ✓ | ✓ |

---

## Visual anatomy

```
              NSCLC     BRCA      CRC       HCC       PDAC
             ┌─────────┬─────────┬─────────┬─────────┬─────────┐
RP2D         │   47%   │   38%   │   13%   │   27%   │   11%   │
(300 mg)     │CR:2 PR:8│         │         │         │         │
             │  n=21   │         │         │         │         │
             ├─────────┼─────────┼─────────┼─────────┼─────────┤
Dose Level 2 │   33%   │   n/a   │   ...   │         │         │
(200 mg)     │         │ (grey)  │         │         │         │
             └─────────┴─────────┴─────────┴─────────┴─────────┘
```

| Element | Description |
|---|---|
| Tile color | ORR (%) — deep blue = 0%, white = ~50%, deep red = 100% |
| Large bold number | ORR % for that cell |
| CR/PR/SD/PD line | Response category counts per cell |
| `n=` | Total patients in that cell |
| Grey tile · `n/a` | No patients enrolled in that cohort × tumor type |
| `n/a` on colored tile | Cell has patients but ORR is undefined |
| White text | ORR ≥ 65% or ORR = 0% (contrast against deep color) |
| Dark text | Mid-range ORR (0% < ORR < 65%), color `#222222` |
| Color bar | Vertical legend, 0–100% with five anchor labels |

### Color scale

| Stop | Hex | ORR |
|---|---|---|
| Deep blue | `#2166AC` | 0% |
| Light blue | `#92C5DE` | 30% |
| White | `#F7F7F7` | 50% |
| Light red | `#F4A582` | 70% |
| Deep red | `#B2182B` | 100% |

---

## Plot variants

### `bor_heatmap.png` — All treatment-arm patients

All 263 treatment-arm patients summarized across four dose cohorts and five tumor types in a single 4 × 5 grid. Each cell displays ORR, response breakdown, and patient count. Cohorts are ordered dose-ascending (bottom to top); tumor types follow the fixed order NSCLC → BRCA → CRC → HCC → PDAC.

- Figure dimensions: 14 × 9 in @ 150 DPI

---

## Dataset — ONCVIZ-001 ADaM v1

```
Study    ONCVIZ-001 · Vizatinib 300 mg QD vs Placebo
Design   Phase II/III Open-Label Randomized Basket Trial
N        400  (TRT = 263 · CTL = 137 · 2:1 ratio)
Records  131,690 across 13 ADaM domains
Seed     42 · fully reproducible
Cutoff   05 Mar 2026
```

### Tumor-stratified response parameters (treatment arm)

| Histology | n | ORR TRT | ORR CTL | OS median | Calibration source |
|---|---|---:|---:|---:|---|
| NSCLC | 88 | 47% | 20% | 19.6 m | KEYNOTE-189 (Gandhi et al., NEJM 2018) |
| HCC | 69 | 27% | 12% | 12.4 m | IMbrave150 (Finn et al., NEJM 2020) |
| CRC | 80 | 13% | 0% | 7.3 m | KEYNOTE-177 (André et al., NEJM 2020) |
| BRCA | 82 | 38% | 12% | 19.8 m | OlympiAD (Robson et al., NEJM 2017) |
| PDAC | 81 | 11% | 0% | 3.6 m | NAPOLI-1 / PRODIGE-4 |

### Domains used by this script

| Domain | Description | Rows | Key variables used |
|---|---|---:|---|
| ADSL | Subject-level | 400 | `USUBJID`, `ARM`, `TUMORTYPE`, `DOSELEVEL`, `BESTRSPC` |

```r
adsl <- read.csv(file.path(DATA_DIR, "ADSL.csv"), stringsAsFactors = FALSE)

trt <- adsl |> filter(ARM == "TREATMENT") |>
  mutate(COHORT = case_when(
    DOSELEVEL == 100 ~ "Dose Level 1\n(100 mg)",
    DOSELEVEL == 200 ~ "Dose Level 2\n(200 mg)",
    DOSELEVEL == 300 ~ "RP2D\n(300 mg)",
    DOSELEVEL == 400 ~ "Dose Level 4\n(400 mg)"
  ))
```

---

## Key implementation details

### Full grid expansion

All cohort × tumor type combinations are pre-generated via `expand.grid()` before joining to the summarized data. This ensures empty cells render as grey tiles rather than disappearing from the plot.

```r
cell_df <- expand.grid(COHORT = COHORT_ORDER, TUMORTYPE = TUMOR_ORDER,
                       stringsAsFactors = FALSE) |>
  left_join(
    trt |> group_by(COHORT, TUMORTYPE) |>
      summarise(n = n(), n_cr = sum(BESTRSPC == "CR"), ...),
    by = c("COHORT","TUMORTYPE")
  ) |>
  mutate(n = replace_na(n, 0L), ...)
```

### Text contrast logic

Text color switches automatically based on ORR value. Tiles at the extremes of the diverging palette (ORR ≥ 65% or ORR = 0%) receive white text; mid-range tiles use dark grey (`#222222`).

```r
txt_color = ifelse(
  !is.na(orr) & (orr >= 65 | orr == 0),
  "white", "#222222"
)
```

### Three-label cell layout

Each populated cell contains three stacked text layers placed with `nudge_y` offsets relative to the tile center:

```r
# ORR % — bold, top
geom_text(aes(label = orr_lbl), fontface = "bold", size = 5.8, nudge_y =  0.12)
# CR/PR/SD/PD counts — small, middle
geom_text(aes(label = counts_lbl),                 size = 2.8, nudge_y = -0.10)
# n= — small, bottom
geom_text(aes(label = n_lbl),                      size = 3.0, nudge_y = -0.28)
```

### Axis ordering

The y-axis is reversed so that Dose Level 1 appears at the bottom and the highest dose at the top, matching the conventional ascending readout of a dose-escalation table.

```r
COHORT    = factor(COHORT,    levels = rev(COHORT_ORDER)),
TUMORTYPE = factor(TUMORTYPE, levels = TUMOR_ORDER)
```

---

## Output files

| File | Description | Dimensions |
|---|---|---|
| `bor_heatmap.R` | Main R script | |
| `ADSL.csv` | Subject-level dataset | 400 rows |
| `bor_heatmap.png` | ORR heatmap | 14 × 9 in · 150 DPI |

---

## When to use

**Appropriate:**
- Displaying ORR across dose cohorts and tumor types in Phase I/II efficacy readouts
- Communicating dose-response patterns across histologies to clinical and regulatory audiences
- Identifying subgroups with consistently high or low ORR at a glance
- Regulatory dossiers and CSR appendices requiring cross-tabulated response summaries
- Basket trial figures where histology-specific ORR patterns need side-by-side comparison
- Paired with waterfall or swimmer plots from the same dataset for a complete efficacy picture

**Limitations:**
- Does not show magnitude of tumor shrinkage — use waterfall plots
- Does not show duration of response — use swimmer plots
- ORR rounding to integer values can collapse near-threshold differences (e.g. 49% vs 51%)
- Cells with very small n carry high uncertainty but are visually indistinguishable from well-powered cells; consider adding confidence interval annotations for cells with n < 5
- A 4 × 5 grid is near the practical density limit; larger grids (more cohorts or histologies) may require faceting or a separate per-histology view

---

## Requirements

```
R >= 4.1
ggplot2 >= 3.4
dplyr  >= 1.1
tidyr  >= 1.3
```

No proprietary packages required. All dependencies are available on CRAN.

---

## References

André T, et al. Pembrolizumab in microsatellite-instability–high advanced colorectal cancer. *N Engl J Med.* 2020;383(23):2207–2218.

Eisenhauer EA, et al. New response evaluation criteria in solid tumours: RECIST 1.1. *Eur J Cancer.* 2009;45(2):228–247.

Finn RS, et al. Atezolizumab plus bevacizumab in unresectable hepatocellular carcinoma. *N Engl J Med.* 2020;382(20):1894–1905.

Gandhi L, et al. Pembrolizumab plus chemotherapy in metastatic non-small-cell lung cancer. *N Engl J Med.* 2018;378(22):2078–2092.

Robson M, et al. Olaparib for metastatic breast cancer in patients with a germline BRCA mutation. *N Engl J Med.* 2017;377(6):523–533.

Seymour L, et al. iRECIST: guidelines for response criteria for use in trials testing immunotherapeutics. *Lancet Oncol.* 2017;18(3):e143–e152.

Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
