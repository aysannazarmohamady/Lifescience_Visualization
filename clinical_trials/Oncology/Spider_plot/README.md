# Spider Plot

A patient-level visualization displaying individual tumor burden trajectories over time as percent change from baseline in sum of longest diameters (SLD), colored by best overall response or tumor type, across a fully synthetic phase II/III basket trial dataset.

**Dataset:** ONCVIZ-001 · N = 263 (Treatment Arm) | **Cutoff:** 05 Mar 2026 | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

Spider plots are standard in oncology publications, yet most open-source implementations share the same structural limitations: small toy datasets, single tumor histology, and no stratification by biomarker. Purpose-built packages solve grammar elegantly but ship with data that behaves nothing like a real trial.

| | Prior art | This work |
|---|---|---|
| Dataset size | 10–30 synthetic patients (typical vignettes) | 263 treatment-arm patients from ONCVIZ-001 |
| Tumor histologies | Single histology or uncalibrated mix | Five independently calibrated profiles |
| Cross-domain consistency | Spider and waterfall from different studies | Same 400 patients, 13 ADaM domains, 131,690 records |
| Stratification | BOR colour only | Panel A: by tumor type (colour = BOR) · Panel B: by TMB status (colour = tumor type) |
| Per-panel statistics | Rarely shown | N, CR count, PR count, ORR% in every panel subtitle |
| Evaluability criteria | Baseline only | Baseline SLD > 0 **and** ≥ 2 post-baseline visits |
| BOR derivation | Usually hard-coded | Derived from ADRS using CR > PR > SD > PD > NE hierarchy |
| Clip value | Varies (+100% to +200%) | Consistently +100%; annotated in caption |
| Reproducibility | No fixed seed, or seed not stated | `set.seed(42)` — bitwise-identical across R ≥ 4.1 |
| Calibration traceability | Unspecified | Anchored to KEYNOTE-189, IMbrave150, OlympiAD, KEYNOTE-177, NAPOLI-1 |

---

## Why the spider plot complements the waterfall and swimmer plots

The three standard oncology response plots answer different questions from the same data. The waterfall plot answers *how much* — the maximum depth of response. The swimmer plot answers *how long* and *what happened next* — the temporal arc through treatment, follow-up, and outcome. The spider plot answers *how did it evolve* — the continuous trajectory of tumor burden at every assessment timepoint. Used together from the same dataset they deliver the complete RECIST efficacy picture that none can provide alone.

| Question | Waterfall | Spider | Swimmer |
|---|---|---|---|
| Depth of response (% SLD change) | ✓ | ✓ (at each visit) | — |
| Tumor burden trajectory over time | — | ✓ | — |
| Time to first response | — | ✓ | ✓ |
| Duration of response | — | ✓ (approximate) | ✓ (exact) |
| Post-treatment follow-up | — | — | ✓ |
| Outcome (ongoing / death) | — | — | ✓ |
| Delayed or sustained responders | — | ✓ | ✓ |
| Responders vs. non-responders at a glance | ✓ | ✓ | ✓ |

---

## Visual anatomy

```
  %
+100 ──────────────────────────────────── (clip ceiling)

  +20 ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄  +20% PD threshold
   0  ──────────────────────────────────  baseline
  -30 ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄  -30% PR threshold

-100
       0    6   12   18   24  ...  months
```

| Element | Description |
|---|---|
| Each line | One patient's tumor burden trajectory from baseline |
| Endpoint dot | Last available on-study measurement per patient |
| Y = 0% | Baseline (no change from treatment start) |
| `--` at −30% | PR threshold (≥30% reduction per RECIST 1.1) |
| `--` at +20% | PD threshold (≥20% increase per RECIST 1.1) |
| Values clipped at +100% | Extreme progressors capped for readability; noted in caption |
| Line colour (Panel A) | Best overall response (CR · PR · SD · PD · NE) |
| Line colour (Panel B) | Tumor type (NSCLC · CRC · HCC · PDAC · BRCA) |
| Panel subtitle | N, CR count, PR count, ORR% for that panel's population |
| Footer | "RECIST 1.1 · Lines clipped at +100% · Data cutoff: 05 Mar 2026" |

---

## Output panels

Two figures are produced, each covering the **Treatment arm** only:

| Output file | Stratification | Panels | Colour by |
|---|---|---|---|
| `spider_plot_A_tumor_type.png` | One panel per tumor type | 5 (NSCLC · CRC · HCC · PDAC · BRCA) | Best overall response |
| `spider_plot_B_tmb.png` | TMB status | 2 (TMB-High · TMB-Low) | Tumor type |

Each panel subtitle reports **N**, number of CRs, number of PRs, and **ORR (%)** for the patients in that panel.

- Figure dimensions (Panel A): 22 × 8.2 in @ 300 DPI
- Figure dimensions (Panel B): 15 × 8.2 in @ 300 DPI

---

## Evaluability criteria

A patient is included if they have **both**:

- A valid baseline sum of longest diameters (`AVAL > 0` at `AVISITN = 0`)
- At least **two** post-baseline tumor measurements (`AVISITN > 0`)

Patients who meet the baseline criterion but have only one post-baseline visit are excluded because a single trajectory segment gives no information about response durability.

---

## Best overall response derivation

Overall response is taken from `ADRS` where `PARAMCD = "OVRLRESP"`. The best response is selected as the record with the lowest rank in the hierarchy:

```
CR = 1  >  PR = 2  >  SD = 3  >  PD = 4  >  NE = 5
```

If no qualifying record exists for a patient, they are assigned **NE**. This derivation is applied dynamically from ADRS rather than read from ADSL, ensuring consistency with the response dataset regardless of pre-derived BOR fields.

```r
rank_map <- c(CR=1, PR=2, SD=3, PD=4, NE=5)

best_resp_from_adrs <- function(uid) {
  rs <- adrs[adrs$USUBJID == uid & adrs$PARAMCD == "OVRLRESP" &
               adrs$AVALC %in% names(rank_map), ]
  if (nrow(rs) == 0) return("NE")
  rs$AVALC[which.min(rank_map[rs$AVALC])]
}
```

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
|-----------|---|--------:|--------:|----------:|---|
| NSCLC | 88 | 47% | 20% | 19.6 m | KEYNOTE-189 (Gandhi et al., NEJM 2018) |
| HCC | 69 | 27% | 12% | 12.4 m | IMbrave150 (Finn et al., NEJM 2020) |
| CRC | 80 | 13% | 0% | 7.3 m | KEYNOTE-177 (André et al., NEJM 2020) |
| BRCA | 82 | 38% | 12% | 19.8 m | OlympiAD (Robson et al., NEJM 2017) |
| PDAC | 81 | 11% | 0% | 3.6 m | NAPOLI-1 / PRODIGE-4 |

### Domains used by this script

| Domain | Description | Rows | Key variables used |
|--------|-------------|-----:|---|
| ADSL | Subject-level | 400 | `USUBJID`, `ARM`, `TUMORTYPE`, `TMBHIGH` |
| ADTR | Tumor measurement | — | `USUBJID`, `AVISITN`, `AVAL`, `ADTN`, `PCHG` |
| ADRS | Overall response per RECIST 1.1 | 1,211 | `PARAMCD="OVRLRESP"`, `AVALC` (CR/PR/SD/PD/NE) |

```r
# Load required ADaM domains
adsl <- read.csv(file.path(DATA_DIR, "ADSL.csv"), stringsAsFactors = FALSE)
adtr <- read.csv(file.path(DATA_DIR, "ADTR.csv"), stringsAsFactors = FALSE)
adrs <- read.csv(file.path(DATA_DIR, "ADRS.csv"), stringsAsFactors = FALSE)

# Subset to treatment arm
trt_ids  <- adsl$USUBJID[adsl$ARM == "TREATMENT"]
adtr_trt <- adtr[adtr$USUBJID %in% trt_ids, ]
```

---

## Key implementation details

### Days-to-months conversion

All time axes use a fixed conversion factor of **30.4375 days/month** (365.25 / 12), consistent with standard ADaM practice. The x-axis upper limit is rounded up to the nearest 3-month boundary and padded by 3 months to accommodate the threshold annotation labels.

```r
DAY2MO        <- 30.4375
adtr_trt$months <- adtr_trt$ADTN / DAY2MO

xmax_data <- max(spider$months, na.rm = TRUE)
XMAX      <- ceiling(xmax_data / 3) * 3 + 3
```

### Clip ceiling

Post-baseline PCHG values are clipped to **+100%** before plotting. The clip is applied to the plotting dataset only; the underlying ADTR values are unmodified. The caption explicitly states the clip value so readers are not misled about the rate of rapid progression.

```r
post_valid$PCHG <- pmin(post_valid$PCHG, CLIP_TOP, na.rm = TRUE)  # CLIP_TOP = 100
```

### Baseline anchoring

Every evaluable patient receives an explicit baseline row at `months = 0`, `PCHG = 0.0` before their post-baseline measurements are appended. This anchors all trajectories at the origin regardless of whether the source data contains a baseline PCHG record.

```r
bl_rows     <- data.frame(USUBJID = valid, months = 0.0, PCHG = 0.0)
spider_base <- bind_rows(bl_rows, post_valid)
```

### Endpoint dot

The last available assessment per patient is highlighted as a filled point using `slice_max(months)`. This makes it immediately clear where each trajectory ends — useful for distinguishing patients who progressed early from those with long stable disease.

```r
geom_point(
  data = sub_df |>
    group_by(USUBJID) |>
    slice_max(months, n = 1, with_ties = FALSE) |>
    ungroup(),
  aes(x = months, y = PCHG),
  size = 1.6, alpha = 0.88
)
```

### Per-panel ORR statistics

Each panel subtitle is computed dynamically from the ADRS-derived BOR for the patients visible in that panel, not from a pre-calculated summary table. This ensures consistency if the evaluable population changes.

```r
n     <- n_distinct(sub_df$USUBJID)
n_cr  <- sum(meta$BESTRESP[meta$USUBJID %in% unique(sub_df$USUBJID)] == "CR")
n_pr  <- sum(meta$BESTRESP[meta$USUBJID %in% unique(sub_df$USUBJID)] == "PR")
orr   <- round((n_cr + n_pr) / n * 100)
stats <- sprintf("N=%d   CR:%d  PR:%d  ORR=%d%%", n, n_cr, n_pr, orr)
```

### Dynamic x-axis

The x-axis upper limit is computed from the data ceiling, so panels with shorter follow-up (e.g. PDAC) do not show unnecessary white space, while panels with longer trajectories (e.g. NSCLC, BRCA) extend automatically.

---

## Key parameters

| Parameter | Value | Description |
|---|---|---|
| `DAY2MO` | 30.4375 | Days-to-months conversion factor |
| `PR_TH` | −30 | Partial response threshold (%) |
| `PD_TH` | +20 | Progressive disease threshold (%) |
| `CLIP_TOP` | +100 | Upper clip value for PCHG before plotting (%) |
| `TUMOR_ORDER` | NSCLC, CRC, HCC, PDAC, BRCA | Left-to-right panel display order in Figure A |
| Min post-baseline visits | 2 | Evaluability threshold |

---

## When to use

**Appropriate:**
- Displaying the longitudinal tumor burden trajectory for each patient alongside a waterfall plot in Phase I/II efficacy readouts
- Identifying patients with delayed responses, pseudoprogression, or prolonged stable disease that a waterfall plot would misclassify
- Stratifying trajectories by biomarker (e.g. TMB, PD-L1) to support exploratory efficacy analyses
- Supplementary figures in clinical publications and regulatory dossiers requiring patient-level time course data

**Limitations:**
- Does not show absolute SLD values — use a raw spaghetti plot of SLD alongside for an unbiased view of variability (Mercier et al., 2019)
- Percentage change from baseline can appear optimistic when baseline SLD is small; unequal variance across timepoints is hidden
- Trajectories become hard to follow beyond ~40–50 patients per panel due to line crossing; consider per-histology or subgroup panels
- Clipping at +100% obscures the true magnitude of rapid progression for extreme progressors
- Does not encode outcome (death, discontinuation) or post-treatment follow-up — use swimmer plot for that
- Does not encode treatment duration as a distinct bar — depth and time are entangled in the trajectory line

---

## Files

| File | Description | Dimensions |
|---|---|---|
| `spider_plot.R` | Main R script — all plots | |
| `Data/V1/ADSL.csv` | Subject-level dataset | 400 rows |
| `Data/V1/ADTR.csv` | Tumor measurement dataset | — |
| `Data/V1/ADRS.csv` | Overall response assessments | 1,211 rows |
| `Outputs/spider_plot_A_tumor_type.png` | Spider plot stratified by tumor type | 22 × 8.2 in · 300 DPI |
| `Outputs/spider_plot_B_tmb.png` | Spider plot stratified by TMB status | 15 × 8.2 in · 300 DPI |

---

## Requirements

```r
library(dplyr)      # >= 1.1
library(ggplot2)    # >= 3.4
library(patchwork)  # >= 1.2
library(gridExtra)  # for grid.arrange
library(grid)       # for textGrob / gpar
```

No proprietary packages required. All dependencies are available on CRAN.

---

## References

André T, et al. Pembrolizumab in microsatellite-instability–high advanced colorectal cancer. *N Engl J Med.* 2020;383(23):2207–2218.

Eisenhauer EA, et al. New response evaluation criteria in solid tumours: RECIST 1.1. *Eur J Cancer.* 2009;45(2):228–247.

Finn RS, et al. Atezolizumab plus bevacizumab in unresectable hepatocellular carcinoma. *N Engl J Med.* 2020;382(20):1894–1905.

Gandhi L, et al. Pembrolizumab plus chemotherapy in metastatic non-small-cell lung cancer. *N Engl J Med.* 2018;378(22):2078–2092.

Mercier F, et al. From waterfall plots to spaghetti plots in early oncology clinical development. *Pharm Stat.* 2019;18(5):526–534.

Robson M, et al. Olaparib for metastatic breast cancer in patients with a germline BRCA mutation. *N Engl J Med.* 2017;377(6):523–533.

Seymour L, et al. iRECIST: guidelines for response criteria for use in trials testing immunotherapeutics. *Lancet Oncol.* 2017;18(3):e143–e152.

Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
