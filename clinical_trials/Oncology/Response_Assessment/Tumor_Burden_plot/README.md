# Tumor Burden Plot

A longitudinal visualization of tumor burden trajectories displaying mean/median Sum of Longest Diameters (SLD), percentage change from baseline, and response group dynamics by Best Overall Response and tumor histology, per RECIST 1.1 criteria, across a fully synthetic phase II/III basket trial dataset.

**Dataset:** ONCVIZ-001 · N = 62 (Treatment Arm, evaluable subset) | **Cutoff:** 05 Mar 2026 | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

Tumor burden trajectory plots are a standard component of oncology efficacy readouts, yet most open-source implementations are limited to single-panel, single-grouping designs with minimal clinical annotation. Existing examples from package vignettes typically use toy datasets of 20–50 patients with no calibration to real-world trial benchmarks. Multi-panel layouts that simultaneously show absolute SLD, percentage change, median trajectories, and histology-stratified trends are rarely found as reproducible, fully documented code.

| | Prior art | This work |
|---|---|---|
| Dataset size | 20–50 synthetic patients (most vignettes) | 62 evaluable treatment-arm patients from ONCVIZ-001 |
| Tumor histologies | Single histology | Five independently calibrated profiles (NSCLC, BRCA, CRC, HCC, PDAC) |
| Cross-domain consistency | Isolated example data | Same 400-patient, 13 ADaM domain, 131,690-record dataset as companion plots |
| Panel coverage | Single metric (mean SLD or % change) | Four complementary panels: mean SLD, % change, median SLD, tumor-type stratification |
| Uncertainty encoding | Point estimates only, or simple error bars | 95% CI ribbons (panels A, D) · ±1 SE ribbons (panel B) · IQR ribbons (panel C) |
| Reference lines | None or horizontal zero only | BL mean/median · RECIST PD threshold (+20%) · RECIST PR threshold (−30%) |
| Number at risk table | Not included | Embedded per-group patient counts at every visit (BL through C12) |
| Reproducibility | No fixed seed, or seed not stated | `set.seed(42)` — bitwise-identical across R ≥ 4.1 |
| Calibration traceability | Unspecified | Anchored to KEYNOTE-189, IMbrave150, OlympiAD, KEYNOTE-177, NAPOLI-1 |

---

## Why the tumor burden plot complements the waterfall and swimmer plots

The waterfall plot answers *how much* — depth of tumor response at the single best assessment. The swimmer plot answers *how long* — duration of response and outcome status for each patient. The tumor burden plot answers *how* and *when* — the evolving trajectory of tumor size across all assessment cycles, both as group-level trends and as median/mean summaries that reveal when responses deepen, plateau, or reverse.

| Question | Waterfall | Swimmer | Tumor Burden |
|---|---|---|---|
| Depth of response (% SLD change at nadir) | ✓ | — | — |
| Absolute SLD trajectory over time | — | — | ✓ |
| % Change trajectory over time | — | — | ✓ |
| Duration of response | — | ✓ | — |
| Uncertainty around group means | — | — | ✓ |
| Response by tumor histology | — | — | ✓ |
| Patient counts at each visit | — | — | ✓ |
| RECIST threshold context (+20% / −30%) | — | — | ✓ |

---

## Visual anatomy

```
Panel A — Mean SLD by BOR (95% CI ribbon)
─────────────────────────────────────────
  100 ┤
      │  ···BL mean = 40.3 mm (dashed reference)
   50 ┤  ████ CR (green)   ████ PR (blue)
      │  ████ SD (gold)    ████ PD (red)
    0 ┴──BL──C1──C2──···──C12

Panel B — Mean % Change by BOR (±1 SE ribbon)
───────────────────────────────────────────────
   40 ┤  - - - - - - - +20% PD threshold (red dashed)
    0 ┤  ──────────────────── zero line
  −30 ┤  - - - - - - - −30% PR threshold (green dashed)
  −80 ┴──C1──C2──C3──···──C12

Panel C — Median SLD by BOR (IQR ribbon)
──────────────────────────────────────────
   75 ┤  ···BL median = 34.3 mm (dashed reference)
      │  Dashed lines = median trajectories
    0 ┴──BL──C1──C2──···──C12

Panel D — Mean SLD by Tumor Type (95% CI ribbon)
──────────────────────────────────────────────────
   80 ┤  ████ NSCLC (blue)   ████ BRCA (red)
      │  ████ CRC (purple)   ████ HCC (gold)
      │  ████ PDAC (green)
    0 ┴──BL──C1──C2──···──C12

Each panel includes a Number at risk table beneath the x-axis,
showing per-group patient counts at every visit.
```

| Element | Description |
|---|---|
| Coloured ribbon | Uncertainty band (95% CI, ±1 SE, or IQR depending on panel) |
| Solid/dashed line | Group mean or median SLD trajectory |
| Filled circle | Mean or median value at each assessment visit |
| Grey dashed horizontal | Baseline mean or median SLD reference line |
| Red dashed horizontal | +20% RECIST 1.1 PD threshold (panel B only) |
| Green dashed horizontal | −30% RECIST 1.1 PR threshold (panel B only) |
| Number at risk table | Patient counts per group at each visit, colour-coded by group |
| Panel A | Mean SLD (mm) by Best Overall Response, 95% CI ribbon |
| Panel B | Mean % change from baseline by BOR, ±1 SE ribbon |
| Panel C | Median SLD (mm) by Best Overall Response, IQR ribbon |
| Panel D | Mean SLD (mm) by tumor histology, 95% CI ribbon |

---

## Panel descriptions

### Panel A — Mean SLD by Best Overall Response (95% CI ribbon)

Mean absolute SLD across all cycles from baseline (BL) to Cycle 12 (C12), grouped by RECIST Best Overall Response. The 95% CI ribbon is computed using the t-distribution (`qt(0.975, df = n − 1) × sd / √n`). A dashed grey horizontal line marks the overall baseline mean (40.3 mm). CR and PR groups show progressive SLD decline; SD remains near baseline; PD shows upward drift.

- Y-axis: 0 – 100 mm
- Visits: BL, C1 – C12

### Panel B — Mean % Change in SLD by BOR (±1 SE ribbon)

Mean percentage change from baseline SLD at each on-treatment visit, with ±1 standard error ribbons. Two annotated reference lines provide RECIST 1.1 clinical context: the +20% PD threshold (red dashed) and the −30% PR threshold (green dashed). This panel makes it immediately visible at which cycle each response group crosses the RECIST response or progression boundary on average.

- Y-axis: −85% to +40%
- Visits: C1 – C12 (baseline excluded as % change is undefined at BL)
- Reference lines: +20% PD · −30% PR

### Panel C — Median SLD by Best Overall Response (IQR ribbon)

Median SLD trajectory with interquartile range (IQR) ribbon, providing a distribution-robust alternative to panel A. Uses dashed line style to visually distinguish median from mean (panel A). A dashed grey horizontal line marks the overall baseline median (34.3 mm). IQR ribbons are wider for SD and PD groups, reflecting greater inter-patient variability in tumor size among non-responders.

- Y-axis: 0 – 75 mm
- Visits: BL, C1 – C12
- Ribbon: 25th – 75th percentile

### Panel D — Mean SLD by Tumor Type (95% CI ribbon)

Mean SLD trajectory stratified by tumor histology across all five basket trial indications. Enables direct comparison of baseline tumor burden, trajectory shape, and absolute tumor size by histology. PDAC and HCC show the highest baseline SLD; BRCA shows the most rapid early reduction among tumor types.

- Y-axis: 0 – 80 mm
- Visits: BL, C1 – C12
- Groups: NSCLC · BRCA · CRC · HCC · PDAC

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

### BOR distribution in this figure (evaluable subset, n = 62)

| BOR | n | Colour |
|-----|--:|--------|
| CR | 15 | `#2e8b57` (sea green) |
| PR | 7 | `#4169e1` (royal blue) |
| SD | 14 | `#d4a017` (golden) |
| PD | 26 | `#b22222` (firebrick) |

### Domains used by this script

| Domain | Description | Key variables used |
|--------|-------------|---|
| ADSL | Subject-level | `USUBJID`, `ARM`, `BESTRSPC` |
| ADTR | Tumor results (RECIST) | `PARAMCD="SUMDIAM"`, `ANL01FL`, `AVISITN`, `AVAL`, `PCHG`, `TUMORTYPE` |

```r
# Load required ADaM domains
adtr <- read.csv("ADTR.csv")
adsl <- read.csv("ADSL.csv")

# Subset to treatment arm
trt <- adsl |> filter(ARM == "TREATMENT") |> select(USUBJID, BESTRSPC)

# Filter ADTR to SLD parameter, analysis flag, treatment patients
df <- adtr |>
  filter(PARAMCD == "SUMDIAM", ANL01FL == "Y",
         USUBJID %in% trt$USUBJID) |>
  left_join(trt, by = "USUBJID")
```

---

## Key implementation details

### Uncertainty ribbon computation

Three distinct uncertainty encodings are used across the four panels to serve different analytical purposes:

```r
# Panel A / D — 95% CI using t-distribution
t_ci <- function(x)
  qt(0.975, df = length(x) - 1) * sd(x, na.rm = TRUE) / sqrt(sum(!is.na(x)))

# Panel B — ±1 Standard Error
se <- sd(PCHG, na.rm = TRUE) / sqrt(sum(!is.na(PCHG)))

# Panel C — Interquartile Range (25th / 75th percentile)
lo = quantile(AVAL, 0.25, na.rm = TRUE)
hi = quantile(AVAL, 0.75, na.rm = TRUE)
```

### Number at risk table

The `nar_plot()` helper function generates a visit-by-group patient count table that attaches beneath each main panel via `patchwork`. Counts are the number of distinct `USUBJID` values with an observation at each visit. The table is colour-coded by group to match the trajectory lines.

```r
nar_plot <- function(data, grp_var, grp_levels, grp_colors,
                     vis = 0:12, x_labs = vis_labs) {
  counts <- data |>
    filter(AVISITN %in% vis) |>
    group_by(.data[[grp_var]], AVISITN) |>
    summarise(n = n_distinct(USUBJID), .groups = "drop") |>
    tidyr::complete(!!sym(grp_var) := grp_levels, AVISITN = vis,
                    fill = list(n = 0))
  ...
}
```

### RECIST threshold annotations (panel B)

Reference lines at +20% (PD) and −30% (PR) are drawn as annotated `geom_hline` layers. Labels are offset 0.7 months beyond the last visit to avoid overprinting the data.

```r
geom_hline(yintercept = 20,  colour = "#e74c3c", linetype = "dashed", linewidth = 0.8)
geom_hline(yintercept = -30, colour = "#27ae60", linetype = "dashed", linewidth = 0.8)
annotate("text", x = 12.7, y = 20,  label = "+20%\nPD threshold", ...)
annotate("text", x = 12.7, y = -30, label = "-30%\nPR threshold", ...)
```

### Visit axis encoding

`AVISITN` (numeric visit index 0–12) is mapped to a factor with levels BL, C1–C12 for display while retaining the numeric scale for correct x-axis spacing and break alignment.

```r
VL = factor(
  dplyr::recode(as.character(AVISITN),
    "0"="BL","1"="C1", ... ,"12"="C12"),
  levels = c("BL","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12"))
```

### Layout assembly with patchwork

Each panel is paired with its number-at-risk table using a 1 : 0.22 height ratio (1 : 0.30 for panel D to accommodate five tumour-type rows). The four panel+table units are assembled into a 2 × 2 grid with a shared title, subtitle, and caption.

```r
wrap <- function(p, nar, nar_h = 0.22)
  p / nar + plot_layout(heights = c(1, nar_h))

final <- (wrap(pA, narA) | wrap(pB, narB)) /
         (wrap(pC, narC) | wrap(pD, narD, nar_h = 0.30))
```

---

## Output file

| File | Description | Dimensions |
|------|-------------|------------|
| `tumor_burden_plot.R` | R script — all four panels | |
| `ADTR.csv` | Tumor results dataset | — |
| `ADSL.csv` | Subject-level dataset | 400 rows |
| `tumor_burden_plot_R.png` | Four-panel figure | 22 × 22 in · 180 DPI |

---

## When to use

**Appropriate:**
- Displaying longitudinal SLD trajectories alongside BOR in Phase I/II/III efficacy readouts
- Communicating the timing and depth of group-level tumor response across cycles
- Supplementary figures in clinical publications, commonly paired with waterfall and spider plots
- Regulatory dossiers requiring cycle-by-cycle summary of tumor burden evolution
- Basket trial figures where histology-specific SLD trajectories need to be compared side by side

**Limitations:**
- Shows group-level averages only — individual patient trajectories require spider plots
- Missing visit data reduces the number at risk at later cycles and may bias group means
- Percentage change (panel B) is undefined at baseline and excludes patients with zero baseline SLD
- IQR ribbons (panel C) can be wide for small or heterogeneous groups, limiting interpretability
- A treatment arm of n = 62 evaluable patients limits the precision of 95% CI estimates, particularly for the PR group (n = 7)

---

## Requirements

```
R >= 4.1
ggplot2 >= 3.4
dplyr  >= 1.1
tidyr  >= 1.3
tibble >= 3.2
patchwork >= 1.2
rlang >= 1.1
```

No proprietary packages required. All dependencies are available on CRAN.

---

## References

Eisenhauer EA, et al. New response evaluation criteria in solid tumours: RECIST 1.1. *Eur J Cancer.* 2009;45(2):228–247.

Finn RS, et al. Atezolizumab plus bevacizumab in unresectable hepatocellular carcinoma. *N Engl J Med.* 2020;382(20):1894–1905.

Gandhi L, et al. Pembrolizumab plus chemotherapy in metastatic non-small-cell lung cancer. *N Engl J Med.* 2018;378(22):2078–2092.

André T, et al. Pembrolizumab in microsatellite-instability–high advanced colorectal cancer. *N Engl J Med.* 2020;383(23):2207–2218.

Robson M, et al. Olaparib for metastatic breast cancer in patients with a germline BRCA mutation. *N Engl J Med.* 2017;377(6):523–533.

Therasse P, et al. New guidelines to evaluate the response to treatment in solid tumors. *J Natl Cancer Inst.* 2000;92(3):205–216.

Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

Pedersen TL. patchwork: The Composer of Plots. R package version 1.2.0. https://patchwork.data-imaginist.com
