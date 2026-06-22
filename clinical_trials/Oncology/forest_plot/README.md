# Forest Plot — Overall Survival Subgroup Analysis

A publication-ready forest plot reporting Cox proportional hazards hazard ratios (HR) with 95% confidence intervals for overall survival across 12 pre-specified subgroups, comparing the treatment arm against control in the ONCVIZ-001 synthetic basket trial. The overall HR is rendered as a diamond; subgroups with CI entirely below HR = 1 are highlighted in red.

**Dataset:** ONCVIZ-001 · Treatment (n = 62) vs. Control (n = 18) · ADTTE / ADSL | **Endpoint:** Overall Survival (OS) | **Language:** R · ggplot2 · survival | **License:** CC BY 4.0

---

## The gap this fills

Forest plots are mandatory in regulatory submissions and meta-analyses for communicating treatment effect heterogeneity across patient subgroups. R implementations often fall into two categories: (1) purpose-built packages (`forestplot`, `meta`) with limited ggplot2 integration and rigid styling, or (2) raw `ggplot2` builds that omit the subgroup header rows, the overall diamond symbol, or the log-scale x-axis with arrow annotations.

| | Prior art | This work |
|---|---|---|
| Overall HR symbol | Circle or square | Diamond (`shape = 18`) in dark blue |
| Significant subgroups | Not highlighted | Red point when 95% CI excludes HR = 1 |
| X-axis scale | Linear or log (unlabeled) | Log₁₀ with clinical fraction breaks (0.1, 0.25, 0.5, 1, 2, 4) |
| Arrow annotation | Not shown | Caption: "← Favors Treatment   Favors Control →" |
| Reference line | Dashed at 1 | `geom_vline` at x = 1 |
| CI truncation | Not shown | CI clipped to (0.05, 10) range; arrows in caption |
| Cox model | Unstratified | Per-subgroup `coxph(Surv(AVAL, EVENT) ~ TRT)` |
| Reproducibility | No fixed seed | `set.seed(4)` · bitwise-identical |

---

## Visual anatomy

```
  Subgroup              HR (95% CI)
  ─────────────────────────────────────────────────────── HR = 1
  Overall           ◆─────────────────────────────── (diamond)
  NSCLC             ●────────────────────
  CRC               ●──────────────────────────
  BRCA              ●────────────────
  HCC               ● (small n)
  < 65 yrs          ●──────────────
  ≥ 65 yrs          ●──────────────────────
  Male              ●──────────────
  Female            ●──────────────────
  TMB-High          ●────────────
  TMB-Low           ●────────────────────
  KRAS Mut          ● (red, CI < 1) ─────
  KRAS WT           ●────────────────────
  ─────────────────────────────────────────────────────── 
   0.1   0.25  0.5    1      2      4
```

| Element | Description |
|---|---|
| Diamond (◆) | Overall HR · `shape = 18 · size = 5 · color = #1a3f8f` |
| Gray circle (●) | Subgroup HR · `shape = 16 · size = 2.5` |
| Red circle (●) | Subgroup with 95% CI entirely below 1 · `color = #c0392b` |
| Horizontal error bars | 95% CI from `coxph` · clipped to (0.05, 10) |
| Dashed vertical line | HR = 1 reference |

---

## Output figure

| File | Dimensions | DPI |
|---|---|---|
| `plots/04_forest.png` | 13 × 9 in | 180 |

---

## Synthetic dataset — ONCVIZ-001 ADaM v1 (ADTTE / ADSL)

```
n patients     80 (TRT = 62 · CTL = 18)
Endpoint       Overall Survival (OS)
Event          Death from any cause (EVENT = 1 − CNSR)
Subgroups      12 (tumor type × 3 · age × 2 · sex × 2 · TMB × 2 · KRAS × 2)
Seed           4 · fully reproducible
```

### Variables used

| Variable | Source | Description |
|---|---|---|
| `USUBJID` | ADTTE / ADSL | Subject identifier |
| `ARM` | ADTTE | "TREATMENT" · "CONTROL" |
| `AVAL` | ADTTE | Time to event (months) |
| `CNSR` | ADTTE | Censoring flag (1 = censored · 0 = event) |
| `EVENT` | Derived | `1 − CNSR` |
| `TRT` | Derived | `as.integer(ARM == "TREATMENT")` |
| `TUMORTYPE` | ADSL | Histology |
| `AGEGR1` | ADSL | Age group ("<65" · ">=65") |
| `SEX` | ADSL | "M" · "F" |
| `TMBHIGH` | ADSL | TMB-high flag ("Y" · "N") |
| `KRASMUT` | ADSL | KRAS mutation flag ("Y" · "N") |

---

## Statistical method

A separate **Cox proportional hazards model** is fitted for each subgroup:

```r
coxph(Surv(AVAL, EVENT) ~ TRT, data = subgroup_df)
```

`TRT = 1` for treatment, `TRT = 0` for control. The exponential of the coefficient gives the HR; 95% CI is derived from `confint()` on the log-HR scale.

**Conditions for valid HR estimation:** Both arms must have ≥ 2 patients with ≥ 1 event; otherwise `NA` is returned and no point is drawn. With n = 18 controls, many subgroup control arms are very small — CIs are wide and should be interpreted with caution.

**Proportional hazards assumption:** Not formally tested in this script. In practice, `cox.zph()` should be run on the overall fit before reporting.

---

## Key parameters

| Parameter | Value | Description |
|---|---|---|
| CI method | `confint(coxph_fit)` | Profile likelihood CI on log-HR scale |
| Reference arm | CONTROL | `TRT = 0` |
| Minimum arm size | 2 patients each | Below this, HR returned as NA |
| CI x-clip | (0.05, 10) | `pmax(lo, 0.05)` · `pmin(hi, 10)` |
| X-axis breaks | 0.1, 0.25, 0.5, 1, 2, 4 | Log scale |
| Overall shape | 18 (diamond) | `size = 5` |
| Subgroup shape | 16 (filled circle) | `size = 2.5` |
| Significant color | `#c0392b` | Red · when `hi < 1` |
| Overall color | `#1a3f8f` | Dark blue |
| Non-significant color | `#555555` | Gray |
| Figure dimensions | 13 × 9 in | Portrait, readable subgroup labels |
| DPI | 180 | Publication quality |

---

## Limitations

- **Small control arm (n = 18):** Subgroup control arms may contain as few as 2–4 patients. CIs will be extremely wide and median OS may not be estimable. This is inherent to Phase I/II basket trial designs.
- **Unstratified Cox model:** The per-subgroup model does not stratify by other covariates (tumor type, age). A stratified or multivariable model is needed if confounding is a concern.
- **No interaction test:** A statistically valid test of treatment-by-subgroup interaction (the correct test for effect modification) is not included. Visually, a CI that excludes HR = 1 for a subgroup does *not* imply a significant interaction — this is a common misinterpretation.
- **Multiple testing:** 12 subgroup-level p-values are implicitly presented. No adjustment is applied; subgroup results should be considered exploratory.

---

## Requirements

```r
library(ggplot2)    # >= 3.4
library(dplyr)      # >= 1.1
library(survival)   # >= 3.5
```

---

## Files

| File | Description |
|---|---|
| `04_forest.R` | Self-contained script · synthetic data generated internally |
| `plots/04_forest.png` | Output figure · 13 × 9 in · 180 DPI |

---

## References

Therneau TM, Grambsch PM. *Modeling Survival Data: Extending the Cox Model.* Springer, 2000.

Pocock SJ, Assmann SE, Enos LE, Kasten LE. Subgroup analysis, covariate adjustment and baseline comparisons in clinical trial reporting: current practice and problems. *Stat Med.* 2002;21(19):2917–2930.

Royston P, Altman DG. Visualizing and assessing discrimination in the logistic regression model. *Stat Med.* 2010;29(24):2508–2520.

ICH E9(R1). Statistical Principles for Clinical Trials: Addendum on Estimands and Sensitivity Analysis. 2019.
