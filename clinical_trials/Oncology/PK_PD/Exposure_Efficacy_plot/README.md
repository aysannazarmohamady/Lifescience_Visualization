# Exposure-Efficacy Plot

An exposure-response visualization relating Cycle 1 Day 1 AUC₀₋₂₄ to best percent change in tumor size, fit with a sigmoid Emax model and complemented by an exposure-quartile categorical analysis, from a fully synthetic phase I/II oncology trial dataset.

**Dataset:** ONCVIZ-001 · N = 62 (PK/efficacy-evaluable) | **Cutoff:** 05 March 2026 | **Language:** R · ggplot2 · patchwork | **License:** CC BY 4.0

---

## The gap this fills

Exposure-response figures are frequently built with a simple linear regression on log-exposure, which is not the pharmacologically motivated model regulators expect, and often decorate the plot with a legend or parameter box sitting directly on top of the data it's meant to explain.

| | Prior art | This work |
|---|---|---|
| Model | Linear regression on log(AUC) | Sigmoid Emax: `E = E0 + Emax·AUC^γ / (EC50^γ + AUC^γ)` |
| Uncertainty | Analytic CI (assumes linearity) | 500-replicate patient-level bootstrap 95% CI band |
| Model parameters | Often omitted or approximate | E0, Emax, EC50, γ (Hill coefficient) reported exactly in the caption |
| Legend / parameter placement | Frequently overlaps data points | Placed entirely outside the plotted data area |
| Categorical view | Often absent | Companion exposure-quartile panel (median exposure vs. mean response ± CI) |
| Null result handling | Sometimes forced into a misleading trend | Flat/non-significant fit reported as observed, not forced |

---

## Visual anatomy

```
 40% |  ●  ●
     |    ●    ___________________
   0 |  ●   ●_/                       ← sigmoid Emax fit, shaded 95% CI
-30% |. . . . . . . . . . . . . . .   ← PR threshold
     |              ●   ●    ●
-100%|         ●  ●    ●  ●     ●
     3        10       30   (log AUC, ng·h/mL)
```

| Element | Description |
|---|---|
| Points | Individual patients, colored by dose cohort |
| Black curve | Sigmoid Emax model fit (nonlinear least squares) |
| Gray band | 95% CI from 500 bootstrap refits (patient-level resampling) |
| Dotted line | PR threshold (−30%), labeled inside the plot in verified empty space |
| Right panel | Exposure-quartile analysis — median AUC per quartile vs. mean response ± 95% CI |
| Legend | Placed below the plot, never over data |
| Model parameters | Folded into the figure caption, not an in-plot text box |

---

## Model specification

```r
emax_fun <- function(auc, e0, emax, ec50, gamma)
  e0 + emax * (auc^gamma) / (ec50^gamma + auc^gamma)

fit <- nls(BEST_PCHG ~ emax_fun(AUC, e0, emax, ec50, gamma), data = df,
           start = list(e0 = mean(df$BEST_PCHG), emax = -80, ec50 = median(df$AUC), gamma = 1),
           lower = c(-100, -300, 0.1, 0.3), upper = c(100, 300, 100, 8),
           algorithm = "port")
```

95% CI band is derived from 500 bootstrap refits (resampling patients with replacement), not a delta-method approximation, since the nonlinear model's parameter uncertainty is not well approximated analytically at this sample size.

**This dataset's actual result:** the fit is flat and poorly constrained (wide CI on Emax and EC50) — reported honestly as no true exposure-response signal, rather than forced to imply a relationship that isn't supported by the data.

---

## Dataset — ONCVIZ-001 ADaM v1

| Domain | Description | Key variables |
|--------|-------------|---|
| ADPK | AUC₀₋₂₄, Cycle 1 Day 1 | `PARAMCD="AUC"`, `AVAL`, `DOSE` |
| ADTR | Tumor measurements | `AVISITN`, `PCHG` (best post-baseline value used) |

```r
auc  <- adpk |> filter(PARAMCD == "AUC", AVISIT == "CYCLE 1 DAY 1") |> distinct(USUBJID, .keep_all = TRUE)
best <- adtr |> filter(AVISITN > 0) |> group_by(USUBJID) |> summarise(BEST_PCHG = min(PCHG, na.rm = TRUE))
```

---

## Output files

| File | Description | Dimensions |
|------|-------------|------------|
| `exposure_efficacy.R` | Main R script | |
| `Out/exposure_efficacy.png` | Emax fit panel + exposure-quartile panel | 12×6.5 in · 300 DPI |

---

## When to use

**Appropriate:**
- Dose/exposure justification sections of regulatory submissions
- Assessing whether higher exposure is associated with greater efficacy (or a plateau)
- Reporting a null exposure-response result transparently

**Limitations:**
- Emax model requires reasonable exposure range coverage — narrow ranges (e.g., single fixed dose) will not constrain EC50/γ well
- Bootstrap CI assumes patient-level resampling is a valid uncertainty proxy — not a substitute for a full population PK/PD model
- Best % change (not landmark response) is the endpoint here — a time-to-event exposure-response model would be needed for OS/PFS endpoints

---

## Requirements

```
R >= 4.1
```

```r
install.packages(c("dplyr", "ggplot2", "patchwork"))
```

---

## References

U.S. Food and Drug Administration. *Exposure-Response Relationships — Study Design, Data Analysis, and Regulatory Applications.* CDER; 2003.

Sheiner LB. Learning versus confirming in clinical drug development. *Clin Pharmacol Ther.* 1997;61(3):275–291.

Wang Y, et al. Model-based drug development: the road to quantitative pharmacology. *J Pharmacokinet Pharmacodyn.* 2009;36(6):427–427.

Gabrielsson J, Weiner D. *Pharmacokinetic and Pharmacodynamic Data Analysis: Concepts and Applications.* 5th ed. Swedish Pharmaceutical Press; 2016.

