# Forest Plot — Overall Survival Subgroup Analysis (ADaM-driven)

A publication-ready forest plot reporting Cox proportional hazards ratios (HR) with 95% confidence intervals for overall survival across 16 pre-specified subgroups, comparing the treatment arm against control in the ONCVIZ-001 synthetic basket trial, with inverse-variance weight shown as marker size and a companion HR/CI/weight text table.

**Dataset:** ONCVIZ-001 · Treatment (n = 62) vs. Control (n = 18) · ADTTE / ADSL | **Endpoint:** Overall Survival (OS) | **Language:** R · ggplot2 · survival | **License:** CC BY 4.0

---

## The gap this fills

Forest plots are mandatory in regulatory submissions for communicating treatment-effect heterogeneity across subgroups. Most ggplot2 builds omit the weight-proportional marker size that read-through publications use to communicate subgroup precision at a glance, and few pair the graphic with a text column reproducing the exact HR/CI/weight values next to each row.

| | Prior art | This work |
|---|---|---|
| Marker size | Fixed | ∝ inverse-variance weight (1/SE² of log HR), with a size legend |
| Text columns | Not shown, or a separate table | HR (95% CI) and Weight % columns aligned to each row |
| Overall HR symbol | Circle or square | Diamond, dark blue |
| Significant subgroups | Not highlighted | Red when 95% CI excludes HR = 1 |
| X-axis scale | Linear or unlabeled log | Log₁₀ with clinical fraction breaks (0.1, 0.25, 0.5, 1, 2, 4) |
| Cox model | Unstratified pooled | Per-subgroup `coxph(Surv(AVAL, EVENT) ~ TRT)` |
| Reproducibility | Inline synthetic data | Reads real `ADTTE`/`ADSL` from `Data/V1` |

---

## Visual anatomy

```
  Subgroup              HR (95% CI)          Weight %
  ─────────────────────────────────────────────────── HR = 1
  Overall           ◆────────────────           0.51 (0.27–0.98)     —
  NSCLC             ●───────────────            0.36 (0.04–3.12)   1.7%
  Age < 65          ●──── (red, sig)            0.24 (0.09–0.61)   8.9%
  ...
  ─────────────────────────────────────────────────────
   0.1   0.25  0.5    1      2      4
   ← Favors Treatment      Favors Control →
```

| Element | Description |
|---|---|
| Diamond (◆) | Overall pooled HR |
| Circle (●), size-scaled | Subgroup HR, radius ∝ inverse-variance weight |
| Red circle | Subgroup with 95% CI entirely below HR = 1 |
| Dashed vertical line | HR = 1 reference |
| Weight legend | Three reference dot sizes shown below the x-axis |

---

## Dataset variables used

| Variable | Source | Description |
|---|---|---|
| `USUBJID`, `ARM` | ADTTE | Subject identifier, arm |
| `AVAL`, `CNSR` | ADTTE | Time to event, censoring flag |
| `TUMORTYPE`, `AGEGR1`, `SEX`, `ECOG`, `TMBHIGH`, `PDL1GRP` | ADTTE | Subgroup stratification variables |

---

## Statistical method

A separate Cox model is fit per subgroup: `coxph(Surv(AVAL, EVENT) ~ TRT, data = subgroup_df)`, with `EVENT = 1 - CNSR`. HR = exp(coefficient); 95% CI from `confint()`. Inverse-variance weight = `1/SE²` of the log-HR, normalized to sum to 100% across subgroup rows (the Overall summary row is excluded from this normalization, consistent with meta-analysis convention).

**Conditions for valid estimation:** both arms need ≥2 patients with ≥1 event; otherwise the row is marked NE (not estimable).

---

## Key parameters

| Parameter | Value |
|---|---|
| CI method | `confint(coxph_fit)`, profile likelihood on log-HR scale |
| CI x-clip | (0.05, 10) for display only |
| X-axis breaks | 0.1, 0.25, 0.5, 1, 2, 4 |
| Marker size range | 2–9 pt radius, linear in weight % |

---

## Limitations

- **Small control arm (n = 18)** means several subgroup control cells contain only 2–5 patients; CIs are correspondingly wide.
- **Unstratified per-subgroup Cox model** — no adjustment for other covariates within each subgroup fit.
- **No formal interaction test.** A CI excluding HR = 1 in one subgroup does not, by itself, constitute a significant treatment-by-subgroup interaction.
- **Multiple testing.** 16 subgroup estimates are shown without adjustment; treat as exploratory/hypothesis-generating.

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
| `forest_plot.R` | Self-contained script · reads `Data/V1/ADTTE.csv` |
| `Out/forest_plot.png` | Output figure · 11 × 8.3 in · 300 DPI |
| `Out/forest_plot_table.csv` | Companion table: label, n, HR, CI, weight % |

---

## References

Therneau TM, Grambsch PM. *Modeling Survival Data: Extending the Cox Model.* Springer; 2000.

Pocock SJ, Assmann SE, Enos LE, Kasten LE. Subgroup analysis, covariate adjustment and baseline comparisons in clinical trial reporting: current practice and problems. *Stat Med.* 2002;21(19):2917–2930.

ICH E9(R1). *Statistical Principles for Clinical Trials: Addendum on Estimands and Sensitivity Analysis.* 2019.
