# Survival Curves

A suite of four Kaplan–Meier survival curve figures — Overall Survival (OS), Progression-Free Survival (PFS), Event-Free Survival (EFS), and Disease-Free Survival (DFS) — rendered to publication standard with risk tables, hazard ratio insets, landmark rate annotations, and a custom `theme_oncoviz` theme, across a fully synthetic phase I/II basket trial dataset.

**Dataset:** ONCVIZ-001 · N = 80 (TRT = 62 · CTL = 18) | **Cutoff:** 05 Mar 2026 | **Language:** R · survival · survminer · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

Kaplan–Meier plots are the single most common figure in oncology clinical trial publications, yet the open-source landscape for producing publication-ready KM figures has a consistent set of structural gaps. The `survminer` package (`ggsurvplot`) solves the grammar problem elegantly and is the de-facto standard in academic R workflows, but its canonical vignettes ship with toy datasets — typically the `lung` or `colon` datasets from the `survival` package — that are too small, single-endpoint, and unrepresentative of a real trial's analytical demands. Purpose-built clinical examples are often single-arm, single-endpoint, and divorced from any consistent cross-domain data standard.

| | Prior art | This work |
|---|---|---|
| Dataset | `lung` / `colon` (survival package), 10–30-patient vignettes | ONCVIZ-001: 80 patients, 13 ADaM domains, 131,690 records |
| Endpoints covered | One endpoint per script | Four endpoints in a single script: OS · PFS · EFS · DFS |
| Cross-domain consistency | Each figure from a different dataset | All curves from the same 400-patient ADaM trial |
| Hazard ratio inset | Rarely shown; typically only a p-value | HR + 95% CI + log-rank p-value in every two-arm figure |
| Landmark rates | Not shown | OS: 12-mo, 24-mo · PFS: 3-mo, 6-mo, 12-mo · EFS: 3-mo, 6-mo, 12-mo |
| DFS single-arm handling | HR not applicable; often omitted with no note | Explicit note: "Single-arm; no control arm — HR analysis not applicable" |
| EFS composite event breakdown | Not shown | Event counts split by progression vs. death without progression |
| Death-without-progression marker | Not shown | Purple triangle annotation with label at patient time-point |
| Landmark vertical guides | Not shown | Dotted vertical lines at clinically meaningful time-points |
| Confidence interval method | `plain` (default) | `log-log` — standard for clinical reporting (Collett, 2015) |
| Risk table styling | Default survminer | Custom `theme_cleantable()` + italic title + bold arm labels |
| Reproducibility | No fixed seed | `set.seed(42)` in upstream data generation; bitwise-identical |
| Calibration traceability | Unspecified | Anchored to KEYNOTE-189, IMbrave150, OlympiAD, KEYNOTE-177, NAPOLI-1 |

---

## Why these four endpoints complement each other

The four curves answer distinct clinical questions from the same patient population, and regulatory dossiers routinely require all of them for a comprehensive efficacy profile.

| Question | OS | PFS | EFS | DFS |
|---|---|---|---|---|
| Does treatment extend life? | ✓ | — | — | — |
| Does treatment delay radiological progression? | — | ✓ | ✓ | — |
| Does treatment reduce all clinically meaningful events (progression, death, toxicity-driven discontinuation)? | — | — | ✓ | — |
| In patients who achieved CR, how long do they remain disease-free? | — | — | — | ✓ |
| Appropriate for two-arm comparison with HR? | ✓ | ✓ | ✓ | — (single-arm) |
| Landmark rates reported | 12-mo, 24-mo | 3-mo, 6-mo, 12-mo | 3-mo, 6-mo, 12-mo | 6-mo, 12-mo, 24-mo |

---

## Visual anatomy

Each two-arm figure (OS, PFS, EFS) shares the following structure:

```
  Probability
  1.00 ──────────────────────────────────────────────────────
       │╲  shaded 95% CI band
  0.75 │ ╲──── Treatment arm (blue)
       │   ╲
  0.50 │    ╲──────── Control arm (red)         ┌──────────────────┐
       │          ╲                              │ HR (95% CI): x.xx│
  0.25 │           ╲                            │ Log-rank: P=x.xxx│
       │            ╲           ┆       ┆       │ Median ...        │
  0.00 │─────────────────────────────────────── │ Landmark rates    │
        0      6     12    18    24    30    36  └──────────────────┘
                 Time from Randomization (Months)
  ─────────────────────────────────────────────────────────────────
  Number at risk
  Control     18    14     9    ...
  Treatment   62    55    48    ...
```

| Element | Description |
|---|---|
| Blue step line | Treatment arm KM estimate |
| Red step line | Control arm KM estimate |
| Shaded band | 95% confidence interval (log-log method, α = 0.10 fill) |
| `\|` censor marks | Patients censored at that time-point (censor.shape = `\|`, size = 5) |
| Dotted vertical lines | Landmark time-points (6, 12, 24 mo for OS; 3, 6, 12 mo for PFS/EFS) |
| Bottom-right inset box | HR + 95% CI · log-rank p-value · median OS/PFS/EFS · landmark rates |
| Risk table | Number at risk per arm, colored by stratum, italic title "Number at risk" |
| Footer caption | Abbreviation key and dataset identifier |

The DFS figure (single-arm) has a simplified layout: one teal step curve, per-patient markers at the bottom (× for recurrence, `\|` for censored), tumor type labels per patient, and a stats box noting that HR analysis is not applicable.

---

## Output figures

| File | Endpoint | Arms | Key stats shown |
|---|---|---|---|
| `km_overall_survival.png` | OS | Treatment vs Control | Log-rank p, median OS per arm |
| `os_curve.png` | OS (full reporting standard) | Treatment vs Control | HR + 95% CI, log-rank p, median OS, 12-mo and 24-mo rates |
| `pfs_curve.png` | PFS | Treatment vs Control | HR + 95% CI, log-rank p, median PFS, 3-mo, 6-mo, 12-mo rates |
| `efs_curve.png` | EFS | Treatment vs Control | HR + 95% CI, log-rank p, median EFS, 3/6/12-mo rates, event breakdown |
| `dfs_curve.png` | DFS | Single-arm (CR patients) | Median DFS, 6/12/24-mo rates, per-patient tumor type annotation |

- Figure dimensions: 13 × 10.5 in @ 200 DPI
- Output directory: `./Outputs/`

---

## Endpoint definitions

| Endpoint | PARAMCD | Event | Censoring |
|---|---|---|---|
| Overall Survival (OS) | `OS` | Death from any cause | Last known alive date |
| Progression-Free Survival (PFS) | `PFS` | Radiological progression (RECIST 1.1) or death | Last tumor assessment without progression |
| Event-Free Survival (EFS) | `EFS` | Progression, death from any cause, or treatment discontinuation due to toxicity | Last event-free assessment |
| Disease-Free Survival (DFS) | `DFS` | Disease recurrence, progression, or death from any cause | Last disease-free assessment |

EFS event type is encoded in `COMPTYPE`: `PROGRESSION` or `DEATH_WITHOUT_PROGRESSION`. The EFS figure marks death-without-progression patients with a purple triangle (▲) at their event time.

---

## Statistical methods

### Kaplan–Meier estimator

Survival functions are estimated using `survfit()` from the `survival` package with `conf.type = "log-log"`. The log-log (complementary log-log) transformation is the regulatory standard for confidence interval construction (ICH E9, FDA guidance on time-to-event endpoints) and is more reliable than the default `plain` method when survival probabilities are near 0 or 1.

```r
km_fit <- survfit(
  Surv(AVALM, EVENT) ~ ARM,
  data      = km_df,
  conf.type = "log-log"
)
```

### Hazard ratio and confidence interval

The hazard ratio is estimated from a Cox proportional hazards model (`coxph`). The control arm is set as the reference level via `relevel()` so that HR < 1 favors treatment throughout.

```r
get_hr <- function(df_param, ..., ref = "CONTROL") {
  df[[arm_col]] <- relevel(factor(df[[arm_col]]), ref = ref)
  fit  <- coxph(Surv(AVALM, EVENT) ~ ARM, data = df)
  smry <- summary(fit)
  list(
    hr   = smry$conf.int[1, "exp(coef)"],
    lo   = smry$conf.int[1, "lower .95"],
    hi   = smry$conf.int[1, "upper .95"],
    pval = smry$sctest["pvalue"]
  )
}
```

### Log-rank test

The log-rank (score) test p-value is extracted from `summary(coxph_fit)$sctest["pvalue"]`. For figures, p-values below 0.001 are displayed as `P < 0.001`; otherwise they are formatted to three decimal places.

### Landmark survival rates

Landmark rates are computed with `summary(survfit, times = t, extend = TRUE)$surv` and multiplied by 100. `extend = TRUE` is required to return a value when no event occurs exactly at the requested time-point; without it, the call returns `NULL` for time-points beyond the last observed event in a small arm.

---

## Key implementation details

### Event derivation from CNSR

ADTTE follows the ADaM convention where `CNSR = 1` means censored and `CNSR = 0` means an event. The binary event indicator expected by `Surv()` is therefore:

```r
adtte$EVENT <- 1L - as.integer(adtte$CNSR)
```

### Time axis units

All time axes are in months (`AVALM`), which is pre-derived in ADTTE using the standard ADaM 30.4375 days/month conversion factor (365.25 / 12). No in-script conversion is required.

### Arm ordering

`ARM` values after `trimws()` are `"TREATMENT"` and `"CONTROL"`. `ggsurvplot` assigns palette colors in alphabetical order of factor levels, so `CONTROL` maps to the first palette color (red, `#D44C36`) and `TREATMENT` to the second (blue, `#1A73E8`). `scale_color_manual` and `scale_fill_manual` are added explicitly after the fact to ensure deterministic color assignment regardless of factor level order.

### DFS single-arm population

The DFS population is restricted to subjects with `PARAMCD = "DFS"` joined to ADSL for `PRIORSURG` and `PHASE`. The cohort consists of complete responders only, as encoded in the upstream data generation; no additional filtering on `BESTRSPC` is applied in this script because the DFS endpoint is already scoped to CR patients in ADTTE. A note is included in the stats inset confirming single-arm status and the inapplicability of HR analysis.

---

## Dataset — ONCVIZ-001 ADaM v1

```
Study    ONCVIZ-001 · Vizatinib 300 mg QD vs Control
Design   Phase I/II Open-Label Randomized Basket Trial
N        400  (TRT = 263 · CTL = 137 · 2:1 ratio)
Survival N    80   (TRT = 62  · CTL = 18  — Phase I/II subset)
Records  131,690 across 13 ADaM domains
Seed     42 · fully reproducible
Cutoff   05 Mar 2026
```

### Domains used by this script

| Domain | Description | Key variables used |
|--------|-------------|---|
| ADTTE | Time-to-event analysis dataset | `USUBJID`, `PARAMCD`, `AVALM`, `CNSR`, `EVENT`, `ARM`, `TUMORTYPE`, `COMPTYPE` |
| ADSL | Subject-level dataset | `USUBJID`, `PHASE`, `DOSELEVEL`, `BESTRSPC`, `PRIORSURG` |

```r
adtte <- read.csv(file.path(DATA_DIR, "ADTTE.csv"), stringsAsFactors = FALSE)
adsl  <- read.csv(file.path(DATA_DIR, "ADSL.csv"),  stringsAsFactors = FALSE)
adtte$EVENT <- 1L - as.integer(adtte$CNSR)
```

### Population per figure

| Figure | PARAMCD | n (TRT) | n (CTL) | Events (TRT) | Events (CTL) |
|--------|---------|--------:|--------:|-------------:|-------------:|
| KM / OS curve | `OS` | 62 | 18 | 31 | 13 |
| PFS curve | `PFS` | 62 | 18 | Derived | Derived |
| EFS curve | `EFS` | 62 | 18 | Derived | Derived |
| DFS curve | `DFS` | Single-arm CR patients | — | Derived | — |

---

## Key parameters

| Parameter | Value | Description |
|---|---|---|
| `conf.type` | `"log-log"` | CI method for `survfit` (complementary log-log, regulatory standard) |
| `conf.int.alpha` | `0.10` | Fill opacity for CI bands |
| `censor.shape` | `"\|"` | Tick mark for censored observations |
| `break.x.by` (OS/KM) | `6` | X-axis tick interval in months |
| `break.x.by` (PFS/EFS) | `3` | X-axis tick interval in months |
| `xlim` (OS/KM) | `c(0, 44)` | X-axis range in months |
| `xlim` (PFS/EFS) | `c(0, 43)` | X-axis range in months |
| OS landmark times | 12, 24 mo | Dotted verticals and inset rates |
| PFS/EFS landmark times | 3, 6, 12 mo | Dotted verticals and inset rates |
| DFS landmark times | 6, 12, 24 mo | Inset rates only |
| HR reference arm | `"CONTROL"` | Set via `relevel(factor(...), ref = "CONTROL")` |
| Color: Treatment | `#1A73E8` | Blue — consistent with all ONCVIZ-001 figures |
| Color: Control | `#D44C36` | Red — consistent with all ONCVIZ-001 figures |
| Color: DFS single-arm | `#20776A` | Teal |
| Color: DWP marker (EFS) | `#6929C4` | Purple triangle |

---

## When to use

**Appropriate:**
- Reporting time-to-event primary and secondary endpoints in Phase I/II or Phase II/III oncology trial publications and regulatory dossiers
- Presenting OS and PFS as co-primary endpoints with a shared visual theme and consistent stats inset format
- EFS when the protocol defines a composite event including toxicity-driven discontinuation alongside progression and death
- DFS for complete responder cohorts where a control arm is absent and HR analysis is not applicable
- Any context requiring both the survival curve and the risk table in a single self-contained figure for journal submission or slide presentation

**Limitations:**
- The Kaplan–Meier estimator assumes non-informative censoring; if censoring is related to patient prognosis (e.g., sicker patients discontinue earlier), survival will be overestimated (Collett, 2015)
- The log-rank test and Cox model assume proportional hazards; immuno-oncology agents frequently produce delayed separation or crossing curves that violate this assumption — in such cases the weighted log-rank test or restricted mean survival time (RMST) is a more appropriate primary analysis (Uno et al., 2014; Royston & Parmar, 2013)
- Median survival cannot be estimated if fewer than 50% of patients have experienced the event by the data cutoff; the table entry will show `NR` (not reached) — the landmark rate at a pre-specified time-point is a better communicable summary in this setting
- With small arm sizes (n = 18 for Control), confidence intervals are wide and the median may be unreliable; this is a known property of Phase I/II designs and is noted in the figure subtitle
- The proportional hazards assumption is not formally tested in this script; users should run `cox.zph()` on the coxph fit before reporting the HR as a valid summary statistic
- Competing risks (e.g., death before progression for PFS) are not handled here; a cumulative incidence function (CIF) via `cmprsk` or `survminer::ggcompetingrisks()` is more appropriate when competing events are non-negligible

---

## Files

| File | Description | Dimensions |
|---|---|---|
| `survival_curves.R` | Main R script — all four curves | |
| `Data/V1/ADTTE.csv` | Time-to-event analysis dataset | 400+ rows |
| `Data/V1/ADSL.csv` | Subject-level dataset | 400 rows |
| `Out/km_overall_survival.png` | KM estimate of OS | 13 × 10.5 in · 200 DPI |
| `Out/os_curve.png` | OS — full clinical reporting standard | 13 × 10.5 in · 200 DPI |
| `Out/pfs_curve.png` | Progression-Free Survival | 13 × 10.5 in · 200 DPI |
| `Out/efs_curve.png` | Event-Free Survival | 13 × 10.5 in · 200 DPI |
| `Out/dfs_curve.png` | Disease-Free Survival (single-arm, CR) | 13 × 10.5 in · 200 DPI |

---

## Requirements

```r
library(survival)    # >= 3.5  — survfit, coxph, Surv
library(survminer)   # >= 0.4  — ggsurvplot, theme_cleantable
library(ggplot2)     # >= 3.4
library(dplyr)       # >= 1.1
library(patchwork)   # >= 1.2
library(grid)        # textGrob, gpar, unit
library(gridExtra)   # grid.arrange
library(scales)      # label helpers
```

No proprietary packages required. All dependencies are available on CRAN.

---

## References

Collett D. *Modelling Survival Data in Medical Research.* 3rd ed. CRC Press; 2015.

Eisenhauer EA, et al. New response evaluation criteria in solid tumours: RECIST 1.1. *Eur J Cancer.* 2009;45(2):228–247.

Finn RS, et al. Atezolizumab plus bevacizumab in unresectable hepatocellular carcinoma. *N Engl J Med.* 2020;382(20):1894–1905.

Gandhi L, et al. Pembrolizumab plus chemotherapy in metastatic non-small-cell lung cancer. *N Engl J Med.* 2018;378(22):2078–2092.

Kassambara A, Kosinski M, Biecek P. survminer: Drawing Survival Curves using 'ggplot2'. R package version 0.4.9. https://CRAN.R-project.org/package=survminer

Royston P, Parmar MKB. Restricted mean survival time: an alternative to the hazard ratio for the design and analysis of randomized trials with a time-to-event outcome. *BMC Med Res Methodol.* 2013;13:152.

Therneau TM. A Package for Survival Analysis in R. R package version 3.5-8. https://CRAN.R-project.org/package=survival

Uno H, et al. Moving beyond the hazard ratio in quantifying the between-group difference in survival analysis. *J Clin Oncol.* 2014;32(22):2380–2385.

Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
