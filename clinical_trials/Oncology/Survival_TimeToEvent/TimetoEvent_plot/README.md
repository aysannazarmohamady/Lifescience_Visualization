# Oncology Survival Plots — Six-Panel Suite

Patient-level and cohort-level survival visualizations generated from RECIST 1.1 and time-to-event ADaM data across a fully synthetic phase I/II basket trial dataset (ONCVIZ-001), covering Time to Response, Time to Progression, Duration of Response, Landmark Overall Survival, Competing Risks CIF, and Restricted Mean Survival Time.

**Dataset:** ONCVIZ-001 · N = 80 | **Cutoff:** 05 March 2026 | **Language:** R · ggplot2 + survminer + cmprsk | **License:** CC BY 4.0

---

## The gap this fills

Most published oncology survival visualization suites rely on a single endpoint or present each figure from an independent dataset, making cross-figure internal consistency impossible to verify.

|  | Prior art | This work |
|---|---|---|
| Endpoint coverage | Single KM curve | 7 PARAMCD endpoints (OS, PFS, TTP, EFS, DOR, TTR, DFS) |
| Cross-domain consistency | Figures from different studies | All 6 plots from identical 80-patient cohort, same seed |
| Competing events | Not available | COMPTYPE encoded (PROGRESSION / DEATH\_WITHOUT\_PROGRESSION / CENSORED) |
| Landmark flags | Post-hoc | LM6MFL, LM12MFL, LM24MFL pre-computed in ADTTE |
| Competing risks method | 1 − KM (biased) | Aalen-Johansen CIF via `cmprsk::cuminc` |
| RMST | Not reported | Bootstrap CI at τ = 12, 18, 24, 30 months |
| Reproducibility | No fixed seed | `set.seed(42)` — bitwise-identical across R ≥ 4.1 |

---

## Visual anatomy

### Image 1 — Time to Response (TTR)

```
100% |___
     |   |____
  50%|........|_____________ ← 50% (median crossover)
     |                 ______
  0% |_________________________
     0    0.5   1.0   1.5 (Months)
```

| Element | Description |
|---|---|
| Left panel | Kaplan-Meier curve, responders only (TTR cohort) |
| Right panel | Individual TTR horizontal bar chart, sorted by response time |
| Blue bars | Treatment arm (Vizatinib) |
| Red bars | Control arm |
| Dashed vertical lines | Median TTR per arm (0.8 mo vs 0.5 mo) |
| Stats box (left) | n per arm · Median TTR |

---

### Image 2 — Time to Progression (TTP)

```
100% |\_
     |  \___
  50%|.......\___________  ← median crossover
     |              \____
  0% |_________________________
              Number at risk table
```

| Element | Description |
|---|---|
| KM curve | Both arms with 95% log-log CI shading |
| Number-at-risk table | Below x-axis, color-coded by arm |
| Efficacy box (top-right) | n · Median TTP · Log-rank p-value |
| Dotted horizontal line | 50% probability (median reference) |
| Dotted vertical lines | Median time per arm |

---

### Image 3 — Duration of Response (DOR)

```
100% |\_
     |  \____             ← CR curve
     |   \____\____       ← PR curve (longer median in this data)
  50%|...........\____    ← Overall dashed
  0% |_________________________
```

| Element | Description |
|---|---|
| Left panel | KM stratified by CR (n=15) · PR (n=7) · Overall (n=22) |
| Right panel | Individual DOR swimmer, sorted by duration |
| Dark green bars / line | Complete Response (CR) |
| Blue bars / line | Partial Response (PR) |
| Dashed overall curve | All responders combined |
| Circle marker | Event (loss of response) |
| Triangle marker | Censored |
| Stats box | Median DOR per response category |

---

### Image 4 — Landmark Analysis: Overall Survival

```
100%|\_
    |  \____
    |6m  12m  24m  ← landmark lines
 0% |________________________
     Bar chart      Diff + 95% CI
```

| Element | Description |
|---|---|
| Top panel | OS KM with vertical landmark lines at 6, 12, 24 months |
| Colored vertical lines | Gold = 6m · Purple = 12m · Teal = 24m |
| Bottom-left | Grouped bar chart: OS rate per arm at each landmark |
| Bottom-right | Absolute difference (Treatment − Control) with 95% bootstrap CI |
| Censor marks | `+` on KM curve at each censored observation |

---

### Image 5 — Competing Risks CIF

```
100%|
    |                Progression ────
 50%|            ___/
    |       ____/
    |  ____/        Death w/o Prog ──
  0%|_________________________
    0     12     24     36 (Months)
```

| Element | Description |
|---|---|
| Two panels | Left = Treatment (n=62) · Right = Control (n=18) |
| Blue step curve | Cumulative incidence of Progression (event of interest) |
| Red step curve | Cumulative incidence of Death without Progression (competing) |
| Shaded CI | 95% bootstrap confidence band (300 resamples, seed 42) |
| Info box | CIF values at 12 and 24 months for both event types |
| Method | Aalen-Johansen estimator via `cmprsk::cuminc` — not 1−KM |

> **Why not 1−KM?** In the presence of competing events, `1 − KM` overestimates the cumulative incidence of the event of interest. The correct estimator is the Aalen-Johansen CIF (Fine & Gray, 1999).

---

### Image 6 — Restricted Mean Survival Time (RMST)

```
100%|\_
    |  \____  ← shaded RMST area under treatment curve
    |         |
  0%|_________|____
              τ=24m
   RMST vs τ plot    RMST diff + CI
```

| Element | Description |
|---|---|
| Top panel | OS KM with shaded area under treatment curve up to τ = 24 months |
| Orange dashed line | Restriction time τ = 24 months |
| RMST box | RMST at τ=24m: Treatment = 15.1 mo · Control = 9.7 mo · Δ = +5.4 mo |
| Bottom-left | RMST vs τ line plot at τ = 12, 18, 24, 30 months |
| Diagonal dashed line | Reference: RMST = τ (perfect survival) |
| Bottom-right | RMST difference (Treatment − Control) + 95% bootstrap CI |
| Method | Left-Riemann integration of KM step function; CI via 300 bootstrap resamples |

---

## Dataset — ONCVIZ-001 ADaM v1

```
Study    ONCVIZ-001 · Vizatinib 300 mg QD vs Control
Design   Phase I/II Open-Label Randomized Basket Trial
N        80  (TRT = 62 · CTL = 18)
Records  26,723 across 13 ADaM domains
Seed     42 · fully reproducible
Cutoff   05 March 2026
```

### Tumor-stratified summary (treatment arm, OS endpoint)

| Histology | n (TRT) | n (CTL) | ORR TRT | ORR CTL | Median OS (TRT) |
|-----------|--------:|--------:|--------:|--------:|----------------:|
| NSCLC | 23 | 2 | 22% | 0% | 37.2 mo |
| HCC | 10 | 2 | 10% | 50% | 2.8 mo |
| CRC | 9 | 7 | 22% | 0% | 19.1 mo |
| BRCA | 8 | 5 | 25% | 20% | NR |
| PDAC | 12 | 2 | 8% | 0% | 8.8 mo |

### Endpoint summary (all arms)

| Endpoint | PARAMCD | TRT (n) | Events | CTL (n) | Events | Median TRT | Median CTL |
|----------|---------|--------:|-------:|--------:|-------:|----------:|----------:|
| Overall Survival | OS | 62 | 31 | 18 | 13 | NR | 5.3 mo |
| Progression-Free Survival | PFS | 62 | — | 18 | — | — | — |
| Time to Progression | TTP | 62 | 58 | 18 | 18 | 2.8 mo | 1.2 mo |
| Event-Free Survival | EFS | 62 | — | 18 | — | — | — |
| Duration of Response | DOR | 22 | 22 | 4 | 4 | 1.7 mo | — |
| Time to Response | TTR | 22 | 22 | 4 | 4 | 0.8 mo | 0.5 mo |
| Disease-Free Survival | DFS | — | — | — | — | — | — |

### ADaM domains

| Domain | Description | Rows | Key variables |
|--------|-------------|-----:|---|
| ADSL | Subject-level | 80 | `USUBJID`, `ARM`, `TUMORTYPE`, `LIVERMETS`, `BESTRSPC`, `PDL1GRP`, `TMBHIGH`, `MSISTS` |
| ADTTE | Time-to-event | 377 | `PARAMCD`, `AVAL`, `AVALM`, `CNSR`, `COMPTYPE`, `LM6MFL`, `LM12MFL`, `LM24MFL` |
| ADRS | Tumor response | 769 | `PARAMCD="OVRLRESP"`, `AVALC` (CR/PR/SD/PD/NE), `BICR_CONF` |
| ADTR | SLD measurements | 769 | `AVISIT`, `AVAL` (mm), `PCHG` (% change from baseline) |
| ADAE | Adverse events | 752 | `AESOC`, `AEDECOD`, `AESEV`, `AESER` |
| ADLB | Laboratory | 11,880 | `PARAMCD`, `AVAL`, `ANRIND` |
| ADEX | Exposure | 2,361 | `EXTRT`, `EXDOSE`, `EXDUR` |
| ADBM | Biomarkers | 3,321 | `PARAMCD`, `AVAL`, `ANRIND` |
| ADMUT | Mutations | 187 | `HUGO_SYMBOL`, `CHROMOSOME`, `VARIANT_CLASS` |
| ADPK | Pharmacokinetics | 1,984 | `PARAMCD`, `AVAL`, `NTIME` |
| ADPR | Patient-reported | 3,888 | `PARAMCD`, `AVAL`, `AVISIT` |
| ADSIG | Biomarker signatures | 263 | `SIG_NAME`, `SIG_SCORE` |
| ADRAND | Randomization | 92 | `RANDDT`, `STRATFL` |

---

## Output files

| File | Description | Dimensions |
|------|-------------|------------|
| `oncoviz_plots.R` | Main R script — all 6 plot variants | |
| `ADSL.csv` | Subject-level dataset | 80 rows |
| `ADTTE.csv` | Time-to-event dataset (7 endpoints) | 377 rows |
| `ADRS.csv` | Overall response assessments | 769 rows |
| `Image1_TTR.png` | Time to Response — KM + Individual bars | 18×8 in · 180 DPI |
| `Image2_TTP.png` | Time to Progression — KM + risk table | 13×11 in · 180 DPI |
| `Image3_DOR.png` | Duration of Response — KM by CR/PR + swimmer | 18×9 in · 180 DPI |
| `Image4_Landmark_OS.png` | Landmark OS — KM + bar + diff | 14×13 in · 180 DPI |
| `Image5_CIF.png` | Competing Risks — CIF both arms | 16×7 in · 180 DPI |
| `Image6_RMST.png` | RMST — KM shaded + RMST vs τ + diff | 15×13 in · 180 DPI |

---

## Evaluability criteria

A patient is **evaluable for response** (DOR / TTR endpoints) if they have:
- A confirmed best overall response of CR or PR in `ADRS` (`AVALC ∈ {CR, PR}`)
- A valid `AVALM` (months) in `ADTTE` with `CNSR` encoded

Patients without a post-baseline response are classified as Not Evaluable (NE) and excluded from ORR denominators in evaluable-population analyses.

```r
# Responder filter applied in Image 1 (TTR) and Image 3 (DOR)
ttr_df <- adtte %>% filter(PARAMCD == "TTR")   # n=26 (TRT=22, CTL=4)
dor_df <- adtte %>% filter(PARAMCD == "DOR",
                            ARM == "TREATMENT") %>%
          left_join(adsl %>% select(USUBJID, BESTRSPC), by = "USUBJID") %>%
          filter(BESTRSPC %in% c("CR", "PR"))   # n=22 (CR=15, PR=7)
```

---

## Statistical methods

| Plot | Method | Function | Notes |
|------|---------|----------|-------|
| KM curves | Kaplan-Meier | `survfit()` + `ggsurvplot()` | `conf.type = "log-log"` |
| HR | Cox proportional hazards | `coxph()` | CONTROL as reference |
| Log-rank test | Log-rank | `survdiff()` via `survminer` | p-value reported |
| CIF | Aalen-Johansen | `cmprsk::cuminc()` | Competing event: death w/o progression |
| RMST | Trapezoid integration of KM | Custom `calc_rmst()` | Left-Riemann step function |
| Bootstrap CI | Percentile bootstrap | `boot::boot()` | R = 300, seed = 42 |
| Landmark rates | Point estimate from KM | `summary(km, times = t)` | `extend = TRUE` for late landmarks |

---

## When to use each plot

| Plot | Use when | Do not use when |
|------|----------|-----------------|
| TTR (Image 1) | Reporting speed of response onset | No responders in cohort |
| TTP (Image 2) | Primary efficacy endpoint is progression control | Crossover design confounds TTP |
| DOR (Image 3) | Characterizing depth/durability of response | Most patients are non-responders (n < 10 responders) |
| Landmark OS (Image 4) | Immuno-oncology trials with delayed separation | Short follow-up < landmark time |
| CIF (Image 5) | Competing events exist (e.g. death before progression) | No competing events encoded |
| RMST (Image 6) | Non-proportional hazards or regulatory requirement | Short follow-up limits τ choice |

---

## Requirements

```
R >= 4.1
```

```r
install.packages(c(
  "survival",    # survfit, coxph, Surv
  "survminer",   # ggsurvplot, theme_cleantable
  "ggplot2",     # all plot layers
  "dplyr",       # data wrangling
  "patchwork",   # plot composition
  "gridExtra",   # grid.arrange
  "grid",        # textGrob, gpar
  "scales",      # percent_format
  "cmprsk",      # cuminc (Aalen-Johansen CIF)
  "boot"         # bootstrap CI for RMST and landmark differences
))
```

---

## Corrections applied vs. original figures

Two plots required correction before publication based on cross-figure consistency review:

| Plot | Original issue | Correction applied |
|------|---------------|-------------------|
| Image 3 — DOR | PR median (1.1 mo) ≥ CR median (1.7 mo) — biologically implausible | Re-estimated from actual ADTTE data: CR = 1.7 mo, PR = 4.3 mo |
| Image 5 — CIF | All CIF values = 0.0% at 12m and 24m — internally inconsistent with TTP curve showing 58/62 progression events | Replaced with Aalen-Johansen CIF: TRT Progression at 12m = 53.2%, at 24m = 67.7% |

---

## References

Eisenhauer EA, et al. New response evaluation criteria in solid tumours: RECIST version 1.1. *Eur J Cancer.* 2009;45(2):228–247.

Fine JP, Gray RJ. A proportional hazards model for the subdistribution of a competing risk. *J Am Stat Assoc.* 1999;94(446):496–509.

Gandhi L, et al. Pembrolizumab plus chemotherapy in metastatic non-small-cell lung cancer. *N Engl J Med.* 2018;378(22):2078–2092.

Finn RS, et al. Atezolizumab plus bevacizumab in unresectable hepatocellular carcinoma. *N Engl J Med.* 2020;382(20):1894–1905.

Robson M, et al. Olaparib for metastatic breast cancer in patients with a germline BRCA mutation. *N Engl J Med.* 2017;377(6):523–533.

Royston P, Parmar MKB. Restricted mean survival time: an alternative to the hazard ratio for the design and analysis of randomized trials with a time-to-event outcome. *BMC Med Res Methodol.* 2013;13:152.

Uno H, et al. On the utility of summary measures of the restricted mean survival time for the design and analysis of clinical trials. *J Clin Oncol.* 2014;32(16):1764.

Wolchok JD, et al. Overall survival with combined nivolumab and ipilimumab in advanced melanoma. *N Engl J Med.* 2017;377(14):1345–1356.

Zhang Z. Survival analysis in the presence of competing risks. *Ann Transl Med.* 2017;5(3):47.
