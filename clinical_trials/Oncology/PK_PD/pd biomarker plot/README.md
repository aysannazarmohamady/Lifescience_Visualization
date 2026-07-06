# PD Biomarker Plot

A longitudinal pharmacodynamic visualization displaying circulating tumor DNA (ctDNA) percent change from baseline over study cycles, stratified by best overall response (BOR), from a fully synthetic phase I/II oncology trial dataset.

**Dataset:** ONCVIZ-001 · N = 62 (treatment arm, PD-evaluable) | **Cutoff:** 05 March 2026 | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

PD biomarker plots are often placed next to PK concentration-time plots with no acknowledgment that the two are usually sampled on completely different time scales — hours (single-dose PK) versus weeks (longitudinal PD). This script makes that distinction explicit rather than implying a shared axis by proximity alone.

| | Prior art | This work |
|---|---|---|
| Baseline normalization | Raw values | % change from baseline (PCHG), consistent with PD summary convention |
| Error statistic | SD or unstated | Mean ± SEM, explicitly labeled |
| Time axis | Cycle number only | Cycle number **and** a secondary top axis in nominal study day |
| PK/PD time-scale mismatch | Implied comparable, often confusing readers | Explicit caption note: PK is intensive single-dose (hours), PD is longitudinal (weeks) |
| Color convention | Ad hoc | Same `RESP_COLORS` palette used across every other plot in this repo |

---

## Visual anatomy

```
Day:    0    21    42    63    84       126        168   ← secondary top axis
Cycle:  0     1     2     3     4        6           8   ← primary bottom axis
+100 |                                            ___●   PD (progressive disease)
   0 |------●-----------------------------------------   ← baseline
 -100 |___●________________________________________      CR (complete response)
```

| Element | Description |
|---|---|
| Bottom x-axis | Study cycle (`AVISITN`) |
| Top x-axis | Nominal study day = cycle × 21 (cycle length confirmed exact from `ADTN`) |
| Line + ribbon | Mean ± SEM ctDNA % change from baseline, one line per BOR group |
| Color | `RESP_COLORS` — navy (CR), light blue (PR), amber (SD), red (PD) |
| Dashed line at 0% | No change from baseline |

---

## Why two time axes, not one shared with the PK plot

PK sampling in this dataset is **intensive single-dose** (Cycle 1 Day 1, hours 0–24). PD sampling is **longitudinal** (every cycle, weeks to months). Forcing both onto one shared axis would either compress the PK curve to a sliver or stretch the PD curve implausibly. Instead:

```r
scale_x_continuous(breaks = visits,
  sec.axis = sec_axis(~ . * CYCLE_DAYS, name = "Nominal study day post first dose (cycle × 21)"))
```

The top axis lets a reader convert cycle → day for cross-referencing without pretending the two plots share a single time base.

---

## Dataset — ONCVIZ-001 ADaM v1

| Domain | Description | Key variables |
|--------|-------------|---|
| ADBM | Biomarker assessments | `PARAMCD="CTDNA"`, `AVISITN`, `PCHG`, `BESTRSPC`, `ARM` |

```r
adbm <- read.csv(file.path(DATA_DIR, "ADBM.csv"), stringsAsFactors = FALSE)
ct <- adbm |> filter(PARAMCD == "CTDNA", ARM == "TREATMENT", AVISITN >= 0)
```

---

## Output files

| File | Description | Dimensions |
|------|-------------|------------|
| `pd_biomarker.R` | Main R script | |
| `Out/pd_biomarker.png` | ctDNA dynamics by BOR, dual time axis | 9×7.2 in · 300 DPI |

---

## When to use

**Appropriate:**
- Demonstrating ctDNA clearance kinetics differ meaningfully by response category
- Supporting an early on-treatment biomarker as a candidate response predictor
- Any PD biomarker with a validated baseline (%CHG interpretable)

**Limitations:**
- N per group shrinks at later cycles due to progression-driven discontinuation; per-timepoint N is not shown on the curve itself (see source data)
- Only 4 BOR categories shown; NE patients excluded
- Not a substitute for a formal landmark or joint PK/PD model — this is a descriptive summary plot

---

## Requirements

```
R >= 4.1
```

```r
install.packages(c("dplyr", "ggplot2"))
```

---

## References

Bratman SV, et al. Personalized circulating tumor DNA analysis as a predictive biomarker in solid tumor patients treated with pembrolizumab. *Nat Cancer.* 2020;1(9):873–881.

Gabrielsson J, Weiner D. *Pharmacokinetic and Pharmacodynamic Data Analysis: Concepts and Applications.* 5th ed. Swedish Pharmaceutical Press; 2016.

