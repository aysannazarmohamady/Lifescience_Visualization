# PK Concentration-Time Curve

A dose-level pharmacokinetic visualization displaying geometric mean (± 95% CI) Vizatinib plasma concentration versus nominal time post-dose at Cycle 1 Day 1, across all four dose-escalation cohorts, from a fully synthetic phase I/II oncology trial dataset.

**Dataset:** ONCVIZ-001 · N = 62 (PK-evaluable) | **Cutoff:** 05 March 2026 | **Language:** R · ggplot2 · patchwork | **License:** CC BY 4.0

---

## The gap this fills

Most PK figure examples either show a single dose group (hiding dose-proportionality) or omit the semi-log panel that reviewers expect for terminal-phase inspection. This script reproduces both panels side by side from one dataset, applies the BLQ-handling rule that most toy examples skip, and never fabricates a mean from a majority-BLQ timepoint.

| | Prior art | This work |
|---|---|---|
| Scale | Linear only, or log only | Linear + semi-log dual panel, same data |
| Summary statistic | Arithmetic mean | Geometric mean ± 95% CI (log-normal PK convention) |
| BLQ handling | Ignored or floored silently | >50% BLQ timepoints excluded from the mean profile, disclosed in caption |
| Dose groups | Single cohort | All 4 dose levels (100/200/300/400 mg) overlaid for dose-proportionality |
| Time basis | Unstated | Explicitly nominal sampling time, stated in caption |

---

## Visual anatomy

```
 2.0 |        ___
     |      _/   \___
 1.0 |    _/          \____        ← geometric mean, one line per dose
     |  _/                  \___
 0.0 |_/                         \_____
     0   1   2   4   6   8   12  24 (h) ← nominal time post-dose
```

| Element | Description |
|---|---|
| Left panel | Linear scale — visual peak/trough shape |
| Right panel | Semi-log scale — terminal elimination slope inspection |
| Point + line | Geometric mean concentration per dose, per nominal timepoint |
| Error bars | 95% CI (log-normal), computed only from quantifiable (>LLOQ) values |
| Missing timepoint | Dropped from the mean profile when >50% of concentrations at that timepoint are BLQ |
| Color | Dose cohort — light blue (100 mg) → dark navy (300 mg, RP2D) → red (400 mg) |

---

## BLQ handling rule

Per standard PK reporting convention, a mean concentration is **not** computed for any nominal timepoint where more than 50% of individual concentrations are below the limit of quantification (BLQ). This matters here: all three non-RP2D dose cohorts are >50% BLQ by 24h post-dose (drug fully eliminated), so the 24h point is dropped from those profiles rather than plotted as a near-zero geometric mean driven by one detectable outlier.

```r
# BLQ exclusion rule — pk_concentration_time.R
pct_blq <- 100 * mean(vals <= 0)
if (t > 0 && pct_blq > 50) {
  # timepoint excluded from the mean profile, logged in the caption
}
```

---

## Dataset — ONCVIZ-001 ADaM v1

```
Study    ONCVIZ-001 · Vizatinib dose-escalation + RP2D expansion
N        62 PK-evaluable patients (100/200/400 mg: n=5 each · 300 mg RP2D: n=47)
Domain   ADPK · PARAMCD = "CONC" · AVISIT = "CYCLE 1 DAY 1"
Seed     42 · fully reproducible
```

| Domain | Description | Key variables |
|--------|-------------|---|
| ADPK | Pharmacokinetic concentrations | `PARAMCD`, `NOMTPT`, `AVAL`, `DOSE`, `ANL01FL` |

```r
adpk <- read.csv(file.path(DATA_DIR, "ADPK.csv"), stringsAsFactors = FALSE)
conc <- adpk |> filter(PARAMCD == "CONC", AVISIT == "CYCLE 1 DAY 1", ANL01FL == "Y")
```

---

## Output files

| File | Description | Dimensions |
|------|-------------|------------|
| `pk_concentration_time.R` | Main R script | |
| `Out/pk_concentration_time.png` | Linear + semi-log dual-panel figure | 12×6 in · 300 DPI |

---

## When to use

**Appropriate:**
- First-in-human / dose-escalation PK reporting
- Assessing dose-proportionality across cohorts
- Regulatory submission figures requiring both linear and semi-log views

**Limitations:**
- Single-dose (Cycle 1) profile only — does not show steady-state accumulation (see Trough Level Plot)
- Nominal, not actual, sampling times are plotted — state this explicitly if actual times differ materially
- Not appropriate once BLQ rate exceeds 50% at most timepoints — switch to a BLQ-rate summary instead

---

## Requirements

```
R >= 4.1
```

```r
install.packages(c("dplyr", "ggplot2", "patchwork", "scales"))
```

---

## References

Gabrielsson J, Weiner D. *Pharmacokinetic and Pharmacodynamic Data Analysis: Concepts and Applications.* 5th ed. Swedish Pharmaceutical Press; 2016.

U.S. Food and Drug Administration. *Guidance for Industry: Population Pharmacokinetics.* CDER; 2022.

Beal SL. Ways to fit a PK model with some data below the quantification limit. *J Pharmacokinet Pharmacodyn.* 2001;28(5):481–504.

Mould DR, Upton RN. Basic concepts in population modeling, simulation, and model-based drug development. *CPT Pharmacometrics Syst Pharmacol.* 2012;1(9):e6.

