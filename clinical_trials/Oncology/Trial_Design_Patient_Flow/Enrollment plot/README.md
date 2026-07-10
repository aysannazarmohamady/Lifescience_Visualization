# Enrollment Over Time Plot

A cumulative enrollment curve overlaid with monthly accrual bars for the ONCVIZ-001 synthetic randomized basket trial, from first patient in through the last randomized patient.

**Dataset:** ONCVIZ-001 · n = 80 randomized · ADRAND / ADSL | **Endpoint:** Accrual over calendar time | **Language:** R · ggplot2 · lubridate | **License:** CC BY 4.0

---

## The gap this fills

Enrollment curves are simple in principle but are often published either as a bare cumulative line (losing the sense of accrual *rate*) or as monthly bars alone (losing the cumulative total). Combining both on dual axes, with milestone annotations, is standard in DSMB and steering-committee decks but rarely shown in public code examples.

| | Prior art | This work |
|---|---|---|
| Display | Cumulative line only, or bars only | Cumulative step curve + monthly accrual bars, dual y-axis |
| Milestones | Unlabeled | "First patient in" and "Target enrollment reached" annotated |
| Site breakdown | Sometimes stacked, often unreadable at this scale | Omitted by default (15 sites, 2–9 pts each — too sparse to stack meaningfully); left as an extension point |

---

## Visual anatomy

```
 Cumulative                                    Monthly
 80 ┤                                    ╱▔    │ ██
    │                                ╱▔▔       │ ██  ██
    │                          ╱▔▔              │ ██  ██  ██
    │                    ╱▔▔                    │
  0 ┤________╱▔▔___________________________     └────────────
    2022-07        2023-01        2023-07
```

| Element | Description |
|---|---|
| Step line + shaded area | Cumulative randomized patients (left axis) |
| Red bars | Patients enrolled per calendar month (right axis) |
| Annotations | First patient in; target enrollment reached |

---

## Dataset variables used

| Variable | Source | Description |
|---|---|---|
| `RANDFL` | ADRAND | Randomization flag |
| `TRTSDT` | ADSL | First dose date, used as the enrollment-date proxy |

---

## Statistical method

None — purely descriptive. Cumulative count is a running total sorted by date; monthly accrual is a `floor_date(TRTSDT, "month")` count. The two series are placed on a shared plot via a secondary axis scaled by `max(cumulative) / max(monthly)`.

---

## Key parameters

| Parameter | Value |
|---|---|
| Enrollment-date proxy | `TRTSDT` (first dose date) |
| Aggregation | Monthly bins |
| Population | Randomized, ITT (n = 80) |

---

## Limitations

- **No planned-accrual curve is available in the source data**, so this figure cannot show actual-vs-planned enrollment — a standard element in real trial reporting. Add a `PLANNED_N` reference series if a target curve becomes available.
- **`TRTSDT` is a proxy for enrollment/randomization date**, not a dedicated `RANDDT` field in this ADaM version; the two are expected to be close but not always identical in a real trial (e.g., a run-in period between randomization and first dose).
- **No per-site breakdown** shown; with 15 sites averaging ~5 patients each, a stacked-by-site version would be visually noisy without aggregating sites into regions.

---

## Requirements

```r
library(ggplot2)    # >= 3.4
library(dplyr)      # >= 1.1
library(lubridate)  # >= 1.9
```

---

## Files

| File | Description |
|---|---|
| `enrollment_plot.R` | Self-contained script · reads `Data/V1/ADRAND.csv`, `Data/V1/ADSL.csv` |
| `Out/enrollment_over_time.png` | Output figure · 10 × 6 in · 300 DPI |

---

## References

ICH E6(R2). *Good Clinical Practice: Integrated Addendum.* Section on trial monitoring and recruitment reporting.

CTTI (Clinical Trials Transformation Initiative). Recommendations on enrollment and recruitment metrics reporting in clinical trials.
