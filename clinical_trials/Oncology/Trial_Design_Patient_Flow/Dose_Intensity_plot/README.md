# Dose Intensity Plot — Relative Dose Intensity by Dose Level

A box-and-jitter plot of per-patient relative dose intensity (RDI), stratified by planned dose level, for the treatment arm of the ONCVIZ-001 synthetic dose-escalation/expansion trial.

**Dataset:** ONCVIZ-001 · Treatment arm, n = 62 across 4 dose levels · ADEX / ADSL | **Endpoint:** Relative Dose Intensity (RDI) | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

Dose intensity is a required exposure metric in oncology CSRs but is frequently reported only as a single summary statistic (mean/median RDI) in a table, losing the distribution and outliers that explain dose modifications.

| | Prior art | This work |
|---|---|---|
| Display | Single summary number per arm | Box plot + jittered patient-level points, per dose level |
| Reference lines | Often absent | 100% planned and 80% clinically-relevant-reduction lines, both labeled |
| Grouping | Pooled across dose levels | Stratified by planned dose level (100/200/300/400 mg) |
| N per group | Implicit | Shown directly in the x-axis label |

---

## Visual anatomy

```
 RDI(%)
 100 ┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄┄  ← Planned dose intensity
     │  ┬     ┬      ┬       ┬
     │ [█]   [█]    [███]   [█]
     │  ┴     ┴      ┴       ┴
  80 ················· ← reduction threshold
     100mg  200mg   300mg   400mg
     (n=5)  (n=5)  (n=47)   (n=5)
```

| Element | Description |
|---|---|
| Box | IQR of per-patient RDI within dose level |
| Points | Individual patients, jittered horizontally |
| Green dashed line | 100% (fully on-protocol dosing) |
| Red dotted line | 80% (clinically relevant reduction threshold) |

---

## Dataset variables used

| Variable | Source | Description |
|---|---|---|
| `USUBJID`, `DOSEINT` | ADEX | Subject, per-cycle dose intensity (actual/planned) |
| `DOSELEVEL` | ADSL | Planned dose level (mg) |
| `ARM` | ADEX | Restricted to "TREATMENT" |

---

## Statistical method

Per-patient RDI is the mean of per-cycle dose intensity across all administered cycles:

```
RDI_patient = mean(DOSEINT across cycles) × 100
```

This mirrors the standard oncology definition: RDI = (actual dose delivered / actual time) ÷ (planned dose / planned time) × 100.

---

## Key parameters

| Parameter | Value |
|---|---|
| Reference lines | 100% (planned), 80% (reduction threshold) |
| Grouping variable | `DOSELEVEL` (100/200/300/400 mg) |
| Aggregation | Mean per-cycle RDI per patient |

---

## Limitations

- **Unequal group sizes** (n = 5, 5, 47, 5) — the 300 mg (RP2D) group dominates visually and statistically; cross-dose-level comparisons should account for this imbalance.
- **No formal trend test** (e.g., Jonckheere-Terpstra) for RDI across dose levels is included.
- **Per-cycle DOSEINT averaging** treats early and late cycles equally; a patient who reduced dose only in later cycles will show the same mean RDI as one who reduced early, obscuring the timing of modifications (see the companion Treatment Exposure swimmer plot for timing detail).

---

## Requirements

```r
library(ggplot2)    # >= 3.4
library(dplyr)      # >= 1.1
```

---

## Files

| File | Description |
|---|---|
| `dose_intensity_plot.R` | Self-contained script · reads `Data/V1/ADEX.csv`, `Data/V1/ADSL.csv` |
| `Out/dose_intensity_plot.png` | Output figure · 8 × 6.5 in · 300 DPI |

---

## References

Hryniuk W, Bush H. The importance of dose intensity in chemotherapy of metastatic breast cancer. *J Clin Oncol.* 1984;2(11):1281–1288.

National Comprehensive Cancer Network / ASCO guidance on dose modification and relative dose intensity reporting in oncology trials.

Wildiers H, et al. Relative dose intensity of chemotherapy and its impact on outcomes in patients with early breast cancer or aggressive lymphoma. *Crit Rev Oncol Hematol.* 2016;101:26–37.
