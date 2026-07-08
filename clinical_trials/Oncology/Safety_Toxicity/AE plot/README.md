# AE Bar Chart — Adverse Events by Preferred Term and Treatment Arm

A patient-level safety summary showing incidence of adverse events (any grade and Grade ≥3) by MedDRA-style Preferred Term, split by treatment arm, for the ONCVIZ-001 synthetic basket trial.

**Dataset:** ONCVIZ-001 · Treatment (n = 62) vs. Control (n = 18) · ADAE / ADSL | **Endpoint:** Adverse event incidence (safety population) | **Language:** R · ggplot2 · patchwork | **License:** CC BY 4.0

---

## The gap this fills

AE tables in clinical study reports are exhaustive but hard to scan; AE bar charts condense them, but many implementations show raw counts instead of incidence, mix grades together, or omit the denominator needed to interpret a percentage.

| | Prior art | This work |
|---|---|---|
| Metric | Raw event counts | % of patients with ≥1 occurrence (incidence), denominator shown in legend |
| Grade split | Single panel, all grades pooled | Two panels: any-grade and Grade ≥3 (CTCAE), side by side |
| Sorting | Alphabetical | Descending incidence in Treatment arm |
| Threshold | None (long tail of rare terms) | ≥5% incidence in either arm |
| Arm comparison | Single-arm or stacked | Mirrored/grouped horizontal bars, shared category axis |
| Legend placement | Ad hoc | Collected once, above both panels, clear of bar area |

---

## Visual anatomy

```
 Fatigue      ████████████████ Treatment        │ ██████████ Grade≥3
 Nausea       ██████████████                     │ ████
 Diarrhoea    █████████████       ██ Control      │ ████████████
              ...                                 │
              100%    50%    0%  |  0%   10%   20%
```

| Element | Description |
|---|---|
| Left panel | Any-grade AE incidence (%), both arms, mirrored bars |
| Right panel | Grade ≥3 AE incidence (%), same PT order |
| Bar color | Blue = Treatment, Red = Control |
| Sort order | Descending incidence in Treatment arm (shared across both panels) |
| Threshold | Only PTs with ≥5% incidence in either arm are shown |

---

## Dataset variables used

| Variable | Source | Description |
|---|---|---|
| `USUBJID` | ADAE / ADSL | Subject identifier |
| `AEPT` | ADAE | Adverse event Preferred Term |
| `AESOC` | ADAE | System Organ Class |
| `AETOXGR` | ADAE | CTCAE toxicity grade (1–5) |
| `ARM` | ADAE / ADSL | "TREATMENT" · "CONTROL" |
| `CTCAE_V` | ADAE | CTCAE version used for grading |

---

## Statistical method

Incidence per Preferred Term and arm is computed at the **patient level** (a patient with multiple occurrences of the same AE is counted once):

```r
n_patients_with_AE / n_patients_in_arm * 100
```

Grade ≥3 incidence uses the same formula restricted to `AETOXGR >= 3`. No formal between-arm statistical test (e.g. Fisher's exact) is applied; the chart is descriptive, consistent with standard CSR safety summary tables.

---

## Key parameters

| Parameter | Value |
|---|---|
| Reporting threshold | ≥5% incidence in either arm |
| Grade cutoff | CTCAE Grade ≥3 for the right panel |
| Sort key | Treatment-arm any-grade incidence, descending |
| Population | Safety population (`SAFFL = "Y"`) |

---

## Limitations

- **No multiplicity adjustment.** With ~25+ Preferred Terms displayed, apparent arm differences should be interpreted descriptively, not as confirmatory signals.
- **Patient-level incidence hides event frequency.** A patient with 5 episodes of nausea counts the same as a patient with 1.
- **Small control arm (n = 18).** Percentages for Control are based on very few patients; single-patient events swing the bar substantially.

---

## Requirements

```r
library(ggplot2)    # >= 3.4
library(dplyr)      # >= 1.1
library(tidyr)       # >= 1.3
library(patchwork)  # >= 1.2
```

---

## Files

| File | Description |
|---|---|
| `ae_bar_plot.R` | Self-contained script · reads `Data/V1/ADAE.csv`, `Data/V1/ADSL.csv` |
| `Out/ae_bar_chart.png` | Output figure · 12 × 7 in · 300 DPI |

---

## References

U.S. Department of Health and Human Services, National Cancer Institute. *Common Terminology Criteria for Adverse Events (CTCAE) v5.0.* 2017.

FDA Center for Drug Evaluation and Research. *Reviewer Guidance: Evaluating the Risks of Drug-Induced Liver Injury and Other Safety Signals in Premarketing Clinical Trials* and related graphical safety-reporting guidance on clinical trial safety data presentation, 2009 draft framework for standard safety graphics.

International Council for Harmonisation. *ICH E3: Structure and Content of Clinical Study Reports.* Section 12 (Safety Evaluation).
