# Trough Level Plot

A steady-state pharmacokinetic visualization displaying pre-dose (trough) Vizatinib concentrations at the RP2D (300 mg QD) across two visits, from a fully synthetic phase I/II oncology trial dataset. Because trough concentrations in this compound are almost entirely below the limit of quantification, the plot design foregrounds that finding honestly rather than forcing a conventional boxplot onto degenerate data.

**Dataset:** ONCVIZ-001 · N = 47 (300 mg RP2D, PK-evaluable) | **Cutoff:** 05 March 2026 | **Language:** R · ggplot2 · patchwork | **License:** CC BY 4.0

---

## The gap this fills

Most trough-level plot templates assume a boxplot is always appropriate. When the BLQ rate is high — as it genuinely is here, given this compound's short half-life relative to a 24h QD interval — a boxplot of mostly-zero values is actively misleading. This script implements the standard alternative: report the BLQ rate explicitly, and only summarize the quantifiable fraction.

| | Prior art | This work |
|---|---|---|
| BLQ handling | Boxplot regardless of BLQ rate | BLQ rate reported explicitly; summary stats suppressed above 50% BLQ |
| Label placement | Often clipped or overflowing bar edges | % BLQ centered inside the bar; % quantifiable placed with clear headroom |
| Quantifiable values | Buried inside a near-empty boxplot | Shown separately, log scale, individual points |
| Visits compared | Single timepoint | Cycle 1 Day 1 vs Cycle 3 Day 1 (steady state) |

---

## Visual anatomy

```
100% |███████████████████████|  ← BLQ (gray), % labeled centered in bar
     |███████████████████████|
     |___________Quantifiable_|  ← thin blue segment, % labeled above
  0% |________________________|
       C1D1              C3D1
```

| Element | Description |
|---|---|
| Left panel | Stacked bar — % BLQ (gray) vs % quantifiable (navy), per visit |
| Right panel | Individual quantifiable trough values only, log scale, jittered points |
| "X% BLQ" label | Centered inside the gray segment (never overflows the bar) |
| "X% quantifiable" label | Placed above the bar with clear headroom |

---

## Why no boxplot

Per standard PK reporting convention, median/IQR summary statistics are not presented for a visit where the BLQ rate exceeds 50%. Both visits here are ~98% BLQ:

```r
# trough_level.R
pct_blq <- 100 * mean(vals <= 0)   # ~97.9% at both C1D1 and C3D1
# → boxplot suppressed; BLQ-rate panel + quantifiable-only strip shown instead
```

**Clinical interpretation:** minimal pre-dose exposure at steady state indicates rapid clearance relative to the 24h QD dosing interval and negligible accumulation — itself a reportable PK finding, not a plotting failure.

---

## Dataset — ONCVIZ-001 ADaM v1

| Domain | Description | Key variables |
|--------|-------------|---|
| ADPK | Pharmacokinetic concentrations | `PARAMCD="TROUGH"`, `DOSE`, `AVISIT`, `AVAL`, `ANL01FL` |

```r
adpk <- read.csv(file.path(DATA_DIR, "ADPK.csv"), stringsAsFactors = FALSE)
tr <- adpk |> filter(PARAMCD == "TROUGH", DOSE == 300, ANL01FL == "Y")
```

---

## Output files

| File | Description | Dimensions |
|------|-------------|------------|
| `trough_level.R` | Main R script | |
| `Out/trough_level.png` | BLQ-rate panel + quantifiable-troughs panel | 11×6 in · 300 DPI |

---

## When to use

**Appropriate:**
- Steady-state accumulation assessment at the recommended phase 2 dose
- Reporting BLQ-dominated PK data honestly rather than suppressing the finding
- Supporting a rapid-clearance / low-accumulation narrative with primary data

**Limitations:**
- Not informative when the BLQ rate is low — use a standard boxplot + jitter instead
- Only two visits compared here; extend the `visit_order` vector for additional cycles
- Does not distinguish assay failure from true sub-LLOQ exposure — check bioanalytical run records if the BLQ rate is unexpectedly high

---

## Requirements

```
R >= 4.1
```

```r
install.packages(c("dplyr", "ggplot2", "patchwork", "tidyr"))
```

---

## References

Beal SL. Ways to fit a PK model with some data below the quantification limit. *J Pharmacokinet Pharmacodyn.* 2001;28(5):481–504.

U.S. Food and Drug Administration. *Guidance for Industry: Bioanalytical Method Validation.* CDER/CVM; 2018.

Gabrielsson J, Weiner D. *Pharmacokinetic and Pharmacodynamic Data Analysis: Concepts and Applications.* 5th ed. Swedish Pharmaceutical Press; 2016.

