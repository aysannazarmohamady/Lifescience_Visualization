# Sum of Longest Diameters (SLD) Over Time — Spaghetti Plot

A patient-level spaghetti plot of percent change in target-lesion sum of longest diameters (SLD) from baseline over time, per RECIST 1.1, colored by best overall response, for all evaluable patients in ONCVIZ-001.

**Dataset:** ONCVIZ-001 · n = 80 evaluable patients · ADTR | **Endpoint:** Tumor burden over time (RECIST 1.1) | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

Longitudinal tumor-burden plots often smooth or interpolate between scan visits, implying continuous measurement between assessments that were never taken, and frequently omit the RECIST 1.1 response thresholds needed to interpret the trajectory.

| | Prior art | This work |
|---|---|---|
| Interpolation | Smoothed curves between visits | Straight-line connection between actual scan timepoints only, markers at each visit |
| Response thresholds | Often absent | +20% (PD) and −30% (PR) reference lines per RECIST 1.1 |
| Color coding | Single color or by arm | By best overall response category (CR/PR/SD/PD) |
| Legend | Ad hoc placement | Collected above the plot, clear of the line fan-out near baseline |

---

## Visual anatomy

```
 %Chg
 +20 ----------------------------  <- PD threshold
   0 -o---o---o---o---o---o--
 -30 ----------------------------  <- PR threshold
      0   6   12   18   24  weeks
```

| Element | Description |
|---|---|
| Line + point per patient | Actual scan timepoints, unsmoothed |
| Line color | Best overall response (CR = dark blue, PR = light blue, SD = orange, PD = red) |
| Dashed red line (+20%) | RECIST 1.1 progressive disease threshold |
| Dashed blue line (-30%) | RECIST 1.1 partial response threshold |

---

## Dataset variables used

| Variable | Source | Description |
|---|---|---|
| `USUBJID`, `ADTN`, `PCHG` | ADTR (`PARAMCD = "SUMDIAM"`) | Subject, study day, % change in SLD from baseline |
| `BESTRSPC` | ADSL | Best overall response category, used for line color |

---

## Statistical method

None — descriptive longitudinal display. `weeks = ADTN / 7`. `PCHG` is taken directly from the ADaM tumor-response dataset (already computed as % change from baseline SLD).

---

## Key parameters

| Parameter | Value |
|---|---|
| PD threshold | +20% from baseline (or from nadir, per full RECIST 1.1 definition) |
| PR threshold | -30% from baseline |
| Time unit | Weeks since baseline scan |

---

## Limitations

- **This figure plots % change from baseline, not from nadir.** RECIST 1.1's progression definition technically requires a >=20% increase *from nadir* (the smallest on-study value), which can differ from % change from baseline for a patient who has an intervening response before progressing. The dashed +20% line here should be read as a baseline-referenced visual guide, not a per-patient PD determination.
- **No censoring for new lesions or non-target progression** is shown; a patient can be classified PD due to new lesions while target-lesion SLD alone appears stable.

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
| `sld_over_time_plot.R` | Self-contained script · reads `Data/V1/ADTR.csv`, `Data/V1/ADSL.csv` |
| `Out/sld_over_time.png` | Output figure · 10 × 7.3 in · 300 DPI |

---

## References

Eisenhauer EA, Therasse P, Bogaerts J, et al. New response evaluation criteria in solid tumours: revised RECIST guideline (version 1.1). *Eur J Cancer.* 2009;45(2):228-247.

Litière S, Isaac G, De Vries EGE, et al. RECIST 1.1 for response evaluation apply not only to trials but also to clinical practice. *J Clin Oncol.* 2019;37(15):1289-1294 (discussion of baseline vs. nadir referencing conventions).
