# Target Lesion Change Plot (Waterfall Plot)

A patient-level waterfall plot of best percent change from baseline in target-lesion sum of longest diameters, sorted descending, colored by best overall response per RECIST 1.1, for all evaluable patients in ONCVIZ-001.

**Dataset:** ONCVIZ-001 · n = 80 evaluable patients · ADTR / ADSL | **Endpoint:** Best % change in target-lesion SLD (RECIST 1.1) | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

The waterfall plot is the single most recognizable oncology efficacy graphic, which makes departures from convention (mislabeled thresholds, unlabeled response colors, labels colliding with bars) especially visible to an expert reader.

| | Prior art | This work |
|---|---|---|
| Sort order | Sometimes by patient ID | Strictly descending by best % change (standard convention) |
| Response color | Often absent | CR/PR/SD/PD color-coded, with n per category in the legend |
| Threshold labels | Placed inside the bar field, prone to overlap | Placed outside the bar range, to the right of the last bar |
| Best % change definition | Sometimes ambiguous (baseline vs. nadir) | Explicitly the minimum post-baseline value (nadir), stated in the caption |

---

## Visual anatomy

```
 +40 │██
     │███
     │████                                              ██
   0 ┼─────────────────────────────────────────────  +20% (PD)
     │                    ▁▂▃▄▅▆▇█                  -30% (PR)
 -100│                                    ██████████████
      ← sorted: greatest increase (left) to greatest decrease (right)
```

| Element | Description |
|---|---|
| Bar | One patient, height = best % change in target-lesion SLD |
| Bar color | CR = dark blue, PR = light blue, SD = orange, PD = red |
| Dashed red line (+20%) | RECIST 1.1 PD threshold |
| Dashed blue line (−30%) | RECIST 1.1 PR threshold |
| Threshold labels | Placed clear of the bar field, to the right of the plot |

---

## Dataset variables used

| Variable | Source | Description |
|---|---|---|
| `USUBJID`, `AVISITN`, `PCHG` | ADTR (`PARAMCD = "SUMDIAM"`) | Subject, visit number, % change from baseline |
| `BESTRSPC` | ADSL | Best overall response category |

---

## Statistical method

For each patient, the **best % change** is the minimum (most negative) post-baseline `PCHG` value across all target-lesion assessments (`AVISITN > 0`) — the standard "best response" nadir definition used in oncology waterfall plots. Patients are then sorted descending by this value.

---

## Key parameters

| Parameter | Value |
|---|---|
| PD threshold | +20% |
| PR threshold | −30% |
| Best-change definition | Minimum post-baseline PCHG (nadir) |
| Evaluable population | Patients with ≥1 post-baseline target-lesion assessment (n = 80) |

---

## Limitations

- **Bar color (BESTRSPC) and the bar height (nadir PCHG) can occasionally disagree** — e.g., a patient whose nadir crosses −30% but whose confirmed best overall response is SD (unconfirmed response, or a later scan reversing the reduction) will show a blue-threshold-crossing bar colored orange. This is expected: `BESTRSPC` reflects the fully adjudicated RECIST 1.1 response (including confirmation and non-target/new-lesion assessment), while the bar height reflects target-lesion measurement alone.
- **No marker for new lesions or non-evaluable patients** (baseline SLD = 0) is included in this version.

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
| `target_lesion_waterfall_plot.R` | Self-contained script · reads `Data/V1/ADTR.csv`, `Data/V1/ADSL.csv` |
| `Out/target_lesion_waterfall.png` | Output figure · 12 × 6.5 in · 300 DPI |

---

## References

Eisenhauer EA, Therasse P, Bogaerts J, et al. New response evaluation criteria in solid tumours: revised RECIST guideline (version 1.1). *Eur J Cancer.* 2009;45(2):228–247.

Karrison TG, Maitland ML, Stadler WM, Ratain MJ. Design of phase II cancer trials using a continuous endpoint of change in tumor size: application to a study of sorafenib and erlotinib in non-small-cell lung cancer. *J Natl Cancer Inst.* 2007;99(19):1455–1461.
