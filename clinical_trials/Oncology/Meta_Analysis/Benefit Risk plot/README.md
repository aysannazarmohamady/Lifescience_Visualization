# Benefit-Risk Plot — Quadrant Framework

A quadrant scatter plot positioning ORR (benefit) against Grade >=3 AE rate (risk), by tumor type and overall arm, for the ONCVIZ-001 synthetic basket trial — an informal multi-criteria decision-analysis (MCDA) style visualization.

**Dataset:** ONCVIZ-001 · 5 tumor types + Overall Treatment/Control · ADSL / ADAE | **Endpoints:** ORR (benefit), Grade >=3 AE rate (risk) | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

Unlike CONSORT, RECIST, or forest plots, there is no single, universally adopted "classic" benefit-risk chart — regulators accept several frameworks (BRAT, PrOACT-URL, MCDA value trees). Most implementations either invent an ungrounded scoring scale or fail to say which framework they follow, making the figure hard to defend to a reviewer.

| | Prior art | This work |
|---|---|---|
| Framework | Often unstated | Explicitly named: quadrant scatter / informal MCDA; alternative frameworks (BRAT, PrOACT-URL) noted in the caption |
| Uncertainty | Point estimates only | 95% Wilson score CIs on both axes |
| Reference lines | Absent | Cross-group median lines delineating quadrants |
| Quadrant labeling | Absent | "More favorable" / "less favorable" corners explicitly labeled, placed clear of data |
| Legend | Cramped | Dedicated panel with vertical spacing between entries |

---

## Visual anatomy

```
 ORR(%)                     o NSCLC
  50 |        More            /
     |      favorable        o CRC
     |                    o BRCA
  25 |  [] Control   --------------- median
     |              o PDAC   o HCC
   0 |                    Less favorable
     20    40    60    80   100
       Grade>=3 AE Rate (%)
```

| Element | Description |
|---|---|
| Circle | Tumor-type subgroup (Treatment arm only) |
| Square | Overall Treatment / Overall Control |
| Error bars | 95% Wilson score CI, both axes |
| Gray reference lines | Cross-group median benefit and median risk |

---

## Dataset variables used

| Variable | Source | Description |
|---|---|---|
| `USUBJID`, `TUMORTYPE`, `ARM`, `BESTRSPC` | ADSL | Tumor type, arm, best overall response (for ORR numerator) |
| `AETOXGR` | ADAE | CTCAE grade, filtered to >=3 (for risk numerator) |

---

## Statistical method

**Benefit** = Objective Response Rate = (CR + PR) / N, with a **Wilson score 95% CI** (more appropriate than the normal approximation for the small subgroup sizes here, n = 8-23 per tumor type). **Risk** = proportion of patients with >=1 Grade >=3 AE, same Wilson CI method. No formal weighting or utility function is applied — this is a **quadrant visualization**, not a scored MCDA output.

```
wilson_ci(k, n):  (p + z^2/2n +/- z*sqrt(p(1-p)/n + z^2/4n^2)) / (1 + z^2/n),  z = 1.96
```

---

## Key parameters

| Parameter | Value |
|---|---|
| CI method | Wilson score, 95% |
| Benefit endpoint | ORR (CR + PR), RECIST 1.1 |
| Risk endpoint | Grade >=3 AE incidence (CTCAE) |
| Reference lines | Cross-group medians (not clinical thresholds) |

---

## Limitations

- **No formal weighting or utility scoring.** This is the simplest defensible benefit-risk visualization (a labeled quadrant scatter); if a regulatory-grade multi-criteria assessment is required, use a named framework such as **BRAT** (Benefit-Risk Action Team) or **PrOACT-URL**, which require explicit criterion weights and stakeholder-elicited value trees — confirm the expected framework with stakeholders before this figure is used for a formal submission.
- **Small subgroup sizes** (n = 8-23 per tumor type) produce wide CIs; the point estimates should not be over-interpreted as tumor-type-specific benefit-risk conclusions.
- **Single-timepoint snapshot.** ORR and Grade >=3 AE rate are both cumulative-to-cutoff proportions; a longitudinal benefit-risk view (e.g., against exposure time) is not shown here.

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
| `benefit_risk_plot.R` | Self-contained script · reads `Data/V1/ADSL.csv`, `Data/V1/ADAE.csv` |
| `Out/benefit_risk_plot.png` | Output figure · 11.5 × 8 in · 300 DPI |

---

## References

Wilson EB. Probable inference, the law of succession, and statistical inference. *J Am Stat Assoc.* 1927;22(158):209-212. (Wilson score interval.)

FDA. *Benefit-Risk Assessment for New Drug and Biological Products* (BRAT framework guidance), Center for Drug Evaluation and Research, 2018.

IMI PROTECT Benefit-Risk Work Package. PrOACT-URL and MCDA methodology reports for benefit-risk assessment of medicines. European Medicines Agency-affiliated IMI PROTECT consortium, 2013.
