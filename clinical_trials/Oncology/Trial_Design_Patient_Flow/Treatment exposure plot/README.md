# Treatment Exposure Plot (Swimmer Plot)

A patient-level swimmer plot showing time on treatment for each patient in the treatment arm of ONCVIZ-001, annotated with best overall response, first CR/PR, progression events, and ongoing-treatment status.

**Dataset:** ONCVIZ-001 · Treatment arm, n = 62 · ADSL / ADRS | **Endpoint:** Treatment duration & response timing | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

Swimmer plots are visually simple but implementations frequently omit one of the three elements regulators expect together: response color-coding, event markers (progression, response), and an explicit ongoing-treatment indicator.

| | Prior art | This work |
|---|---|---|
| Bar color | Often uncolored or by dose only | Colored by best overall response (CR/PR/SD/PD) |
| Response marker | Frequently omitted | Star marker at first confirmed CR/PR |
| Progression marker | Frequently omitted | X marker at first PD assessment |
| Ongoing status | Bar simply ends | Arrowhead appended for patients still on treatment |
| Sort order | Patient ID | Exposure duration (ascending), a CONSORT/CSR convention |

---

## Visual anatomy

```
 Pt80  ████████████████████████████████████████████████▶  (ongoing, CR)
 Pt79  ██████████████████████████████████████████████
 Pt78  ★ ✕ ██████████████████████████████████████▶
  ...
       0        10        20        30        40   months
```

| Element | Description |
|---|---|
| Bar color | CR = dark blue, PR = light blue, SD = orange, PD = red |
| ★ (gold star) | First confirmed CR or PR |
| ✕ (black cross) | First progressive-disease assessment |
| Arrowhead | Patient remains on treatment as of data cutoff |

---

## Dataset variables used

| Variable | Source | Description |
|---|---|---|
| `USUBJID`, `TRTDURD`, `EOSSTT`, `BESTRSPC` | ADSL | Duration on treatment, disposition, best overall response |
| `AVALC`, `ADTN` | ADRS (`PARAMCD = "OVRLRESP"`) | Response category and day, per assessment |

---

## Statistical method

None — descriptive display of patient-level timelines. Duration in months = `TRTDURD / 30.44`. Response/progression markers use the earliest assessment day meeting each category.

---

## Key parameters

| Parameter | Value |
|---|---|
| Sort order | Ascending treatment duration |
| Time unit | Months (days ÷ 30.44) |
| Population | Treatment arm only (n = 62) |

---

## Limitations

- **Tumor-assessment data window is shorter than treatment duration.** Response/progression markers in this dataset cluster within the first ~4 months (through Cycle 3 Day 1) even for patients whose bars extend past 40 months, because longitudinal `ADRS` records were only generated through the early assessment cycles. This is a property of the synthetic dataset's assessment schedule, not a plotting defect — later response changes for long-duration patients are not captured here.
- **No confirmation window enforcement.** "First CR/PR" uses the first record coded CR/PR in `ADRS`, without independently re-verifying RECIST 1.1's confirmation-scan requirement.

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
| `treatment_exposure_plot.R` | Self-contained script · reads `Data/V1/ADSL.csv`, `Data/V1/ADRS.csv` |
| `Out/treatment_exposure_plot.png` | Output figure · 9 in wide × variable height (~0.14 in/patient) · 300 DPI |

---

## References

Eisenhauer EA, Therasse P, Bogaerts J, et al. New response evaluation criteria in solid tumours: revised RECIST guideline (version 1.1). *Eur J Cancer.* 2009;45(2):228–247.

FDA Oncology Center of Excellence. Considerations for graphical presentation of time-to-event and exposure data in oncology submissions (swimmer plot conventions in NDA/BLA review documents).
