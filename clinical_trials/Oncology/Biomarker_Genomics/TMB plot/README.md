# TMB Plot — Tumor Mutational Burden by Tumor Type and Response

A patient-level dot plot of tumor mutational burden (TMB, mut/Mb) grouped by tumor type and colored by best overall response, with the clinical TMB-High threshold marked, for the ONCVIZ-001 synthetic basket trial.

**Dataset:** ONCVIZ-001 · n = 72 patients with mutation calls · ADMUT | **Endpoint:** Tumor Mutational Burden (mut/Mb) | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

TMB is a standard immuno-oncology biomarker, but it's usually shown either as a bare distribution (histogram/violin, losing the per-patient response link) or folded into a larger oncoprint where the continuous value is hard to read precisely. This plot keeps TMB as a continuous, patient-resolved value while layering in the two things a reviewer asks first: which tumor types run higher, and does response track with TMB.

| | Prior art | This work |
|---|---|---|
| Display | Histogram/violin per arm, or buried in oncoprint | Per-patient jittered dot plot, one point per patient |
| Grouping | Pooled or by arm only | By tumor type, ordered by descending median TMB |
| Response link | Not shown | Point color = best overall response (CR/PR/SD/PD) |
| Clinical threshold | Absent | TMB-High reference line (10 mut/Mb), with % of patients above it in the caption |
| Legend placement | Inside plot (risk of covering points) | Placed fully outside the plotting area — no data can ever sit underneath it |

---

## Visual anatomy

```
 TMB
(mut/Mb)
  50 |                              o PR (BRCA)
     |
  20 |        o PD                             *  <- legend
     |- - - - - - - - - - - - - - - - - -  TMB-High (10)        (outside axes)
   0 |  HCC    BRCA    NSCLC   PDAC    CRC
```

| Element | Description |
|---|---|
| Point | One patient, jittered horizontally within its tumor type |
| Point color | Best overall response (CR = dark blue, PR = light blue, SD = orange, PD = red) |
| Dashed horizontal line | TMB-High clinical threshold (10 mut/Mb) |
| X-axis order | Tumor types sorted by descending median TMB |

---

## Dataset variables used

| Variable | Source | Description |
|---|---|---|
| `USUBJID`, `TUMORTYPE`, `ARM` | ADMUT | Subject, tumor type, treatment arm |
| `TMB`, `TMBHIGH` | ADMUT | Tumor mutational burden (mut/Mb), TMB-High flag |
| `BESTRSPC` | ADMUT | Best overall response, used for point color |

---

## Statistical method

None — descriptive, one point per patient. `ADMUT` contains one row per called variant, so the script first reduces to one row per patient (`distinct(USUBJID, TUMORTYPE, ARM, TMB, TMBHIGH, BESTRSPC)`) before plotting, since TMB is a per-patient, not per-variant, value.

---

## Key parameters

| Parameter | Value |
|---|---|
| TMB-High threshold | 10 mut/Mb (standard clinical cutoff; adjust `tmb_high_thresh` if a study-specific cutoff applies) |
| Tumor-type order | Descending median TMB |
| Jitter width | 0.15 (horizontal only; y-values are unmodified) |

---

## Limitations

- **10 mut/Mb is a commonly used but not universal TMB-High cutoff** — some assays/trials use study-specific thresholds (e.g., FoundationOne CDx pan-tumor approval used 10 mut/Mb, but panel-dependent calibration varies). Confirm the assay-specific cutoff before using this threshold in a real analysis.
- **Small per-tumor-type N (8-25 patients)** — apparent differences in TMB distribution or response association across tumor types are descriptive only, not tested for significance here.
- **`BESTRSPC` on a per-patient TMB point does not imply causality** — this figure supports hypothesis generation about TMB-response association, not a biomarker-efficacy claim.

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
| `tmb_plot.R` | Self-contained script · reads `Data/V1/ADMUT.csv` |
| `Out/tmb_plot.png` | Output figure · 8.5 × 5.5 in · 300 DPI |

---

## References

Chalmers ZR, Connelly CF, Fabrizio D, et al. Analysis of 100,000 human cancer genomes reveals the landscape of tumor mutational burden. *Genome Med.* 2017;9(1):34.

Marabelle A, Fakih M, Lopez J, et al. Association of tumour mutational burden with outcomes in patients with select advanced solid tumours treated with pembrolizumab in KEYNOTE-158. *Lancet Oncol.* 2020;21(10):1353-1365.

FDA. Approval of pembrolizumab for TMB-High (>=10 mut/Mb) solid tumors (FoundationOne CDx companion diagnostic), 2020.
