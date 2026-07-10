# Waterfall Plot + PK Overlay

A patient-level visualization displaying best percentage change from baseline in tumor size per RECIST 1.1 criteria, colored by best overall response and sorted descending, with a per-patient Cycle 1 AUC exposure heatmap strip aligned beneath the bars, from a fully synthetic phase I/II oncology trial dataset.

**Dataset:** ONCVIZ-001 · N = 62 (treatment arm) | **Cutoff:** 05 March 2026 | **Language:** R · ggplot2 · patchwork | **License:** CC BY 4.0

---

## The gap this fills

Waterfall plots with a PK overlay commonly implement the overlay as a second y-axis line graph crossing through the bars — visually confusing and not the convention actually used in oncology PK/PD reporting (e.g., Clin Cancer Res, JCO). Threshold labels and legends are also frequently placed directly on top of bars, obscuring the response category color underneath.

| | Prior art | This work |
|---|---|---|
| PK overlay style | Second y-axis line crossing bars | Heatmap strip beneath bars, aligned to the same patient order |
| Threshold labels | Placed on the dashed line at any x, often over bars | Placed in a verified-empty interior region (x=26, where bars sit near 0%) |
| Legend | Frequently overlapping tall bars | Placed above the plot, entirely outside the bar area |
| Response color | Ad hoc per figure | Same `RESP_COLORS` palette as every other plot in this repo and its Waterfall_plot precedent |
| Missing data | Silently dropped | Patients with no post-baseline assessment or no Cycle 1 AUC explicitly counted in caption |

---

## Visual anatomy

```
 40% |██
     |███  ██
   0 |░░░░░░░██░  ← PD threshold (+20%), label in empty interior space
     |        ░░████░░
 -30%|. . . . . . . . .  ← PR threshold (−30%), label in empty interior space
     |               ░░░████
-100%|                    ░████
     ──────────────────────────
      ▓▓▓▓░░██▓▓░░████▓▓░░██▓▓   ← PK overlay strip: AUC heatmap, same patient order
```

| Element | Description |
|---|---|
| Bar | Best % change in tumor size from baseline, one bar per patient |
| Bar color | RECIST 1.1 best overall response — navy (CR), light blue (PR), amber (SD), red (PD) |
| Dashed lines | RECIST thresholds at −30% (PR) and +20% (PD) |
| Threshold labels | Positioned at x≈26 (verified empty region), never over a bar |
| Legend | Above the plot, response counts shown per category |
| Bottom strip | AUC₀₋₂₄ heatmap (log color scale), aligned to the same patient rank order as the bars above |

---

## Why a heatmap strip, not a second y-axis

A second y-axis line graph crossing through waterfall bars is visually ambiguous — readers cannot tell whether a line crossing a bar means anything about that bar's value. The heatmap-strip convention (color intensity per patient, same x-order as the bars) is the standard seen in combined PK/PD oncology figures and avoids this ambiguity entirely.

```r
p_pk <- ggplot(wf_plot, aes(x = x, y = 1, fill = AUC)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("#FFFFCC", "#FD8D3C", "#800026"), trans = "log10")
```

---

## Label placement rule

Threshold labels are placed at a data coordinate confirmed empty by inspection of the actual sorted values (bars near 0% around rank 20–35), rather than at a fixed axis-relative position that risks colliding with tall bars on either end of the distribution:

```r
annotate("label", x = 26, y = PD_TH, label = "PD threshold (+20%)", ...)
annotate("label", x = 26, y = PR_TH, label = "PR threshold (\u201230%)", ...)
```

---

## Data caveat

A small number of patients with `BESTRSPC = "CR"` in ADSL have shallow target-lesion `PCHG` values (e.g., −2.3%) that would not meet RECIST 1.1 CR criteria based on the target-lesion sum alone (CR requires disappearance of all lesions). This is a known synthetic-data cross-domain inconsistency between the categorical BOR field and the target-lesion measurement domain — not a plotting artifact. Real ADaM datasets should have these two fields concordant by construction.

---

## Dataset — ONCVIZ-001 ADaM v1

| Domain | Description | Key variables |
|--------|-------------|---|
| ADSL | Subject-level | `USUBJID`, `ARM`, `BESTRSPC` |
| ADTR | Tumor measurements | `AVISITN`, `PCHG` (best post-baseline value used) |
| ADPK | AUC₀₋₂₄, Cycle 1 Day 1 | `PARAMCD="AUC"`, `AVAL` |

```r
trt  <- adsl |> filter(ARM == "TREATMENT") |> select(USUBJID, BESTRSPC)
best <- adtr |> filter(AVISITN > 0) |> group_by(USUBJID) |> summarise(BEST_PCHG = min(PCHG, na.rm = TRUE))
auc  <- adpk |> filter(PARAMCD == "AUC", AVISIT == "CYCLE 1 DAY 1") |> distinct(USUBJID, .keep_all = TRUE)
```

---

## Output files

| File | Description | Dimensions |
|------|-------------|------------|
| `waterfall_pk_overlay.R` | Main R script | |
| `Out/waterfall_pk_overlay.png` | Waterfall + PK heatmap overlay | 13×7.5 in · 300 DPI |

---

## When to use

**Appropriate:**
- Communicating RECIST response distribution alongside individual PK exposure in one figure
- Phase I/II efficacy signal reporting with a PK exposure hypothesis
- Supplementary figures in clinical publications combining PK and tumor response

**Limitations:**
- Does not show duration of response — use swimmer plots
- Rank-order x-axis is not a measurement scale — do not interpret bar spacing as time or exposure difference
- Requires a per-patient scalar PK exposure metric (AUC or Cmax); not designed for full concentration-time overlays

---

## Requirements

```
R >= 4.1
```

```r
install.packages(c("dplyr", "ggplot2", "patchwork"))
```

---

## References

Eisenhauer EA, et al. New response evaluation criteria in solid tumours: RECIST 1.1. *Eur J Cancer.* 2009;45(2):228–247.

U.S. Food and Drug Administration. *Exposure-Response Relationships — Study Design, Data Analysis, and Regulatory Applications.* CDER; 2003.

Gabrielsson J, Weiner D. *Pharmacokinetic and Pharmacodynamic Data Analysis: Concepts and Applications.* 5th ed. Swedish Pharmaceutical Press; 2016.

