# Toxicity Heatmap — Worst Grade by System Organ Class

A patient-level heatmap of worst CTCAE grade reached per System Organ Class (SOC), Treatment arm, for the ONCVIZ-001 synthetic basket trial — a complementary safety view to the AE Bar Chart, organized by grade severity rather than by individual Preferred Term.

**Dataset:** ONCVIZ-001 · Treatment arm, n = 62 · ADAE / ADSL | **Endpoint:** Worst-grade AE incidence by SOC (CTCAE) | **Language:** R · ggplot2 · tidyr | **License:** CC BY 4.0

---

## The gap this fills

The AE Bar Chart in this repo reports incidence per individual Preferred Term, which is precise but does not answer a common first question in a safety review: *which organ systems carry the most severe toxicity, and at what grade?* A heatmap answers that in one glance, at the SOC level, without the ~25-term list.

| | Prior art | This work |
|---|---|---|
| Granularity | Preferred Term (AE Bar Chart, existing) | System Organ Class (7 SOCs) |
| Metric shown | Any-grade / Grade ≥3 split | Full grade distribution (G1–G5) per SOC |
| Aggregation | Any occurrence counted | **Worst grade per patient per SOC** (a patient is counted once, at their most severe grade for that SOC) |
| Cell annotation | N/A | % of patients labeled directly on each cell |
| Row order | N/A | Sorted by total safety burden (ascending) |

---

## Visual anatomy

```
                        G1   G2   G3   G4   G5
 General disorders      18   47   31    5
 GI disorders            11   48   37    3
 Skin and subcutaneous   31   37   16    8
 ...
```

| Element | Description |
|---|---|
| Row | System Organ Class (AESOC) |
| Column | CTCAE grade (G1–G5) |
| Cell color | % of Treatment-arm patients whose **worst** grade for that SOC equals the column |
| Cell text | Percentage value (blank if 0%) |

---

## Dataset variables used

| Variable | Source | Description |
|---|---|---|
| `USUBJID`, `AESOC`, `AETOXGR`, `ARM` | ADAE | Subject, System Organ Class, CTCAE grade, arm |
| `CTCAE_V` | ADAE | CTCAE version used for grading |

---

## Statistical method

For each patient and SOC, the **worst (maximum) grade** across all AE records in that SOC is taken:

```r
worst_grade = max(AETOXGR)  # per USUBJID x AESOC
pct = 100 * n_patients_at_worst_grade / N_treatment_arm
```

This avoids double-counting a patient who had, e.g., both a Grade 1 and a Grade 3 event in the same SOC — they contribute only to the Grade 3 cell.

---

## Key parameters

| Parameter | Value |
|---|---|
| Population | Treatment arm only (n = 62) |
| Grade levels | CTCAE Grade 1–5 |
| Row sort | Ascending total % burden across all grades |
| Color scale | Sequential (white → dark red), linear in % |

---

## Limitations

- **Treatment-arm only.** A Control-arm panel is not included by default; extend the script by adding a `facet_wrap(~ARM)` if an arm comparison heatmap is needed.
- **Worst-grade aggregation hides multiplicity.** A patient with five Grade 2 GI events and one Grade 1 event appears only once, in the Grade 2 cell — frequency information is intentionally traded for a cleaner severity summary.
- **SOC-level granularity loses the specific Preferred Term.** Use alongside the AE Bar Chart (`Safety_Toxicity/AE plot/`) when Preferred Term detail is needed.

---

## Requirements

```r
library(ggplot2)    # >= 3.4
library(dplyr)      # >= 1.1
library(tidyr)      # >= 1.3
```

---

## Files

| File | Description |
|---|---|
| `toxicity_heatmap.R` | Self-contained script · reads `Data/V1/ADAE.csv`, `Data/V1/ADSL.csv` |
| `Out/toxicity_heatmap.png` | Output figure · 8 × 6 in · 300 DPI |

---

## References

U.S. Department of Health and Human Services, National Cancer Institute. *Common Terminology Criteria for Adverse Events (CTCAE) v5.0.* 2017.

ICH E3: *Structure and Content of Clinical Study Reports*, Section 12 (Safety Evaluation) — SOC-level toxicity summarization conventions.
