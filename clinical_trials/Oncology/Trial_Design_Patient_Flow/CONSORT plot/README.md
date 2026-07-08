# CONSORT Diagram — Participant Flow

A CONSORT 2010-compliant participant flow diagram tracing patients from screening through analysis for the ONCVIZ-001 synthetic randomized basket trial: enrollment, allocation, follow-up, and analysis populations.

**Dataset:** ONCVIZ-001 · Screened (n = 92) → Randomized (n = 80) · ADRAND / ADSL | **Endpoint:** Study flow / disposition | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

CONSORT flow diagrams are near-universally hand-built in PowerPoint or Visio, which breaks reproducibility — the numbers in the boxes are typed once and never re-derived if the underlying dataset changes. Programmatic ggplot2 implementations exist but often improvise box placement, leaving the Excluded box floating without a properly touching connector, or misplacing the Randomized box relative to the flow above it.

| | Prior art | This work |
|---|---|---|
| Box connectors | Approximate, visual gaps common | Every arrow starts and ends exactly on a box edge (computed, not eyeballed) |
| Excluded box position | Side branch or ad hoc placement | Centered in the same vertical column as Assessed/Randomized, per user layout preference |
| Numbers | Manually typed | Derived directly from `ADRAND`/`ADSL` flags (`SFFL`, `RANDFL`, `ITTFL`, `SAFFL`, `DCSREAS`) |
| Reconciliation | Not guaranteed | Subgroup counts sum to parent box by construction |
| Reproducibility | Static image | Re-run against updated data reproduces an updated diagram |

---

## Visual anatomy

```
        Assessed for eligibility (n = 92)
                     │
                     ▼
        Excluded (n = 12) — reasons listed
                     │
                     ▼
             Randomized (n = 80)
              /                \
             ▼                  ▼
   Allocated: Treatment    Allocated: Control
      (n = 62)                 (n = 18)
             │                  │
             ▼                  ▼
   Discontinued (n = 38)   Discontinued (n = 16)
             │                  │
             ▼                  ▼
     Analyzed, ITT/Safety   Analyzed, ITT/Safety
```

| Element | Description |
|---|---|
| Boxes | Rounded rectangles, black border, white fill, no color coding (per CONSORT convention) |
| Arrows | Vertical/diagonal, always terminate exactly on a box's edge |
| Reason lists | Bulleted, sorted by descending frequency |

---

## Dataset variables used

| Variable | Source | Description |
|---|---|---|
| `SFFL` / `SFREASN` | ADRAND | Screen-failure flag and reason |
| `RANDFL` | ADRAND | Randomization flag |
| `ARM` | ADSL | "TREATMENT" · "CONTROL" |
| `TRTSDT` | ADSL | First dose date (used as "received intervention" proxy) |
| `EOSSTT` / `DCSREAS` | ADSL | End-of-study status and discontinuation reason |
| `ITTFL` / `SAFFL` | ADSL | Intent-to-treat and safety population flags |

---

## Statistical method

None — this is a descriptive flow diagram. All figures are direct counts or `sum()`/`count()` aggregations of the flags above; no modeling is involved.

---

## Key parameters

| Parameter | Value |
|---|---|
| Screened | n = 92 |
| Screen-failed | n = 12 |
| Randomized | n = 80 (62 Treatment : 18 Control, ~3.4:1) |
| Canvas | 0–100 arbitrary coordinate grid, `theme_void()` |

---

## Limitations

- **"Received intervention" is proxied by a non-missing `TRTSDT`.** ADaM does not include a dedicated `RECIEVED` flag in this dataset version; a patient who received a single partial dose is still counted as "received."
- **No ineligible-post-randomization category.** This trial design does not model post-randomization eligibility exclusions.
- **Static reasons list.** If a discontinuation reason has zero patients in one arm, it is simply omitted from that arm's box (asymmetric reason lists between arms are expected and correct here).

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
| `consort_plot.R` | Self-contained script · reads `Data/V1/ADRAND.csv`, `Data/V1/ADSL.csv` |
| `Out/consort_diagram.png` | Output figure · 11 × 16 in · 300 DPI |

---

## References

Schulz KF, Altman DG, Moher D; CONSORT Group. CONSORT 2010 Statement: updated guidelines for reporting parallel group randomised trials. *BMJ.* 2010;340:c332. Also published in *Ann Intern Med.* 2010;152(11):726–732.

Moher D, Hopewell S, Schulz KF, et al. CONSORT 2010 Explanation and Elaboration: updated guidelines for reporting parallel group randomised trials. *BMJ.* 2010;340:c869.
