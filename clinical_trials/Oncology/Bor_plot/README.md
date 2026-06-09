# Best Overall Response (BOR) Bar Chart

A cohort-level visualization displaying the distribution of RECIST 1.1–defined response categories (CR, PR, SD, PD) across dose cohorts — as both a stacked proportional bar (per-cohort breakdown) and a grouped bar (cross-cohort comparison) — with an overall ORR summary panel, derived from a fully synthetic phase I dose-escalation dataset.

**Dataset:** ONCVIZ-001 · N = 50 (Treatment Arm) | **Cutoff:** 05 Mar 2026 | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

Best Overall Response bar charts are ubiquitous in Phase I/II oncology publications, yet most open-source implementations reduce the figure to a single stacked bar over the full arm, losing the dose-response signal that is the primary scientific question of a dose-escalation trial. The few multi-cohort examples available use small toy datasets without calibration to real trial benchmarks, and none present both the per-cohort proportional view and the cross-cohort response-category comparison in a single coordinated figure.

| | Prior art | This work |
|---|---|---|
| Panel count | Single stacked bar | Two coordinated panels (stacked + grouped) |
| Granularity | Full-arm aggregate only | Per-cohort breakdown (4 dose levels) |
| Cohort labeling | Numeric dose only | Dose + RP2D annotation |
| Summary statistics | ORR text only | Formatted summary box with per-category counts |
| Patient counts | Rarely shown | `n=` labels above each cohort bar |
| Legend strategy | Duplicate per panel | Shared response and cohort legends, positioned once |
| Reproducibility | No fixed seed | Derived from `set.seed(42)` ADSL — bitwise-identical |

---

## Why the BOR bar chart complements the waterfall plot

The waterfall plot answers *how much* — continuous depth of tumor shrinkage per patient. The BOR bar chart answers *how many* — the aggregate distribution of categorized responses, at the cohort and arm level. Used together from the same dataset they deliver the complete RECIST efficacy picture.

| Question | Waterfall | BOR Bar |
|---|---|---|
| Individual % SLD change | ✓ | — |
| Proportion achieving CR / PR / SD / PD | — | ✓ |
| Dose–response gradient across cohorts | Partial | ✓ |
| ORR and category counts at a glance | — | ✓ |
| Responders vs. non-responders per cohort | — | ✓ |

---

## Visual anatomy

```
 Left panel (stacked)              Right panel (grouped)
 ─────────────────────────         ──────────────────────────────
  n=12   n=13   n=15   n=10         80%  ┤
  100% ┤ ░PD░  ░PD░  ░PD░  ░PD░         │  ■ Cohort 1
       │ ░SD░  ░SD░  ░SD░  ░SD░         │  ■ Cohort 2
       │ ░PR░  ░PR░  ░PR░  ░PR░   40%  ┤  ■ RP2D
       │ ░CR░  ░CR░  ░CR░  ░CR░         │  ■ Cohort 4
    0% ┤──────────────────────           └──────────────────────
        C1    C2   RP2D   C4            CR    PR    SD    PD
        (100) (200)(300) (400 mg)
 ═══════════════════════════════════════════════════════════════
  Total 50     ORR 48%  (24)     CR 8 (16%)   PR 16 (32%)  …
```

| Element | Description |
|---|---|
| Stacked bar (left panel) | Proportion of patients per response category within each cohort |
| Grouped bar (right panel) | Each response category plotted side-by-side across cohorts |
| `n=` label | Number of patients per cohort, displayed above each stacked bar |
| Percentage labels | White bold labels inside stacked segments ≥ 7% |
| Summary box | Navy bar showing total N, ORR%, and per-category count + percentage |
| Response legend | CR · PR · SD · PD with category colors (shared, right of left panel) |
| Cohort legend | Four dose cohorts with cohort colors (shared, right of right panel) |

---

## Objective Response Rate (ORR)

ORR is defined as the proportion of patients achieving CR or PR:

```
ORR = (CR + PR) / N × 100
```

This is the primary efficacy endpoint in most Phase II oncology trials and a key secondary endpoint in Phase I dose-escalation studies. The summary box at the bottom of the figure reports ORR as both a percentage and absolute count.

---

## Dose cohorts

| Cohort | Dose | Label |
|---|---|---|
| Cohort 1 | 100 mg | Standard escalation |
| Cohort 2 | 200 mg | Standard escalation |
| RP2D | 300 mg | Recommended Phase 2 Dose |
| Cohort 4 | 400 mg | Expansion / de-escalation |

---

## Dataset — ONCVIZ-001 ADaM v1

```
Study    ONCVIZ-001 · Vizatinib 300 mg QD vs Placebo
Design   Phase I Dose-Escalation
N        50  (Treatment Arm — ADSL filtered to ARM == "TREATMENT")
Seed     42 · fully reproducible
Cutoff   05 Mar 2026
```

### Domain used by this script

| Domain | Description | Key variables used |
|---|---|---|
| ADSL | Subject-level | `USUBJID`, `ARM`, `DOSELEVEL`, `BESTRSPC` |

```r
adsl <- read.csv(file.path(DATA_DIR, "ADSL.csv"), stringsAsFactors = FALSE)

trt <- adsl |>
  filter(ARM == "TREATMENT") |>
  mutate(COHORT = case_when(
    DOSELEVEL == 100 ~ "Cohort 1\n(100 mg)",
    DOSELEVEL == 200 ~ "Cohort 2\n(200 mg)",
    DOSELEVEL == 300 ~ "RP2D\n(300 mg)",
    DOSELEVEL == 400 ~ "Cohort 4\n(400 mg)"
  ))
```

---

## Key implementation details

### Two-panel layout with shared legends

The figure uses two coordinated `ggplot2` panels combined via `gridExtra::arrangeGrob`. Legends are extracted with a helper (`get_legend`) that pulls the `guide-box` grob from a temporarily modified copy of each plot — leaving both panels themselves legend-free and placing the response legend between the panels and the cohort legend to their right.

```r
get_legend <- function(p) {
  g   <- ggplotGrob(p)
  idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  g$grobs[[idx]]
}
```

### Percentage label suppression

Labels inside stacked segments are suppressed below 7% to avoid overprinting in narrow slices:

```r
geom_text(aes(label = ifelse(pct >= 7, paste0(pct, "%"), "")),
          position = position_stack(vjust = 0.5), ...)
```

### Summary box

The navy summary bar is constructed as a `rectGrob` / `textGrob` pair assembled via `arrangeGrob`, rendering a formatted string that includes total N, ORR (% and count), and per-category counts. Monospaced font (`fontfamily = "mono"`) keeps the columns aligned.

```r
summary_txt <- paste(
  sprintf("Total  %d", n_all),
  sprintf("ORR  %d%%  (%d)", orr_pct, orr_n),
  paste(resp_summary$lbl, collapse = "     "),
  sep = "     "
)
```

### Color palette

Response and cohort colors are defined as named character vectors and passed to `scale_fill_manual` independently in each panel:

```r
RESP_COLORS   <- c(CR="#1A9641", PR="#4A90C4", SD="#E8A020", PD="#C0392B")
COHORT_COLORS <- c("Cohort 1\n(100 mg)"="#E66101",
                   "Cohort 2\n(200 mg)"="#5AAE61",
                   "RP2D\n(300 mg)"    ="#3288BD",
                   "Cohort 4\n(400 mg)"="#9970AB")
```

### Y-axis headroom

The stacked bar panel uses `limits = c(0, 114)` to create headroom above 100% for the `n=` labels (placed at `y = 103`). The grouped bar panel uses `limits = c(0, 85)` scaled to the observed maximum per-category percentage.

---

## Output files

| File | Description | Dimensions |
|---|---|---|
| `bor_plot.R` | Main R script | |
| `ADSL.csv` | Subject-level dataset | 50 rows |
| `bor_plot.png` | Combined two-panel figure + summary box | 20 × 9 in · 150 DPI |

---

## When to use

**Appropriate:**
- Phase I/II dose-escalation trials: communicating dose–response relationship in best overall response
- Regulatory dossiers and clinical publications requiring aggregate RECIST response summaries
- Efficacy overview slides when the waterfall plot shows individual-patient data alongside this cohort-level summary
- Basket trials where the response distribution across tumor types or cohorts needs to be visible simultaneously

**Limitations:**
- Does not show individual patient tumor change — use waterfall plots
- Does not show response duration — use swimmer plots
- Percentage-based stacking can obscure absolute patient counts in small cohorts; always display `n=` labels
- ORR confidence intervals are not shown; add `binom.test` or Wilson intervals for publications

---

## Requirements

```
R >= 4.1
ggplot2 >= 3.4
dplyr >= 1.1
patchwork >= 1.2   # loaded but layout handled via gridExtra in this script
grid
gridExtra
```

No proprietary packages required. All dependencies are available on CRAN.

---

## References

Eisenhauer EA, et al. New response evaluation criteria in solid tumours: RECIST 1.1. *Eur J Cancer.* 2009;45(2):228–247.

Seymour L, et al. iRECIST: guidelines for response criteria for use in trials testing immunotherapeutics. *Lancet Oncol.* 2017;18(3):e143–e152.

Therasse P, et al. New guidelines to evaluate the response to treatment in solid tumors. *J Natl Cancer Inst.* 2000;92(3):205–216.

Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.

Zhu AX, et al. Pembrolizumab in combination with gemcitabine and cisplatin compared with gemcitabine and cisplatin alone for patients with advanced biliary tract cancer (KEYNOTE-966). *Lancet.* 2023;401(10391):1853–1865.
