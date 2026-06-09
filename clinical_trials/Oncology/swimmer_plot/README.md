# Swimmer Plot

A patient-level visualization displaying individual treatment duration, longitudinal response assessments, post-treatment follow-up, and outcome status (ongoing, death, discontinuation) per RECIST 1.1 criteria, with liver metastasis status, across a fully synthetic phase II/III basket trial dataset.

**Dataset:** ONCVIZ-001 · N = 263 (Treatment Arm) | **Cutoff:** 05 Mar 2026 | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

Swimmer plots are well established in oncology publications, yet existing open-source implementations share a consistent set of structural limitations. Most public examples use toy datasets of 15–36 patients with a single tumor type and no calibration to real trial benchmarks. The few multi-histology examples available draw patients from incompatible studies, making cross-domain comparison meaningless. Purpose-built packages such as `ggswim` (CHOP, 2025) and `swimplot` (PMH, 2022) solve the grammar problem elegantly, but ship with de-identified illustrative data that behaves nothing like a real oncology trial.

| | Prior art | This work |
|---|---|---|
| Dataset size | 15–36 synthetic patients (swimplot, ggswim vignettes) | 263 treatment-arm patients from ONCVIZ-001 |
| Tumor histologies | Single histology or uncalibrated mix | Five independently calibrated profiles |
| Cross-domain consistency | Swimmer and waterfall from different studies | Same 400 patients, 13 ADaM domains, 131,690 records |
| Follow-up encoding | Treatment bar only | Separate follow-up bar extending to OS date |
| Outcome symbols | Rarely distinguished | Arrow (ongoing) · ✕ (death) · ✕ red (discontinuation) |
| Response markers | Color only, or not time-located | Shape + color markers at exact assessment month |
| Liver mets annotation | Not shown | Y/N column per patient |
| Sorting | Duration only | BOR then descending duration within each BOR group |
| Reproducibility | No fixed seed, or seed not stated | `set.seed(42)` — bitwise-identical across R ≥ 4.1 |
| Calibration traceability | Unspecified | Anchored to KEYNOTE-189, IMbrave150, OlympiAD, KEYNOTE-177, NAPOLI-1 |

---

## Why the swimmer plot complements the waterfall plot

The waterfall plot answers *how much* — depth of tumor response at best assessment. The swimmer plot answers *how long* and *what happened next* — the temporal arc of each patient's journey from treatment start through follow-up or death. Used together from the same dataset they provide the complete RECIST efficacy picture that neither plot can deliver alone.

| Question | Waterfall | Swimmer |
|---|---|---|
| Depth of response (% SLD change) | ✓ | — |
| Duration of response | — | ✓ |
| Time to first response | — | ✓ |
| Post-treatment follow-up | — | ✓ |
| Outcome (ongoing / death) | — | ✓ |
| Liver metastasis status | — | ✓ |
| Responders vs. non-responders at a glance | ✓ | ✓ |

---

## Visual anatomy

```
RECIST   Liver
Response mets
  CR       N  │  [████████████████░░░░░]──────────────────── ▶  ongoing
  PD       Y  │  [████████████░░░░░░░░░░░░░░] ✕                  death
  SD       N  │  [████████████████████] ✕ (red)                   discontinuation
              │
              └────────────────────────────────────────────────
                   0m     6m     12m    18m    24m    36m
                           Months from Treatment Start
```

| Element | Description |
|---|---|
| Dark blue solid bar | Time on treatment |
| Grey translucent bar | Post-treatment follow-up to OS date (off-treatment, alive) |
| ▶ Filled triangle | Patient ongoing at data cutoff |
| ✕ Black × | Death |
| ✕ Red × | Discontinuation (not death) |
| Filled circle (●) | Complete Response (CR) — green |
| Filled square (■) | Partial Response (PR) — blue |
| Filled triangle (▲) | Stable Disease (SD) — orange |
| Filled diamond (◆) | Progressive Disease (PD) — red |
| Left text column | RECIST Response label, color-coded by BOR |
| Second text column | Liver metastasis status (Y / N), red for Y, blue for N |
| Summary box | Evaluable / Non-evaluable / Total Enrolled counts |
| Milestone lines | Dashed verticals at 6-month intervals |
| Alternating row shading | Guides the eye across long patient IDs |
| Footer | "RECIST 1.1 · Data cutoff: 05 Mar 2026" |

---

## Plot variants

### `swimmer_all_treatment.png` — All treatment-arm patients

All 263 treatment-arm patients in a single panel, sorted by BOR (CR → PR → SD → PD → NE) and by descending treatment duration within each BOR group. Evaluable and non-evaluable patients are separated by a horizontal black line. Provides a comprehensive overview of the full treatment arm.

- Figure dimensions: 16 × 11 in @ 120 DPI

### Per-histology panels (five files)

One plot per tumor type, each following the same layout as the all-patients view but scoped to a single histology. Allows direct comparison of response depth, duration, and outcome patterns within a tumor type without the density of the full-arm plot.

| File | Histology | Patients plotted |
|------|-----------|----------------:|
| `swimmer_nsclc.png` | NSCLC | 23 |
| `swimmer_brca.png` | BRCA | 8 |
| `swimmer_hcc.png` | HCC | 10 |
| `swimmer_crc.png` | CRC | 9 |
| `swimmer_pdac.png` | PDAC | 12 |

- Figure dimensions: 16 × 9 in @ 120 DPI each

> **Note:** Per-histology plots render only the evaluable subset of each tumor type rather than the full enrollment count. Patients with `EOSSTT` records but no evaluable RECIST assessment (`eval = FALSE`) are excluded from the per-histology outputs as currently filtered.

---

## Dataset — ONCVIZ-001 ADaM v1

```
Study    ONCVIZ-001 · Vizatinib 300 mg QD vs Placebo
Design   Phase II/III Open-Label Randomized Basket Trial
N        400  (TRT = 263 · CTL = 137 · 2:1 ratio)
Records  131,690 across 13 ADaM domains
Seed     42 · fully reproducible
Cutoff   05 Mar 2026
```

### Tumor-stratified response parameters (treatment arm)

| Histology | n | ORR TRT | ORR CTL | OS median | Calibration source |
|-----------|---|--------:|--------:|----------:|---|
| NSCLC | 88 | 47% | 20% | 19.6 m | KEYNOTE-189 (Gandhi et al., NEJM 2018) |
| HCC | 69 | 27% | 12% | 12.4 m | IMbrave150 (Finn et al., NEJM 2020) |
| CRC | 80 | 13% | 0% | 7.3 m | KEYNOTE-177 (André et al., NEJM 2020) |
| BRCA | 82 | 38% | 12% | 19.8 m | OlympiAD (Robson et al., NEJM 2017) |
| PDAC | 81 | 11% | 0% | 3.6 m | NAPOLI-1 / PRODIGE-4 |

### Domains used by this script

| Domain | Description | Rows | Key variables used |
|--------|-------------|-----:|---|
| ADSL | Subject-level | 400 | `USUBJID`, `ARM`, `TUMORTYPE`, `TRTSDT`, `TRTEDT`, `TRTDURD`, `EOSSTT`, `DCSREAS`, `OSCR`, `OSDTC`, `LIVERMETS` |
| ADRS | Tumor response per RECIST 1.1 | 1,211 | `PARAMCD="OVRLRESP"`, `ADT`, `AVALC` (CR/PR/SD/PD/NE) |

```r
# Load required ADaM domains
adsl <- read.csv(file.path(DATA_DIR, "ADSL.csv"), stringsAsFactors = FALSE)
adrs <- read.csv(file.path(DATA_DIR, "ADRS.csv"), stringsAsFactors = FALSE)

adsl$TRTSDT <- as.Date(adsl$TRTSDT)
adsl$TRTEDT <- as.Date(adsl$TRTEDT)
adsl$OSDTC  <- suppressWarnings(as.Date(adsl$OSDTC))

# Subset to treatment arm
trt <- adsl[adsl$ARM == "TREATMENT", ]
```

---

## Key implementation details

### Treatment bar vs. follow-up bar

Two separate bar layers encode distinct clinical phases. The treatment bar (dark blue, `#1F4E79`) runs from day 0 to treatment end or data cutoff for ongoing patients. The follow-up bar (grey `#808080`, narrower at 42% of bar height) runs from treatment end to the OS date when a patient discontinued but survived beyond treatment end. This separation makes it visually immediate whether a patient died on treatment, died after treatment, or is still alive.

```r
# Treatment duration
trt_mo <- max((as.numeric(trtedt - trtsdt)) / DAY2MO, 0.3)

# Follow-up bar (only when patient died after treatment end)
fu_mo <- NA_real_
if (!ong && !is.na(osdtc) && osdtc > trtedt)
  fu_mo <- as.numeric(osdtc - trtsdt) / DAY2MO
```

### Response markers

Each ADRS record is placed as a shaped, filled marker at its exact assessment month relative to treatment start. Only the first occurrence of each response category per patient is plotted to avoid overprinting. Shape encodes response category: circle = CR, square = PR, triangle = SD, diamond = PD.

```r
# Deduplicate: first occurrence of each response category per patient
seen <- character(0)
for (j in seq_len(nrow(pt_rs))) {
  rv <- pt_rs$AVALC[j]
  if (rv %in% names(RESP_COLORS) && !rv %in% seen) {
    # record marker; mark as seen
    seen <- c(seen, rv)
  }
}
```

### Outcome symbol placement

End-of-bar symbols are positioned past the rightmost bar edge (treatment bar or follow-up bar, whichever extends further): 0.55 months past for ongoing patients (filled triangle), 0.35 months past for death (black ×) and discontinuation (red ×).

```r
end_sym <- df |>
  mutate(
    x_e   = ifelse(!is.na(fu_mo) & fu_mo > trt_mo, fu_mo, trt_mo),
    sym_x = ifelse(end_t == "ongoing", x_e + 0.55, x_e + 0.35)
  )
```

### Liver metastasis annotation

The `LIVERMETS` variable from ADSL is displayed as a second left-margin column. Values are color-coded: red for Y, blue for N.

```r
LM_COLORS <- c(Y = "#C0392B", N = "#2B6CB0")
```

### Sorting

Patients are sorted by BOR rank (CR=0, PR=1, SD=2, PD=3, NE=4), then by descending treatment duration within each BOR group. Non-evaluable patients appear below a separating horizontal black line.

### Dynamic axis scaling

The x-axis upper limit and milestone interval are computed from the data:

```r
max_mo <- ceiling((max_data + 2) / 3) * 3
step   <- if (max_mo > 36) 6 else 3
```

This means per-histology plots that have shorter follow-up (e.g. PDAC) use 6-month intervals while plots with longer follow-up (e.g. BRCA, NSCLC) automatically extend to 84+ months.

---

## Output files

| File | Description | Dimensions |
|------|-------------|------------|
| `swimmer_plot.R` | Main R script — all plots | |
| `ADSL.csv` | Subject-level dataset | 400 rows |
| `ADRS.csv` | Overall response assessments | 1,211 rows |
| `swimmer_all_treatment.png` | All treatment-arm patients | 16 × 11 in · 120 DPI |
| `swimmer_nsclc.png` | NSCLC subgroup | 16 × 9 in · 120 DPI |
| `swimmer_brca.png` | BRCA subgroup | 16 × 9 in · 120 DPI |
| `swimmer_hcc.png` | HCC subgroup | 16 × 9 in · 120 DPI |
| `swimmer_crc.png` | CRC subgroup | 16 × 9 in · 120 DPI |
| `swimmer_pdac.png` | PDAC subgroup | 16 × 9 in · 120 DPI |

---

## When to use

**Appropriate:**
- Displaying duration of response alongside BOR in Phase I/II efficacy readouts
- Communicating the proportion of patients who remain on treatment at landmark timepoints
- Supplementary figures in clinical publications (commonly paired with waterfall plots)
- Regulatory dossiers requiring patient-level response timelines
- Basket trial figures where histology-specific duration patterns need to be visible side by side

**Limitations:**
- Does not show magnitude of tumor shrinkage — use waterfall plots
- Does not show continuous tumor size trajectories — use spider plots
- Becomes hard to read beyond ~100 patients per panel; consider subgroup filtering
- Sorting by BOR then duration means individual patients cannot be easily tracked across multiple figures
- Per-histology plots currently render the evaluable subset only; adjust `eval` filter in `sort_swimmer()` if full enrollment counts are needed

---

## Requirements

```
R >= 4.1
ggplot2 >= 3.4
dplyr >= 1.1
patchwork >= 1.2
grid
gridExtra
```

No proprietary packages required. All dependencies are available on CRAN.

---

## References

André T, et al. Pembrolizumab in microsatellite-instability–high advanced colorectal cancer. *N Engl J Med.* 2020;383(23):2207–2218.

Eisenhauer EA, et al. New response evaluation criteria in solid tumours: RECIST 1.1. *Eur J Cancer.* 2009;45(2):228–247.

Finn RS, et al. Atezolizumab plus bevacizumab in unresectable hepatocellular carcinoma. *N Engl J Med.* 2020;382(20):1894–1905.

Gandhi L, et al. Pembrolizumab plus chemotherapy in metastatic non-small-cell lung cancer. *N Engl J Med.* 2018;378(22):2078–2092.

Hanna R. ggswim: Create swimmer plots with ggplot2. Presented at R/Medicine 2025; Children's Hospital of Philadelphia. https://github.com/CHOP-CGTInformatics/ggswim

Robson M, et al. Olaparib for metastatic breast cancer in patients with a germline BRCA mutation. *N Engl J Med.* 2017;377(6):523–533.

Seymour L, et al. iRECIST: guidelines for response criteria for use in trials testing immunotherapeutics. *Lancet Oncol.* 2017;18(3):e143–e152.

Shalhout SZ, Miller DM. Graphical representation of survival: swimmer plots for clinical trials in oncology. *The Miller Lab.* 2020. https://www.themillerlab.io/posts/swimmer_plots/

Wickham H. ggplot2: Elegant Graphics for Data Analysis. Springer-Verlag New York, 2016.
