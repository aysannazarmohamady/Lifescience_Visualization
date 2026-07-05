# ctDNA Dynamics Plot — Longitudinal ctDNA VAF Change by Best Overall Response

A publication-ready faceted longitudinal plot showing circulating tumor DNA (ctDNA) VAF percent change from baseline across up to 7 on-treatment timepoints (Day 0–168) in 62 treatment-arm patients, stratified by best overall response (BOR). Individual patient trajectories are shown as thin semi-transparent lines; group mean ± 95% CI is overlaid as a bold line with ribbon.

**Dataset:** ONCVIZ-001 · Treatment Arm · N = 62 patients · ADBM (PARAMCD = CTDNA) | **Language:** R · ggplot2 | **License:** CC BY 4.0

---

## The gap this fills

ctDNA dynamics plots are increasingly required in oncology publications to demonstrate early pharmacodynamic response and correlate molecular response with radiological BOR. Most implementations show either (1) a single spaghetti plot of all patients without BOR stratification, or (2) group means only without individual trajectory context. The combination of individual trajectories (thin lines, low alpha) with group mean ± 95% CI ribbon per BOR category in a synchronized facet layout is the emerging standard for regulatory and high-impact journal submissions.

| | Prior art | This work |
|---|---|---|
| Individual trajectories | Often omitted | Thin semi-transparent lines per patient |
| Group summary | Mean only | Mean ± 95% CI ribbon |
| BOR stratification | Separate figures | Synchronized `facet_wrap(nrow=1)` |
| BOR order | Alphabetical | Clinical ordering: CR → PR → SD → PD |
| Facet label | BOR code only | BOR code + "(n = x)" |
| Reference line | Not shown | Horizontal dashed at 0% change |
| Timepoints | Day or cycle | Days (0, 21, 42, 63, 84, 126, 168) |
| Reproducibility | No fixed seed | `set.seed(13)` · bitwise-identical |

---

## Visual anatomy

```
  ctDNA % change from baseline
  +180 ┤       ╱ (individual PD patients)
  +100 ┤      ╱
     0 ┤──────────────── (zero reference, dashed)
   −50 ┤  ╲              ╲ (PR group mean)
  −100 ┤   ╲──────────────
  −150 ┤    ╲─────── (CR group mean)
       └────────────────────── Day
        0  21  42  63  84 126 168

  [CR (n=x)] [PR (n=x)] [SD (n=x)] [PD (n=x)]
  ← facets, synchronized x-axis →
```

| Element | Description |
|---|---|
| Thin lines | Individual patient trajectories · `alpha = 0.18` |
| Ribbon | Group 95% CI · `geom_ribbon(alpha = 0.45)` |
| Bold line | Group mean · `linewidth = 2.2` |
| Points on bold line | Mean at each timepoint · `size = 2.8` |
| Dashed horizontal | 0% change reference |
| Facet | BOR category · labeled with n |

---

## Output figure

| File | Dimensions | DPI |
|---|---|---|
| `plots/13_ctdna.png` | 17 × 7 in | 180 |

---

## Synthetic dataset — ONCVIZ-001 ADaM v1 (ADBM)

```
Arm            Treatment (n = 62)
PARAMCD        CTDNA
Timepoints     Day 0 · 21 · 42 · 63 · 84 · 126 · 168 (7 visits)
Endpoint       VAF percent change from baseline (PCHG)
BOR groups     CR · PR · SD · PD
Seed           13 · fully reproducible
```

### BOR group parameters (synthetic trajectory generation)

| BOR | Base trend | Interpretation |
|---|---|---|
| CR | −80% × (day/168) | Steep sustained decline → complete molecular response |
| PR | −50% × (day/168) | Moderate decline → partial response |
| SD | +5% × (day/168) | Near-flat → stable disease |
| PD | +60% × (day/168) | Rising → progressive disease |

Each patient's trajectory adds `N(0, 20)` noise to the base trend per timepoint.

### BOR colors

| BOR | Color |
|---|---|
| CR | `#1a3f8f` (Dark blue) |
| PR | `#4f8fd4` (Light blue) |
| SD | `#e8a020` (Amber) |
| PD | `#c0392b` (Red) |

### Variables used

| Variable | Source | Description |
|---|---|---|
| `USUBJID` | ADBM | Subject identifier |
| `PARAMCD` | ADBM | Filtered to "CTDNA" |
| `ARM` | ADBM | Filtered to "TREATMENT" |
| `ADTN` | ADBM | Study day (numeric) |
| `PCHG` | ADBM | % change from baseline ctDNA VAF |
| `BESTRSPC` | ADBM / ADSL | Best overall response (CR · PR · SD · PD) |
| `ANL01FL` | ADBM | Analysis flag (filtered to "Y") |

---

## Statistical method

Group summary statistics are computed per BOR per timepoint:

```r
ctdna_sum <- ctdna %>%
  group_by(BESTRSPC, ADTN) %>%
  summarise(
    mean = mean(PCHG, na.rm = TRUE),
    se   = sd(PCHG, na.rm = TRUE) / sqrt(n()),
    n    = n()
  ) %>%
  mutate(
    lo = mean − qt(0.975, n − 1) × se,
    hi = mean + qt(0.975, n − 1) × se
  )
```

The 95% CI uses the t-distribution (two-sided, α = 0.05). With small group sizes (n < 5 per BOR per timepoint), the CI will be very wide — this is expected and clinically meaningful.

---

## Key parameters

| Parameter | Value | Description |
|---|---|---|
| Individual line alpha | 0.18 | Near-transparent individual trajectories |
| Individual line width | 0.8 | Thin |
| Ribbon alpha | 0.45 | Semi-transparent CI band |
| Mean line width | 2.2 | Bold |
| Mean point size | 2.8 | Visible at publication size |
| Reference line | y = 0 | Dashed gray · `alpha = 0.6` |
| Y-axis limits | `c(−155, 185)` | Accommodates extreme responders and progressors |
| Timepoints (x-axis) | 0, 21, 42, 63, 84, 126, 168 | Days |
| Facet layout | `nrow = 1` | All four BOR groups side by side |
| Figure dimensions | 17 × 7 in | Wide landscape for four facets |
| DPI | 180 | Publication quality |

---

## Limitations

- **Missing timepoints:** Real ctDNA data has irregular follow-up — patients who progress or die early have fewer timepoints. The mean and CI at later timepoints will be biased toward patients who survive long enough to be assessed (informative censoring of the biomarker). Linear mixed-effects models are more appropriate for formal longitudinal analysis.
- **Single-arm analysis:** This figure shows the treatment arm only. A parallel control-arm panel would be needed to assess whether ctDNA decline is treatment-specific or reflects natural disease course.
- **No landmark analysis:** The figure does not formally test whether early ctDNA response (e.g., >50% decline at Day 42) predicts OS or PFS — this requires a landmark survival analysis.
- **Noise level:** The `N(0, 20)` noise added to synthetic trajectories may underestimate real-world ctDNA variability due to assay noise, clonal evolution, and inter-tumor heterogeneity.

---

## Requirements

```r
library(ggplot2)   # >= 3.4
library(dplyr)     # >= 1.1
```

---

## Files

| File | Description |
|---|---|
| `13_ctdna.R` | Self-contained script · synthetic data generated internally |
| `plots/13_ctdna.png` | Output figure · 17 × 7 in · 180 DPI |

---

## References

Tie J, et al. Circulating tumor DNA analysis detects minimal residual disease and predicts recurrence in patients with stage II colon cancer. *Sci Transl Med.* 2016;8(346):346ra92.

Moding EJ, et al. Circulating tumor DNA dynamics predict benefit from consolidation immunotherapy in locally advanced non-small cell lung cancer. *Nat Cancer.* 2020;1(2):176–183.

Gandara DR, et al. Blood-based tumor mutational burden as a predictor of clinical benefit in non-small-cell lung cancer patients treated with atezolizumab. *Nat Med.* 2018;24(9):1441–1448.

Merker JD, et al. Circulating Tumor DNA Analysis in Patients With Cancer: American Society of Clinical Oncology and College of American Pathologists Joint Review. *J Clin Oncol.* 2018;36(16):1631–1641.
