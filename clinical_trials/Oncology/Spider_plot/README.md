# Spider Plot

A spider plot (also known as a spaghetti plot) is a patient-level visualization used in oncology clinical trials to display the change in target lesion size over time for each individual patient, allowing assessment of tumor response trajectories throughout treatment.

## When to use

- Phase I/II oncology trials
- Displaying individual tumor burden trajectories over time
- Visualizing durability of response alongside waterfall plot
- Identifying patients with delayed or sustained responses

## What this plot shows

| Element | Description |
|---|---|
| Each line | One patient's tumor burden change over time |
| Endpoint dot | Last available on-study measurement per patient |
| Y = 0% | Baseline (no change) |
| `--` at −30% | PR threshold (≥30% reduction) |
| `--` at +20% | PD threshold (≥20% increase) |
| Values clipped at +100% | Extreme progressors capped for display |
| Line colour (Panel A) | Best overall response (CR / PR / SD / PD / NE) |
| Line colour (Panel B) | Tumor type (NSCLC / CRC / HCC / PDAC / BRCA) |

## Output panels

Two figures are produced, each covering the **Treatment arm** only:

| Output file | Stratification | Colour by |
|---|---|---|
| `spider_plot_A_tumor_type.png` | One panel per tumor type (NSCLC · CRC · HCC · PDAC · BRCA) | Best overall response |
| `spider_plot_B_tmb.png` | TMB-High vs TMB-Low (two panels) | Tumor type |

Each panel subtitle reports **N**, number of CRs, number of PRs, and **ORR (%)**.

## Evaluability criteria

A patient is included if they have:

- A valid baseline sum of longest diameters (`AVAL > 0` at `AVISITN = 0`)
- At least **two** post-baseline tumor measurements (`AVISITN > 0`)

## Best overall response derivation

Overall response is taken from `ADRS` where `PARAMCD = "OVRLRESP"`. The best response is the record with the highest rank in the hierarchy **CR > PR > SD > PD > NE**. If no qualifying record exists the patient is assigned **NE**.

## Files

| File | Description |
|---|---|
| `spider_plot.R` | Study-agnostic R script |
| `Data/V1/ADSL.csv` | Synthetic ADaM subject-level dataset |
| `Data/V1/ADTR.csv` | Synthetic ADaM tumor measurement dataset |
| `Data/V1/ADRS.csv` | Synthetic ADaM response dataset |
| `Outputs/spider_plot_A_tumor_type.png` | Spider plot stratified by tumor type |
| `Outputs/spider_plot_B_tmb.png` | Spider plot stratified by TMB status |

## Key parameters

| Parameter | Value | Description |
|---|---|---|
| `DAY2MO` | 30.4375 | Days-to-months conversion factor |
| `PR_TH` | −30 | Partial response threshold (%) |
| `PD_TH` | +20 | Progressive disease threshold (%) |
| `CLIP_TOP` | +100 | Upper clip value for y-axis (%) |
| `TUMOR_ORDER` | NSCLC, CRC, HCC, PDAC, BRCA | Panel display order |

## Requirements

```r
library(dplyr)
library(ggplot2)
library(patchwork)
library(gridExtra)  # for grid.arrange
```

## References

- Eisenhauer EA, et al. (2009). New response evaluation criteria in solid tumours: Revised RECIST guideline (version 1.1). *European Journal of Cancer*.
- Seymour L, et al. (2017). iRECIST: guidelines for response criteria for use in trials testing immunotherapeutics. *Lancet Oncology*.
- JNCI (2016). Current and Evolving Methods to Visualize Biological Data in Cancer Research.
