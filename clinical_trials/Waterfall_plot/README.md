# Waterfall Plot

A waterfall plot is a patient-level visualization used in oncology clinical trials to display the best percentage change from baseline in target lesion size for each individual patient, ordered from greatest reduction to greatest increase.

## When to use

- Phase I/II oncology trials
- Displaying best overall tumor response per patient
- Communicating RECIST-based response thresholds (PR, PD)
- Comparing response distribution across cohorts or dose levels

## What this plot shows

| Element | Description |
|---|---|
| Bar below 0% | Tumor reduction from baseline |
| Bar above 0% | Tumor increase from baseline |
| `--` at -30% | PR threshold (≥30% reduction) |
| `--` at +20% | PD threshold (≥20% increase) |
| Dark blue bar | Complete Response (CR) |
| Light blue bar | Partial Response (PR) |
| Yellow bar | Stable Disease (SD) |
| Red bar | Progressive Disease (PD) |
| Grey bar | Not Evaluable |

## Evaluability criteria

A patient is **Evaluable** if they have:
- A valid baseline sum of target lesion diameters (`> 0`)
- At least one post-baseline tumor measurement

Patients without baseline or post-baseline measurements are classified as **Not Evaluable** and displayed as a grey bar at 0%.

## Files

| File | Description |
|---|---|
| `waterfall_plot.py` | Study-agnostic Python script |
| `data/ADSL.csv` | Synthetic ADaM subject-level dataset (N=50) |
| `data/ADRS.csv` | Synthetic ADaM response dataset |
| `data/ADTR.csv` | Synthetic ADaM tumor measurement dataset |
| `output/waterfall_evaluable.png` | Evaluable patients output |
| `output/waterfall_all.png` | All patients output |
| `output/waterfall_cohort1.png` | Cohort 1 output |
| `output/waterfall_cohort2.png` | Cohort 2 output |
| `output/waterfall_cohort3.png` | Cohort 3 output |
| `output/waterfall_cohort4.png` | Cohort 4 output |

## Requirements

```
pandas
numpy
matplotlib
```

## References

- Seymour L, et al. (2017). iRECIST: guidelines for response criteria for use in trials testing immunotherapeutics. *Lancet Oncology*.
- JNCI (2016). Current and Evolving Methods to Visualize Biological Data in Cancer Research.
