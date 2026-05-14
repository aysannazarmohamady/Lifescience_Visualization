# Lifescience Visualization

Data tells the story. This repository makes it visible.

Clinical trials generate complex, multi-dimensional data about patients, tumors, responses, and outcomes. Turning that data into clear, reproducible, publication-ready visualizations is both a scientific and a communication challenge.

This repository is a growing collection of Python-based visualizations built for oncology and life science research — each one grounded in real analytical standards, implemented cleanly, and demonstrated on synthetic ADaM datasets that mirror real study data.

---

## Visualizations

### Clinical trials

| Plot | Description |
|---|---|
| [Swimmer plot](clinical_trials/swimmer_plot/) | Individual patient treatment timelines, responses, and endpoints |
| [Waterfall plot](clinical_trials/Waterfall_plot/) | Best % change from baseline in target lesion size |
| [Spider plot](clinical_trials/Spider_plot/) | Tumor burden trajectories over time |
| Kaplan-Meier | Survival and time-to-event analysis |
| Forest plot | Subgroup treatment effect estimates |

### Coming soon

Safety · PK/PD · Biomarker · Genomics · Single cell

---

## Synthetic data

All visualizations run on a shared synthetic ADaM dataset — a simulated Phase I oncology trial, 50 patients, 4 dose cohorts.

| File | Description |
|---|---|
| `ADSL.csv` | Subject-level: demographics, treatment dates, status |
| `ADRS.csv` | Response: CR, PR, SD, PD per assessment visit |
| `ADTR.csv` | Tumor measurements: sum of diameters over time |

---

## Requirements

```
pandas
numpy
matplotlib
```

---

## Structure

Each visualization lives in its own folder with a README, a standalone Python script, a Colab notebook, synthetic data, and example outputs — ready to run as-is or adapt to real study data.
