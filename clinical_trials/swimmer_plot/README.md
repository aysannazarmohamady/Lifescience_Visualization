# Swimmer Plot

A swimmer plot is a patient-level visualization used in oncology clinical trials to display treatment duration, response assessments, and study endpoints for each individual patient on a single chart.

## When to use

- Phase I/II oncology trials
- Displaying individual patient timelines (treatment start → end → follow-up)
- Showing response assessments over time (CR, PR, SD, PD)
- Communicating study endpoints (death, withdrawal, still on treatment)

## What this plot shows

| Element | Description |
|---|---|
| Blue bar | Duration on treatment |
| Grey bar | Off treatment, still on study follow-up |
| `▶` arrow | Patient still on treatment at data cutoff |
| `✕` black | Death or early termination |
| `✕` red | Withdrawal of consent / lost to follow-up |
| `●` green | Complete Response (CR) |
| `■` blue | Partial Response (PR) |
| `▲` yellow | Stable Disease (SD) |
| `◆` red | Progressive Disease (PD) |

## Files

| File | Description |
|---|---|
| `swimmer_plot.py` | Study-agnostic Python script |
| `swimmer_plot_colab.ipynb` | Google Colab notebook version |
| `data/sample_data.xlsx` | Synthetic example data |
| `output/swimmer_plot_example.png` | Example output image |
