# Lifescience Visualization
 
Data tells the story. This repository makes it visible.
 
---
 
## Why this exists
 
Clinical trials produce some of the most consequential data in medicine — yet the standard for how that data is communicated has lagged far behind its complexity. Regulatory submissions, scientific publications, and clinical team reviews still rely heavily on dense summary tables that require significant effort to interpret and offer little intuition about what is actually happening at the patient level.
 
Visualization changes that. A well-designed swimmer plot communicates treatment duration, response, and outcome for every patient in a trial simultaneously — something no table can do. A waterfall plot makes the distribution of tumor response immediately legible. A Kaplan-Meier curve has become one of the most recognized figures in oncology precisely because it translates survival probability into something physicians and regulators can reason about directly.
 
This repository is a systematic, production-ready implementation of the full spectrum of visualizations used in clinical trial research — from individual patient timelines to population-level pharmacokinetics, from adverse event profiling to subgroup forest plots. Each implementation is built on synthetic ADaM-standard datasets that mirror the structure of real study data, making them directly adaptable to actual trials.
 
The goal is not a gallery. It is a rigorous, reproducible reference — grounded in published analytical standards, implemented cleanly in Python, and structured so that each visualization can be understood, validated, and deployed independently.
 
---

## Implemented

| Plot | Description |
|---|---|
| [Swimmer plot](clinical_trials/swimmer_plot/) | Individual patient treatment timelines, responses, and endpoints |
| [Waterfall plot](clinical_trials/Waterfall_plot/) | Best % change from baseline in target lesion size |
| [Spider plot](clinical_trials/Spider_plot/) | Tumor burden trajectories over time |
| [Response heatmap](clinical_trials/Heatmap_plot/) | ORR by cohort and tumor type |
| Kaplan-Meier | Survival and time-to-event analysis |
| Forest plot | Subgroup treatment effect estimates |

---
 
## Scope
 
Clinical trial data spans five analytical domains. This repository covers all of them systematically. Each domain reflects a distinct set of scientific questions that visualization is uniquely positioned to answer.alization lives in its own folder with a README, a standalone Python script, synthetic data, and example outputs — ready to run as-is or adapt to real study data.

### Response & Efficacy
 
The central question of any oncology trial is whether the treatment works — and for whom. Response visualizations operate at the level of individual patients, making heterogeneity visible and enabling the kind of nuanced interpretation that aggregate statistics obscure.
 
Swimmer plots show each patient's full treatment timeline alongside response assessments and study endpoints, making it possible to see patterns — late responders, durable remissions, early discontinuations — that summary tables cannot capture. Waterfall plots rank patients by their best tumor response and communicate the distribution of efficacy at a glance. Spider plots add the temporal dimension, showing not just how much tumors changed but when and how consistently.
 
*Swimmer plot · Waterfall plot · Spider plot · Response heatmap · Best overall response bar chart · Dose-response curve · Sanctuary site plot · Time-to-response vs. duration-of-response · Concordance plot · RECIST transition plot · Target lesion trajectory plot*
 
### Survival & Time-to-Event
 
Most Phase III oncology trials are powered to detect a survival endpoint — overall survival or progression-free survival. The Kaplan-Meier curve has been the standard representation of survival data since 1958, and it remains indispensable precisely because it conveys both the magnitude and the uncertainty of the treatment effect over time.
 
Beyond the standard curve, this domain includes the diagnostics needed to validate survival models — proportional hazards testing via log-log plots and Schoenfeld residuals — as well as newer methods like RMST that are increasingly required when the proportional hazards assumption does not hold, and composite endpoint visualizations like the win ratio that have gained traction in regulatory settings.
 
*Kaplan-Meier curve · Competing risk plot · Landmark analysis plot · Hazard ratio forest plot · Log-log survival plot · Schoenfeld residual plot · RMST plot · Cumulative hazard plot · Net benefit curve · Win ratio plot*
 
### Safety & Tolerability
 
Safety is not a secondary concern in clinical trials — in Phase I, it is the primary endpoint. Adverse event data is high-dimensional, hierarchically coded via MedDRA, and graded by severity, making it one of the most analytically demanding domains to visualize well.
 
The eDISH plot was developed specifically to detect drug-induced liver injury by tracking ALT and bilirubin simultaneously — a pattern that emerges from the combination of values that neither alone would reveal. The tendril plot maps the temporal evolution of adverse events across the trial timeline, making it possible to see whether toxicities are early and resolving or late and accumulating. Butterfly plots allow direct comparison of AE profiles between treatment arms without requiring the reader to mentally subtract one table from another.
 
*AE dot plot · Butterfly plot · Circular / radar plot · eDISH plot · Shift plot · Lab longitudinal plot · Volcano plot · Tendril plot · AE onset time plot · QTc dispersion plot · Lab outlier boxplot · Hy's Law quadrant chart · Prevalence-over-time safety plot*
 
### Patient-Level & Subgroup Data
 
Individual patient data and subgroup analyses serve different but complementary purposes. Patient-level visualizations — dose modification timelines, concomitant medication overlays, missing data pattern maps — are essential for clinical operations and data quality review. Subgroup visualizations are central to regulatory submissions and scientific communication.
 
The forest plot is one of the four most common figure types in published randomized controlled trials, alongside Kaplan-Meier curves, flow diagrams, and repeated measures graphs. It distills treatment effect estimates and confidence intervals across demographic and clinical subgroups into a form that makes consistency — or heterogeneity — immediately apparent. The patient disposition Sankey diagram makes trial flow legible in a way that the standard CONSORT flowchart cannot match when dropout reasons are complex or numerous.
 
*Slide plot · Forest plot · Oncoprint · Patient profile dashboard · Subgroup interaction plot · Spaghetti plot · Patient disposition Sankey diagram · Dose intensity step plot · Concomitant medication timeline · Missing data pattern plot*
 
### PK/PD & Biomarkers
 
Pharmacokinetic and pharmacodynamic data underpin dose selection, exposure-response relationships, and biomarker-driven patient stratification. Population PK modeling and visual predictive checks are standard components of NDA and BLA submissions. Exposure-response plots are increasingly required by regulators as evidence that the selected dose is justified not just by safety but by efficacy.
 
The ROC curve and biomarker kinetic plots connect molecular data to clinical outcomes, supporting the development of companion diagnostics and patient enrichment strategies that have become central to modern oncology drug development.
 
*Concentration-time curve · Exposure-response plot · Visual Predictive Check (VPC) · Hysteresis loop plot · Biomarker kinetic plot · ROC curve · IVIVC plot · Population PK parameter distribution · Target engagement plot*
 
---
 
## Synthetic data
 
All visualizations are demonstrated on a shared synthetic ADaM dataset simulating a Phase I oncology trial — 50 patients, 4 dose cohorts, multiple tumor types. The dataset follows ADaM structure and naming conventions used in regulatory submissions.
 
| File | Description |
|---|---|
| `ADSL.csv` | Subject-level: demographics, treatment dates, study status, disposition |
| `ADRS.csv` | Response: CR, PR, SD, PD per assessment visit with dates |
| `ADTR.csv` | Tumor measurements: sum of target lesion diameters over time |
 
---
 
## Standards and references
 
The visualizations in this repository are grounded in published methodology:
 
- Chia PL, et al. (2016). Current and Evolving Methods to Visualize Biological Data in Cancer Research. *JNCI*, 108(8).
- Seymour L, et al. (2017). iRECIST: guidelines for response criteria for use in trials testing immunotherapeutics. *Lancet Oncology*, 18(3).
- FDA (2019). Exposure-Response Relationships — Study Design, Data Analysis, and Regulatory Applications.
- ICH E9 (R1) (2019). Addendum on Estimands and Sensitivity Analysis in Clinical Trials.
---
 
## Requirements
 
```
pandas
numpy
matplotlib
```
 
---
 
## Structure
 
Each visualization lives in its own folder:
 
```
plot_name/
├── README.md                  ← analytical context, design decisions, references
├── plot_name.py               ← study-agnostic Python script
├── plot_name_colab.ipynb      ← interactive Colab notebook
├── data/                      ← synthetic ADaM input files
└── output/                    ← example outputs at publication resolution
```
 
