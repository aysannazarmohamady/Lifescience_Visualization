# Variant 2 of the Phase 1 benchmark (second implementation per plot type)

Each of the five plot types has a second, independently structured R script, a
runner that renders one clean image plus one flawed image per taxonomy error
(4 per plot type), and a Python verifier. The verifier requires that the clean
image passes every check and that each flawed image fails exactly the check of
its own error and differs from the clean image (md5).

Run from the repository root (data in `Data/V1/`, unzipped from `oncoviz_data_v1.zip`):

| Plot type | Script | Runner | Verifier |
|---|---|---|---|
| Waterfall | `clinical_trials/Oncology/Response_Assessment/Waterfall_plot/waterfall_plot_v2.R` | `run_waterfall_v2.R` | `verify_waterfall_v2.py` |
| Swimmer | `clinical_trials/Oncology/Response_Assessment/swimmer_plot/swimmer_plot_v2.R` | `run_swimmer_v2.R` | `verify_swimmer_v2.py` |
| Forest | `clinical_trials/Oncology/Meta_Analysis/forest_plot/forest_plot_v2.R` | `run_forest_v2.R` | `verify_forest_v2.py` |
| CONSORT | `clinical_trials/Oncology/Trial_Design_Patient_Flow/CONSORT_Diagram_plot/consort_diagram_plot_v2.R` | `run_consort_v2.R` | `verify_consort_v2.py` |
| Kaplan-Meier | `clinical_trials/Oncology/Survival_TimeToEvent/SurvivalCurves_plot/survivalcurves_plot_v2.R` | `run_km_v2.R` | `verify_km_v2.py` |

```
Rscript Phase1_Benchmark/variant2/run_<plot>_v2.R
python3 -I Phase1_Benchmark/variant2/verify_<plot>_v2.py     # prints ALL OK
```

Findings recorded while building variant 2
- Waterfall v1: `waterfall_plot.R` crashed (missing `tumor` column); fixed with one line.
  The clean v1 data also fail `check_waterfall_color_matches_pchg` (37/62 ADRS colours
  differ from the category implied by best % change); v2 derives the plotted category
  from % change so the clean reference is consistent.
- Repository checkers that are too weak to detect their own error and were
  strengthened in the verifiers: `check_swimmer_marker_matches_assessment_date`,
  `check_swimmer_bar_length_matches_dates`, `check_consort_exclusion_reasons_present`
  (any newline passes), CN03 (relabelling one control patient keeps the arm total
  at 80; detectable only against the per-arm counts in ADRAND).
- Kaplan-Meier v1 (survminer) does not run with ggplot2 3.4.4; v2 does not use survminer.
- All variant-2 files were drafted with LLM assistance; an independent author must
  rewrite or review them if the paper claims hand-written benchmark scripts.
