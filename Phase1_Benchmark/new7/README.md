# Benchmark for the seven remaining plot types (new7)

One clean image and four flawed images per plot type (28 errors, see `taxonomy_new7.csv`).
Each flawed image differs from the clean one by exactly one error. Errors are
arguments of the render function (or a small data injector in the runner); the
plot code itself sits next to the original plot as `*_bench.R`.

| Plot | Script | Runner / verifier |
|---|---|---|
| Enrollment over time | `Trial_Design_Patient_Flow/Enrollment_Over_Time_plot/enrollment_over_time_bench.R` | `run_enrollment_bench.R`, `verify_enrollment_bench.py` |
| Spider | `Response_Assessment/Spider_plot/spider_plot_bench.R` | `run_spider_bench.R`, `verify_spider_bench.py` |
| Volcano | `Biomarker_Genomics/Volcano_plot/volcano_plot_bench.R` | `run_volcano_bench.R`, `verify_volcano_bench.py` |
| AE bar chart | `Safety_Toxicity/AE_Bar_Chart_plot/ae_bar_chart_bench.R` | `run_ae_bench.R`, `verify_ae_bench.py` |
| OncoPrint | `Biomarker_Genomics/OncoPrint/oncoprint_bench.R` | `run_onco_bench.R`, `verify_onco_bench.py` |
| PK concentration-time | `PK_PD/PK_Concentration_Time_plot/pk_concentration_time_bench.R` | `run_pk_bench.R`, `verify_pk_bench.py` |
| Cumulative incidence | `Survival_TimeToEvent/TimetoEvent_plot/cif_bench.R` | `run_cif_bench.R`, `verify_cif_bench.py` |

First unzip the data once: `cd Data/V1 && unzip -o oncoviz_data_v1.zip`. Then run everything from the repository root: `bash Phase1_Benchmark/new7/run_all.sh`
(R with ggplot2, cmprsk; Python with pandas, scipy). Images go to `Outputs/` (git-ignored).
Each verifier recomputes ground truth from the clean ADaM data, reads the figure's
plotted values from side files written by the renderer, and asserts that the clean image
passes every check and each flawed image fails only its own check.

Notes: scripts are LLM-drafted and need independent review. The volcano plot is
exploratory/unadjusted because the synthetic biomarker data have no adjusted-significant
signal. In the cumulative-incidence data there is only one competing event (treatment arm).
