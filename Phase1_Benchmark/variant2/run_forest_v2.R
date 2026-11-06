# Render the clean variant-2 forest plot and one flawed image per error (FP01-FP04).
# Run from the repository root: Rscript Phase1_Benchmark/variant2/run_forest_v2.R
# DRAFT NOTE: drafted with LLM assistance.
source("clinical_trials/Oncology/Meta_Analysis/forest_plot/forest_plot_v2.R")
OUT <- "./Outputs/forest_v2"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
fdf <- prepare_forest_v2(read.csv("./Data/V1/ADTTE.csv", stringsAsFactors = FALSE))
render_forest_v2(fdf, file.path(OUT, "FP_v2_clean.png"))
render_forest_v2(fdf, file.path(OUT, "FP_v2_FP01.png"), ci_ok = FALSE)
render_forest_v2(fdf, file.path(OUT, "FP_v2_FP02.png"), null_x = 2)
render_forest_v2(fdf, file.path(OUT, "FP_v2_FP03.png"), log_axis = FALSE)
render_forest_v2(fdf, file.path(OUT, "FP_v2_FP04.png"), n_shift = 25)
