# Render clean + CI01-CI04 flawed cumulative-incidence images. Run from the repo root:
#   Rscript Phase1_Benchmark/new7/run_cif_bench.R
# DRAFT NOTE: drafted with LLM assistance.
source("clinical_trials/Oncology/Survival_TimeToEvent/TimetoEvent_plot/cif_bench.R")
OUT <- "./Outputs/cif_bench"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
d <- prepare_cif_bench(read.csv("./Data/V1/ADTTE.csv", stringsAsFactors = FALSE))
render_cif_bench(d, file.path(OUT, "CI_clean.png"))
render_cif_bench(d, file.path(OUT, "CI_CI01.png"), one_minus_km = TRUE)
render_cif_bench(d, file.path(OUT, "CI_CI02.png"), wrong_event = TRUE)
render_cif_bench(d, file.path(OUT, "CI_CI03.png"), days_axis = TRUE)
render_cif_bench(d, file.path(OUT, "CI_CI04.png"), no_risk = TRUE)
