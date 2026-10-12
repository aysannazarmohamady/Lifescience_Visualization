# Render the clean variant-2 Kaplan-Meier plot and one flawed image per error (KM01-KM04).
# Run from the repository root: Rscript Phase1_Benchmark/variant2/run_km_v2.R
# DRAFT NOTE: drafted with LLM assistance.
source("clinical_trials/Oncology/Survival_TimeToEvent/SurvivalCurves_plot/survivalcurves_plot_v2.R")
source("Phase1_Benchmark/injectors/injectors_kaplan_meier.R")   # km02 data injector
OUT <- "./Outputs/km_v2"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
adtte <- read.csv("./Data/V1/ADTTE.csv", stringsAsFactors = FALSE)
render_km_v2(adtte, file.path(OUT, "KM_v2_clean.png"))
render_km_v2(adtte, file.path(OUT, "KM_v2_KM01.png"), corrupt = "KM01")
render_km_v2(km02_zero_censor_flag(adtte), file.path(OUT, "KM_v2_KM02.png"))
render_km_v2(adtte, file.path(OUT, "KM_v2_KM03.png"), corrupt = "KM03")
render_km_v2(adtte, file.path(OUT, "KM_v2_KM04.png"), corrupt = "KM04")
