# Render the clean variant-2 CONSORT diagram and one flawed image per error (CN01-CN04).
# Run from the repository root: Rscript Phase1_Benchmark/variant2/run_consort_v2.R
# DRAFT NOTE: drafted with LLM assistance.
source("clinical_trials/Oncology/Trial_Design_Patient_Flow/CONSORT_Diagram_plot/consort_diagram_plot_v2.R")
source("Phase1_Benchmark/injectors/injectors_consort.R")   # cn01 / cn03 data injectors
OUT <- "./Outputs/consort_v2"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
adrand <- read.csv("./Data/V1/ADRAND.csv", stringsAsFactors = FALSE)
adsl   <- read.csv("./Data/V1/ADSL.csv",   stringsAsFactors = FALSE)
render_consort_v2(adrand, adsl, file.path(OUT, "CN_v2_clean.png"))
render_consort_v2(cn01_add_unaccounted_patient(adrand), adsl, file.path(OUT, "CN_v2_CN01.png"))
render_consort_v2(adrand, adsl, file.path(OUT, "CN_v2_CN02.png"), excl_reasons = FALSE)
render_consort_v2(adrand, cn03_relabel_arm(adsl), file.path(OUT, "CN_v2_CN03.png"))
render_consort_v2(adrand, adsl, file.path(OUT, "CN_v2_CN04.png"), disc_detail = FALSE)
