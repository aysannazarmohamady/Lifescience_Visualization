# Render clean + TX01-TX04 flawed AE bar-chart images. Run from the repository root:
#   Rscript Phase1_Benchmark/new7/run_ae_bench.R
# DRAFT NOTE: drafted with LLM assistance.
source("clinical_trials/Oncology/Safety_Toxicity/AE_Bar_Chart_plot/ae_bar_chart_bench.R")
OUT <- "./Outputs/ae_bench"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) read.csv(file.path("./Data/V1", f), stringsAsFactors = FALSE)
adae <- rd("ADAE.csv"); adsl <- rd("ADSL.csv")
tab <- prepare_ae_bench(adae, adsl)
render_ae_bench(tab, file.path(OUT, "TX_clean.png"))
render_ae_bench(tab, file.path(OUT, "TX_TX01.png"), bad_colours = TRUE)
render_ae_bench(tab, file.path(OUT, "TX_TX02.png"), drop_g45 = TRUE)
render_ae_bench(tab, file.path(OUT, "TX_TX03.png"), denom = nrow(adsl))   # N=80 instead of 62
render_ae_bench(tab, file.path(OUT, "TX_TX04.png"), shuffle = TRUE)
