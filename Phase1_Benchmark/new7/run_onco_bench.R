# Render clean + OP01-OP04 flawed OncoPrint images. Run from the repository root:
#   Rscript Phase1_Benchmark/new7/run_onco_bench.R
# DRAFT NOTE: drafted with LLM assistance.
source("clinical_trials/Oncology/Biomarker_Genomics/OncoPrint/oncoprint_bench.R")
OUT <- "./Outputs/onco_bench"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) read.csv(file.path("./Data/V1", f), stringsAsFactors = FALSE)
prep <- prepare_onco_bench(rd("ADMUT.csv"), rd("ADSL.csv"))
render_onco_bench(prep, file.path(OUT, "OP_clean.png"))
render_onco_bench(prep, file.path(OUT, "OP_OP01.png"), freq_denom = length(unique(rd("ADMUT.csv")$USUBJID)))
render_onco_bench(prep, file.path(OUT, "OP_OP02.png"), legend_swap = TRUE)
render_onco_bench(prep, file.path(OUT, "OP_OP03.png"), shuffle = TRUE)
render_onco_bench(prep, file.path(OUT, "OP_OP04.png"), topbar_off = TRUE)
