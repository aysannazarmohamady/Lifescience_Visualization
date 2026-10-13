# Render clean + SP01-SP04 flawed spider images. Run from the repository root:
#   Rscript Phase1_Benchmark/new7/run_spider_bench.R
# DRAFT NOTE: drafted with LLM assistance.
source("clinical_trials/Oncology/Response_Assessment/Spider_plot/spider_plot_bench.R")
OUT <- "./Outputs/spider_bench"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
rd <- function(f) read.csv(file.path("./Data/V1", f), stringsAsFactors = FALSE)
adtr <- rd("ADTR.csv"); adrs <- rd("ADRS.csv"); adsl <- rd("ADSL.csv")

# SP02 data injector: 10 patients' whole trajectory is shifted up, so it no
# longer starts at 0% at baseline.
sp02_shift_baseline <- function(adtr, n = 10, shift = 12) {
  ids <- sort(unique(adtr$USUBJID))[seq_len(n)]
  adtr$PCHG[adtr$USUBJID %in% ids] <- adtr$PCHG[adtr$USUBJID %in% ids] + shift
  adtr
}

clean <- prepare_spider_bench(adtr, adrs, adsl)
render_spider_bench(clean, file.path(OUT, "SP_clean.png"))
render_spider_bench(clean, file.path(OUT, "SP_SP01.png"), ref_lines = FALSE)
render_spider_bench(prepare_spider_bench(sp02_shift_baseline(adtr), adrs, adsl), file.path(OUT, "SP_SP02.png"))
render_spider_bench(clean, file.path(OUT, "SP_SP03.png"), truncate = FALSE)
render_spider_bench(clean, file.path(OUT, "SP_SP04.png"), mixed_units = TRUE)
