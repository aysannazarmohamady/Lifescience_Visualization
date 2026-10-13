# Render clean + EN01-EN04 flawed enrollment images. Run from the repository root:
#   Rscript Phase1_Benchmark/new7/run_enrollment_bench.R
# DRAFT NOTE: drafted with LLM assistance.
source("clinical_trials/Oncology/Trial_Design_Patient_Flow/Enrollment_Over_Time_plot/enrollment_over_time_bench.R")
OUT <- "./Outputs/enrollment_bench"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
adrand <- read.csv("./Data/V1/ADRAND.csv", stringsAsFactors = FALSE)
adsl   <- read.csv("./Data/V1/ADSL.csv", stringsAsFactors = FALSE)

# EN03 data injector: the last 5 randomized patients get a first-dose date after closure.
en03_shift_after_closure <- function(adsl, adrand, n = 5, new_date = "2026-04-15") {
  last <- adrand$USUBJID[adrand$RANDFL == "Y"]
  last <- last[order(as.Date(adsl$TRTSDT[match(last, adsl$USUBJID)]), decreasing = TRUE)][seq_len(n)]
  adsl$TRTSDT[adsl$USUBJID %in% last] <- new_date
  adsl
}

prep <- prepare_enrollment_bench(adrand, adsl)
render_enrollment_bench(prep, file.path(OUT, "EN_clean.png"))
render_enrollment_bench(prep, file.path(OUT, "EN_EN01.png"), monotone = FALSE)
render_enrollment_bench(prep, file.path(OUT, "EN_EN02.png"), arm_label_shift = -3)
render_enrollment_bench(prepare_enrollment_bench(adrand, en03_shift_after_closure(adsl, adrand)),
                        file.path(OUT, "EN_EN03.png"))
render_enrollment_bench(prep, file.path(OUT, "EN_EN04.png"), target_line = FALSE)
