# Render the clean variant-2 waterfall and one flawed image per taxonomy error
# (WF01-WF04). Run from the repository root:
#   Rscript Phase1_Benchmark/variant2/run_waterfall_v2.R
# DRAFT NOTE: drafted with LLM assistance; see waterfall_plot_v2.R.

source("clinical_trials/Oncology/Response_Assessment/Waterfall_plot/waterfall_plot_v2.R")

DATA_DIR <- "./Data/V1"
OUT <- "./Outputs/waterfall_v2"
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

adsl <- read.csv(file.path(DATA_DIR, "ADSL.csv"), stringsAsFactors = FALSE)
adtr <- read.csv(file.path(DATA_DIR, "ADTR.csv"), stringsAsFactors = FALSE)
adrs <- read.csv(file.path(DATA_DIR, "ADRS.csv"), stringsAsFactors = FALSE)
df <- prepare_waterfall_v2(adsl, adtr, adrs)

# Data-level injectors work on the plotted category of one patient.
# WF01: a treatment-arm PR patient is plotted as SD.
wf01_v2 <- function(df) {
  target <- sort(df$USUBJID[df$resp == "PR"])[1]
  df$resp[df$USUBJID == target] <- "SD"
  attr(df, "target") <- target
  df
}
# WF04: a different PR patient is plotted as CR (best response contradicts the data).
wf04_v2 <- function(df, exclude = NULL) {
  target <- sort(setdiff(df$USUBJID[df$resp == "PR"], exclude))[2]
  df$resp[df$USUBJID == target] <- "CR"
  attr(df, "target") <- target
  df
}

render_waterfall_v2(df, file.path(OUT, "WF_v2_clean.png"))
d1 <- wf01_v2(df);  render_waterfall_v2(d1, file.path(OUT, "WF_v2_WF01.png"))
d4 <- wf04_v2(df);  render_waterfall_v2(d4, file.path(OUT, "WF_v2_WF04.png"))
render_waterfall_v2(df, file.path(OUT, "WF_v2_WF02.png"), pd_line = FALSE)
render_waterfall_v2(df, file.path(OUT, "WF_v2_WF03.png"), sorted = FALSE)
writeLines(c(paste("WF01 target:", attr(d1, "target")),
             paste("WF04 target:", attr(d4, "target"))),
           file.path(OUT, "injection_targets.txt"))
message("Rendered clean + WF01-WF04 to ", OUT)
