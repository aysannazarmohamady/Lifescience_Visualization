# Render the clean variant-2 swimmer plot and one flawed image per error (SW01-SW04).
# Run from the repository root: Rscript Phase1_Benchmark/variant2/run_swimmer_v2.R
# DRAFT NOTE: drafted with LLM assistance.
source("clinical_trials/Oncology/Response_Assessment/swimmer_plot/swimmer_plot_v2.R")
source("Phase1_Benchmark/injectors/injectors_swimmer.R")   # sw01 / sw02 data injectors

OUT <- "./Outputs/swimmer_v2"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
adsl <- read.csv("./Data/V1/ADSL.csv", stringsAsFactors = FALSE)
adrs <- read.csv("./Data/V1/ADRS.csv", stringsAsFactors = FALSE)

render_swimmer_v2(prepare_swimmer_v2(adsl, adrs), file.path(OUT, "SW_v2_clean.png"))
render_swimmer_v2(prepare_swimmer_v2(adsl, sw01_shift_response_date(adrs, adsl)),
                  file.path(OUT, "SW_v2_SW01.png"))
render_swimmer_v2(prepare_swimmer_v2(sw02_shift_trtedt(adsl), adrs),
                  file.path(OUT, "SW_v2_SW02.png"))
render_swimmer_v2(prepare_swimmer_v2(adsl, adrs), file.path(OUT, "SW_v2_SW03.png"),
                  ongoing_arrow = FALSE)
render_swimmer_v2(prepare_swimmer_v2(adsl, adrs), file.path(OUT, "SW_v2_SW04.png"),
                  sorted = FALSE)
