# Render clean + PK01-PK04 flawed PK images. Run from the repository root:
#   Rscript Phase1_Benchmark/new7/run_pk_bench.R
# DRAFT NOTE: drafted with LLM assistance.
source("clinical_trials/Oncology/PK_PD/PK_Concentration_Time_plot/pk_concentration_time_bench.R")
OUT <- "./Outputs/pk_bench"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
s <- prepare_pk_bench(read.csv("./Data/V1/ADPK.csv", stringsAsFactors = FALSE))
render_pk_bench(s, file.path(OUT, "PK_clean.png"))
render_pk_bench(s, file.path(OUT, "PK_PK01.png"), linear_y = TRUE)
render_pk_bench(s, file.path(OUT, "PK_PK02.png"), arith = TRUE)
render_pk_bench(s, file.path(OUT, "PK_PK03.png"), index_x = TRUE)
render_pk_bench(s, file.path(OUT, "PK_PK04.png"), swap_labels = TRUE)
