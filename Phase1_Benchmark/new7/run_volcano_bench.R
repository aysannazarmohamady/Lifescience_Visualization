# Render clean + VP01-VP04 flawed volcano images. Run from the repository root:
#   Rscript Phase1_Benchmark/new7/run_volcano_bench.R
# DRAFT NOTE: drafted with LLM assistance.
source("clinical_trials/Oncology/Biomarker_Genomics/Volcano_plot/volcano_plot_bench.R")
OUT <- "./Outputs/volcano_bench"; dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
adbm <- read.csv("./Data/V1/ADBM.csv", stringsAsFactors = FALSE)
r <- prepare_volcano_bench(adbm)
cat("features:", nrow(r), " significant:", sum(r$sig), "\n")
render_volcano_bench(r, file.path(OUT, "VP_clean.png"))
render_volcano_bench(r, file.path(OUT, "VP_VP01.png"), thr_alpha = 0.25)   # wrong threshold line
render_volcano_bench(r, file.path(OUT, "VP_VP02.png"), ln_p = TRUE)        # -ln(p) on a "log10" axis
render_volcano_bench(r, file.path(OUT, "VP_VP03.png"), fc_raw = TRUE)      # ratio instead of log2
render_volcano_bench(r, file.path(OUT, "VP_VP04.png"), wrong_highlight = 3)# non-significant highlighted
