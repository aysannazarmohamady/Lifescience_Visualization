# =============================================================================
# ONCVIZ-001 -- Forest plot, variant 2 (benchmark second implementation)
#
# How this differs from forest_plot.R (variant 1):
#   * endpoint is PFS (variant 1: OS); subgroup set differs (adds stage, liver
#     metastases; ECOG 0 vs >=1; PD-L1 High vs not High);
#   * a printed "HR (95% CI)" text column and n per row; point size is constant
#     (no weights);
#   * CI bounds are truncated at the axis limits (no silent dropping of CIs);
#   * the four taxonomy errors are switches (null_x, log_axis, ci_ok, n_shift)
#     instead of edited copies of the script.
#
# DRAFT NOTE: drafted with LLM assistance; an independent author must rewrite
# or review it if the paper claims hand-written benchmark scripts.
# Usage (repository root): Rscript <this file>   (draws the clean figure)
# =============================================================================

suppressPackageStartupMessages({ library(ggplot2); library(dplyr); library(survival) })

FP_LIM <- c(0.1, 10)

fp_cox <- function(d) {
  if (length(unique(d$TRT)) < 2 || sum(d$EVENT) < 2) return(c(NA, NA, NA))
  fit <- tryCatch(coxph(Surv(AVAL, EVENT) ~ TRT, data = d), error = function(e) NULL)
  if (is.null(fit)) return(c(NA, NA, NA))
  ci <- exp(confint(fit)); c(exp(coef(fit)), ci[1], ci[2])
}

prepare_forest_v2 <- function(adtte) {
  d <- adtte |> filter(PARAMCD == "PFS") |>
    mutate(EVENT = 1L - CNSR, TRT = as.integer(ARM == "TREATMENT"),
           ECOG_GRP = ifelse(ECOG == 0, "ECOG 0", "ECOG >=1"),
           PDL1_GRP = ifelse(PDL1GRP == "HIGH", "PD-L1 High", "PD-L1 not High"),
           STAGE_GRP = ifelse(STAGE == "STAGE IV", "Stage IV", "Stage I-III"))
  sg <- list(
    "All patients" = d,
    "NSCLC" = filter(d, TUMORTYPE == "NSCLC"), "BRCA" = filter(d, TUMORTYPE == "BRCA"),
    "HCC" = filter(d, TUMORTYPE == "HCC"), "CRC" = filter(d, TUMORTYPE == "CRC"),
    "PDAC" = filter(d, TUMORTYPE == "PDAC"),
    "Age <65" = filter(d, AGEGR1 == "<65"), "Age >=65" = filter(d, AGEGR1 == ">=65"),
    "Male" = filter(d, SEX == "M"), "Female" = filter(d, SEX == "F"),
    "ECOG 0" = filter(d, ECOG_GRP == "ECOG 0"), "ECOG >=1" = filter(d, ECOG_GRP == "ECOG >=1"),
    "Stage I-III" = filter(d, STAGE_GRP == "Stage I-III"), "Stage IV" = filter(d, STAGE_GRP == "Stage IV"),
    "Liver mets: yes" = filter(d, LIVERMETS == "Y"), "Liver mets: no" = filter(d, LIVERMETS == "N"),
    "TMB-High" = filter(d, TMBHIGH == "Y"), "TMB-Low" = filter(d, TMBHIGH == "N"),
    "PD-L1 High" = filter(d, PDL1_GRP == "PD-L1 High"),
    "PD-L1 not High" = filter(d, PDL1_GRP == "PD-L1 not High"))
  bind_rows(lapply(names(sg), function(lab) {
    x <- sg[[lab]]; r <- fp_cox(x)
    data.frame(label = lab, hr = r[1], lo = r[2], hi = r[3],
               n = length(unique(x$USUBJID)), overall = lab == "All patients")
  })) |> filter(!is.na(hr)) |>
    mutate(lo_c = pmax(lo, FP_LIM[1]), hi_c = pmin(hi, FP_LIM[2]))
}

# Switches: ci_ok = FALSE moves the first subgroup point outside its CI (FP01);
# null_x = reference line value (FP02); log_axis = FALSE uses a linear axis
# (FP03); n_shift adds to the n shown for the first subgroup (FP04).
render_forest_v2 <- function(fdf, out_png, ci_ok = TRUE, null_x = 1,
                             log_axis = TRUE, n_shift = 0) {
  if (!ci_ok) {
    i <- which(!fdf$overall & fdf$hi_c <= 5)[1]
    fdf$hr[i] <- fdf$hi_c[i] * 1.6
  }
  fdf$n_shown <- fdf$n
  if (n_shift != 0) fdf$n_shown[which(!fdf$overall)[1]] <- fdf$n[which(!fdf$overall)[1]] + n_shift
  fdf$y <- nrow(fdf):1
  fdf$txt <- sprintf("%.2f (%.2f-%.2f)", fdf$hr, fdf$lo, fdf$hi)
  xs <- if (log_axis) scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10))
        else scale_x_continuous(breaks = c(0.1, 0.25, 0.5, 1, 2, 5, 10))
  p <- ggplot(fdf, aes(y = y)) +
    geom_vline(xintercept = null_x, linetype = "dotted", colour = "#666666") +
    geom_segment(aes(x = lo_c, xend = hi_c, yend = y), linewidth = 0.7, colour = "#444444") +
    geom_point(aes(x = hr, shape = overall, fill = overall), size = 3, colour = "#222222") +
    scale_shape_manual(values = c(`FALSE` = 21, `TRUE` = 23), guide = "none") +
    scale_fill_manual(values = c(`FALSE` = "#6BAED6", `TRUE` = "#08306B"), guide = "none") +
    xs + scale_y_continuous(breaks = fdf$y,
                            labels = sprintf("%s  (n=%d)", fdf$label, fdf$n_shown)) +
    geom_text(aes(x = FP_LIM[2] * 1.15, label = txt), hjust = 0, size = 3) +
    coord_cartesian(xlim = FP_LIM, clip = "off") +
    labs(x = "Hazard ratio, Treatment vs Control (95% CI)", y = NULL,
         title = "Progression-free survival by subgroup (Cox model)") +
    theme_classic(base_size = 11) +
    theme(plot.margin = margin(8, 110, 8, 8))
  ggsave(out_png, p, width = 11, height = 8, dpi = 150, device = "png")
  base <- sub("\\.png$", "", out_png)
  write.csv(fdf[, c("label", "hr", "lo_c", "hi_c", "n_shown")],
            paste0(base, "_table.csv"), row.names = FALSE)
  writeLines(c(paste("null_x", null_x), paste("scale", if (log_axis) "log10" else "linear")),
             paste0(base, "_meta.txt"))
  invisible(fdf)
}

if (sys.nframe() == 0) {
  out <- "./Outputs/forest_v2"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
  adtte <- read.csv("./Data/V1/ADTTE.csv", stringsAsFactors = FALSE)
  render_forest_v2(prepare_forest_v2(adtte), file.path(out, "FP_v2_clean.png"))
}
