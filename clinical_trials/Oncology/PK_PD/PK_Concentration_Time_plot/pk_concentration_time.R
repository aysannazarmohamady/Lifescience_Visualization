# pk_concentration_time.R
# PK Concentration-Time Curve — Cycle 1 Day 1, all dose levels (dose-proportionality display)
# Standards applied: linear + semi-log dual panel, geometric mean +/- 95% CI, nominal time
# basis (stated), >50% BLQ timepoints excluded from the mean profile (standard PK reporting
# convention), colorblind-safe dose palette, N per dose group shown, PK-evaluable population.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
  library(scales)
})

SCRIPT_DIR <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
if (length(SCRIPT_DIR) == 0 || SCRIPT_DIR == "") SCRIPT_DIR <- getwd()
source(file.path(SCRIPT_DIR, "..", "pkpd_common.R"))

DATA_DIR   <- file.path(SCRIPT_DIR, "..", "Data", "V1")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "Out")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

adpk <- read.csv(file.path(DATA_DIR, "ADPK.csv"), stringsAsFactors = FALSE)

conc <- adpk |>
  filter(PARAMCD == "CONC", AVISIT == "CYCLE 1 DAY 1", ANL01FL == "Y")

doses      <- sort(unique(conc$DOSE))
timepoints <- sort(unique(conc$NOMTPT))
lloq       <- min(conc$AVAL[conc$AVAL > 0], na.rm = TRUE)

# Build per-dose, per-timepoint geometric mean summary, applying the >50% BLQ exclusion rule
summ <- list()
excluded <- character(0)
for (d in doses) {
  sub <- conc |> filter(DOSE == d)
  n_subj <- n_distinct(sub$USUBJID)
  for (t in timepoints) {
    vals <- sub$AVAL[sub$NOMTPT == t]
    pct_blq <- 100 * mean(vals <= 0)
    if (t > 0 && pct_blq > 50) {
      summ[[length(summ) + 1]] <- data.frame(DOSE = d, NOMTPT = t, GM = NA, LO = NA, HI = NA,
                                              N = length(vals), n_subj = n_subj)
      excluded <- c(excluded, sprintf("%d mg @ %gh (%.0f%% BLQ)", d, t, pct_blq))
      next
    }
    if (t == 0) {
      stat <- c(gm = mean(vals), lo = NA_real_, hi = NA_real_)
    } else {
      stat <- gm_ci(vals)
    }
    summ[[length(summ) + 1]] <- data.frame(DOSE = d, NOMTPT = t, GM = stat[["gm"]],
                                            LO = stat[["lo"]], HI = stat[["hi"]],
                                            N = length(vals), n_subj = n_subj)
  }
}
summ <- bind_rows(summ) |> filter(!is.na(GM))
summ$DOSE_f <- factor(summ$DOSE, levels = doses)

make_panel <- function(logscale) {
  d <- summ
  if (logscale) d <- d |> mutate(GM = ifelse(GM <= 0, lloq / 2, GM))

  p <- ggplot(d, aes(x = NOMTPT, y = GM, color = DOSE_f)) +
    geom_line(linewidth = 0.9) +
    geom_point(size = 1.8) +
    geom_errorbar(aes(ymin = LO, ymax = HI), width = 0.6, linewidth = 0.5, na.rm = TRUE) +
    scale_color_manual(values = DOSE_COLORS, name = "Dose (C1D1)",
                       labels = function(x) paste0(x, " mg")) +
    scale_x_continuous(breaks = timepoints) +
    labs(x = "Nominal time post-dose (h)",
         y = paste0("Vizatinib plasma concentration (ng/mL)", if (logscale) " \u2014 log scale" else ""),
         title = if (logscale) "Semi-log scale" else "Linear scale") +
    theme_pkpd()

  if (logscale) p <- p + scale_y_log10()
  p
}

p_lin <- make_panel(FALSE)
p_log <- make_panel(TRUE) + theme(legend.position = "none")

excl_txt <- if (length(excluded)) paste(excluded, collapse = "; ") else "none"
cap <- paste0(
  "Geometric mean \u00b1 95% CI; nominal sampling times plotted. Timepoints with >50% of concentrations below the ",
  "limit of quantification (BLQ) are not shown, per standard PK reporting convention: ", excl_txt, ". ",
  sprintf("LLOQ = %.3f ng/mL. ", lloq),
  "Dose-escalation cohorts (100/200/400 mg, n=5 each) and RP2D expansion cohort (300 mg, n=47) combined at Cycle 1 Day 1, PK-evaluable population."
)

combined <- (p_lin | p_log) +
  plot_annotation(
    title = "PK Concentration-Time Curve \u2014 Vizatinib (ONCVIZ-001), Cycle 1 Day 1, PK-evaluable population",
    caption = paste(strwrap(cap, width = 160), collapse = "\n"),
    theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5))
  )

ggsave(file.path(OUTPUT_DIR, "pk_concentration_time.png"), combined,
       width = 12, height = 6, dpi = 300)
message("Saved pk_concentration_time.png to ", OUTPUT_DIR)
