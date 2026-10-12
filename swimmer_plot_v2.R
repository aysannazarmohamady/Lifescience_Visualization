# =============================================================================
# ONCVIZ-001 -- Swimmer plot, variant 2 (benchmark second implementation)
#
# How this differs from swimmer_plot.R (variant 1):
#   * horizontal segments (geom_segment) instead of tiles, time axis in DAYS
#     instead of months;
#   * one figure with a facet per tumor type (variant 1: one panel for all
#     patients, plus one figure per tumor type);
#   * rows sorted by treatment duration (descending) within tumor type,
#     not by best overall response;
#   * ongoing treatment is a segment ending in an arrow head;
#   * the two code_patch errors (SW03, SW04) are switches (ongoing_arrow,
#     sorted) instead of edited copies of the script.
#
# DRAFT NOTE: drafted with LLM assistance; an independent author must rewrite
# or review it if the paper claims hand-written benchmark scripts.
# Usage (repository root): Rscript <this file>   (draws the clean figure)
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

SW_RESP_COLORS <- c(CR = "#1A9641", PR = "#2B6CB0", SD = "#E8A020", PD = "#C0392B")
SW_SHAPES      <- c(CR = 16, PR = 15, SD = 17, PD = 18)
SW_TUMORS      <- c("NSCLC", "BRCA", "HCC", "CRC", "PDAC")

# One row per treatment-arm patient (bar length in days) and one row per
# plotted response marker (day of assessment since treatment start).
prepare_swimmer_v2 <- function(adsl, adrs) {
  trt <- adsl[adsl$ARM == "TREATMENT", ] |>
    mutate(TRTSDT = as.Date(TRTSDT), TRTEDT = as.Date(TRTEDT),
           bar_days = as.numeric(TRTEDT - TRTSDT),
           ongoing  = EOSSTT == "ONGOING",
           tumor    = factor(TUMORTYPE, levels = SW_TUMORS)) |>
    select(USUBJID, tumor, TRTSDT, bar_days, ongoing)

  markers <- adrs |>
    filter(PARAMCD == "OVRLRESP", AVALC %in% names(SW_SHAPES), !is.na(ADT), ADT != "") |>
    mutate(ADT = as.Date(ADT)) |>
    inner_join(trt[, c("USUBJID", "TRTSDT")], by = "USUBJID") |>
    mutate(day = as.numeric(ADT - TRTSDT)) |>
    arrange(USUBJID, day) |>
    group_by(USUBJID, AVALC) |> slice(1) |> ungroup() |>
    select(USUBJID, resp = AVALC, day)
  list(patients = trt, markers = markers)
}

# Writes <png>, <base>_rows.csv (row order and bar length) and
# <base>_markers.csv, <base>_ongoing.txt so checkers can verify the figure.
render_swimmer_v2 <- function(prep, out_png, ongoing_arrow = TRUE, sorted = TRUE,
                              seed = 1) {
  pts <- prep$patients
  if (sorted) {
    pts <- pts |> arrange(tumor, desc(bar_days))
  } else {
    set.seed(seed)
    pts <- pts |> arrange(tumor) |> group_by(tumor) |>
      slice(sample(n())) |> ungroup()
  }
  pts <- pts |> group_by(tumor) |> mutate(row = n():1) |> ungroup()
  # keep file order = plotted top-to-bottom order
  pts <- pts |> arrange(tumor, desc(row))

  mk <- prep$markers |> inner_join(pts[, c("USUBJID", "tumor", "row")], by = "USUBJID")
  arrows_df <- if (ongoing_arrow) pts[pts$ongoing, ] else pts[0, ]

  p <- ggplot() +
    geom_segment(data = pts, aes(x = 0, xend = bar_days, y = row, yend = row),
                 linewidth = 2.6, colour = "#34495E", lineend = "butt") +
    geom_segment(data = arrows_df,
                 aes(x = bar_days, xend = bar_days + 60, y = row, yend = row),
                 linewidth = 0.8, colour = "#34495E",
                 arrow = arrow(length = unit(0.12, "inches"), type = "closed")) +
    geom_point(data = mk, aes(x = day, y = row, shape = resp, colour = resp),
               size = 2.4) +
    scale_shape_manual(values = SW_SHAPES, name = "Response") +
    scale_colour_manual(values = SW_RESP_COLORS, name = "Response") +
    scale_x_continuous(name = "Days since treatment start", breaks = seq(0, 1500, 100)) +
    facet_grid(tumor ~ ., scales = "free_y", space = "free_y", switch = "y") +
    labs(title = sprintf("Swimmer plot (treatment arm, N=%d); arrow = treatment ongoing",
                         nrow(pts)),
         y = NULL) +
    theme_minimal(base_size = 10) +
    theme(axis.text.y = element_blank(), panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(), legend.position = "top",
          strip.placement = "outside", strip.text.y.left = element_text(angle = 0, face = "bold"))
  ggsave(out_png, p, width = 12, height = 12, dpi = 120, device = "png")

  base <- sub("\\.png$", "", out_png)
  write.csv(pts[, c("USUBJID", "tumor", "row", "bar_days")],
            paste0(base, "_rows.csv"), row.names = FALSE)
  write.csv(mk[, c("USUBJID", "resp", "day")], paste0(base, "_markers.csv"), row.names = FALSE)
  writeLines(arrows_df$USUBJID, paste0(base, "_ongoing.txt"))
  invisible(pts)
}

if (sys.nframe() == 0) {
  out <- "./Outputs/swimmer_v2"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
  adsl <- read.csv("./Data/V1/ADSL.csv", stringsAsFactors = FALSE)
  adrs <- read.csv("./Data/V1/ADRS.csv", stringsAsFactors = FALSE)
  render_swimmer_v2(prepare_swimmer_v2(adsl, adrs), file.path(out, "SW_v2_clean.png"))
}
