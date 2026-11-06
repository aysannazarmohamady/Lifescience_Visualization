# =============================================================================
# ONCVIZ-001 -- Waterfall plot, variant 2 (benchmark second implementation)
#
# How this differs from waterfall_plot.R (variant 1):
#   * one panel with all treatment-arm patients (variant 1 draws one panel per
#     tumor type plus a liver-metastasis strip);
#   * a tumor-type colour strip under the bars instead of the liver-mets strip;
#   * the plotted response category is derived from best % change with the
#     RECIST thresholds (CR <= -100, PR <= -30, PD >= 20, otherwise SD), so a
#     clean figure always satisfies the colour-vs-percent-change rule. The ADRS
#     overall response is kept in column resp_adrs for reference only (it also
#     reflects non-target and new lesions, so it can disagree with the bar);
#   * best % change is derived with dplyr joins, not row-wise helper functions;
#   * the two code_patch errors (WF02, WF03) are switches (pd_line, sorted)
#     rather than edited copies of the script.
#
# DRAFT NOTE: this file was drafted with LLM assistance. If the revised paper
# claims hand-written benchmark scripts, an independent author must rewrite or
# review it before release.
#
# Usage (from the repository root):
#   Rscript clinical_trials/Oncology/Response_Assessment/Waterfall_plot/waterfall_plot_v2.R
# =============================================================================

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

RESP_COLORS  <- c(CR = "#1A3A7C", PR = "#4A90C4", SD = "#E8A020",
                  PD = "#C0392B", NE = "#AAAAAA")
TUMOR_COLORS <- c(NSCLC = "#264653", BRCA = "#2A9D8F", HCC = "#E9C46A",
                  CRC = "#F4A261", PDAC = "#E76F51")
RESP_RANK    <- c(CR = 1, PR = 2, SD = 3, PD = 4, NE = 5)
PR_TH <- -30
PD_TH <- 20

# Treatment-arm table: one row per patient with best post-baseline % change in
# sum of lesion diameters (ADTR), the RECIST category implied by it (resp), and
# the ADRS overall response (resp_adrs, reference only).
recist_category <- function(pct, pr_th = PR_TH, pd_th = PD_TH) {
  ifelse(pct <= -100, "CR",
         ifelse(pct <= pr_th, "PR", ifelse(pct >= pd_th, "PD", "SD")))
}

prepare_waterfall_v2 <- function(adsl, adtr, adrs) {
  trt <- adsl[adsl$ARM == "TREATMENT", c("USUBJID", "TUMORTYPE", "DOSELEVEL")]

  best_resp <- adrs |>
    filter(PARAMCD == "OVRLRESP", AVALC %in% names(RESP_RANK)) |>
    mutate(rank = RESP_RANK[AVALC]) |>
    group_by(USUBJID) |>
    slice_min(rank, n = 1, with_ties = FALSE) |>
    ungroup() |>
    select(USUBJID, resp_adrs = AVALC)

  best_pct <- adtr |>
    filter(AVISITN > 0) |>
    mutate(PCHG = suppressWarnings(as.numeric(PCHG))) |>
    filter(!is.na(PCHG)) |>
    group_by(USUBJID) |>
    summarise(pct = min(PCHG), .groups = "drop")

  trt |>
    left_join(best_resp, by = "USUBJID") |>
    left_join(best_pct,  by = "USUBJID") |>
    mutate(resp_adrs = ifelse(is.na(resp_adrs), "NE", resp_adrs),
           pct   = ifelse(is.na(pct), 0, pct),
           resp  = recist_category(pct),
           short = sub("^[^-]+-[^-]+-", "", USUBJID))
}

# Render one waterfall. Writes <out_png> plus two small files describing what
# was drawn, so checkers can verify the figure without reading pixels:
#   <out_png minus .png>_order.csv   plotted bar order (x, USUBJID, pct, resp)
#   <out_png minus .png>_hlines.txt  y-values of every reference line drawn
render_waterfall_v2 <- function(df, out_png, pd_line = TRUE, sorted = TRUE,
                                seed = 1) {
  if (sorted) {
    df <- df |> arrange(desc(pct))
  } else {
    set.seed(seed)
    df <- df[sample(nrow(df)), ]
  }
  df <- df |> mutate(x = row_number())
  n  <- nrow(df)
  orr <- round(sum(df$resp %in% c("CR", "PR")) / n * 100)

  hlines <- c(0, PR_TH, if (pd_line) PD_TH)

  bars <- ggplot(df, aes(x = x, y = pct, fill = resp)) +
    geom_col(width = 0.8) +
    geom_hline(yintercept = 0, linewidth = 0.9, color = "black") +
    geom_hline(yintercept = PR_TH, linewidth = 0.6, color = "#555555",
               linetype = "dashed") +
    (if (pd_line) geom_hline(yintercept = PD_TH, linewidth = 0.6,
                             color = "#555555", linetype = "dashed")) +
    scale_fill_manual(values = RESP_COLORS, name = "Response category (from best % change)",
                      drop = FALSE) +
    scale_x_continuous(breaks = df$x, labels = df$short, expand = c(0.01, 0)) +
    scale_y_continuous(breaks = seq(-100, 100, 20),
                       labels = function(v) paste0(v, "%"),
                       limits = c(-105, 110)) +
    labs(title = sprintf("Best %% change from baseline in tumor size, treatment arm (N=%d, ORR %d%%)",
                         n, orr),
         x = NULL, y = "Best % change in SLD") +
    theme_minimal(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                     size = 6.5, family = "mono"),
          panel.grid.minor = element_blank(),
          panel.grid.major.x = element_blank(),
          legend.position = "top")

  strip <- ggplot(df, aes(x = x, y = 1, fill = TUMORTYPE)) +
    geom_tile(height = 1) +
    scale_fill_manual(values = TUMOR_COLORS, name = "Tumor type") +
    scale_x_continuous(expand = c(0.01, 0)) +
    theme_void(base_size = 10) +
    theme(legend.position = "bottom")

  p <- bars / strip + plot_layout(heights = c(14, 0.7))
  ggsave(out_png, p, width = 14, height = 7, dpi = 150, device = "png")

  base <- sub("\\.png$", "", out_png)
  write.csv(df[, c("x", "USUBJID", "pct", "resp", "resp_adrs")],
            paste0(base, "_order.csv"), row.names = FALSE)
  writeLines(as.character(hlines), paste0(base, "_hlines.txt"))
  invisible(df)
}

# Run as a script: draw the clean figure.
if (sys.nframe() == 0) {
  DATA_DIR   <- "./Data/V1"
  OUTPUT_DIR <- "./Outputs/waterfall_v2"
  dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)
  adsl <- read.csv(file.path(DATA_DIR, "ADSL.csv"), stringsAsFactors = FALSE)
  adtr <- read.csv(file.path(DATA_DIR, "ADTR.csv"), stringsAsFactors = FALSE)
  adrs <- read.csv(file.path(DATA_DIR, "ADRS.csv"), stringsAsFactors = FALSE)
  render_waterfall_v2(prepare_waterfall_v2(adsl, adtr, adrs),
                      file.path(OUTPUT_DIR, "WF_v2_clean.png"))
  message("Saved clean variant-2 waterfall to ", OUTPUT_DIR)
}
