# =============================================================================
# ONCVIZ-001 -- CONSORT diagram, variant 2 (benchmark second implementation)
#
# How this differs from consort_diagram_plot.R (variant 1):
#   * landscape LEFT-TO-RIGHT layout (variant 1 is a portrait top-to-bottom
#     flow): Enrollment -> Allocation -> Discontinuation -> Analysis, one row
#     per arm;
#   * screening and randomization numbers come from ADRAND, arm numbers from
#     ADSL; all box text is written to <png minus .png>_boxes.csv so checkers
#     can read what the figure says without reading pixels;
#   * the two code_patch errors (CN02, CN04) are switches (excl_reasons,
#     disc_detail) instead of edited copies of the script.
#
# DRAFT NOTE: drafted with LLM assistance; an independent author must rewrite
# or review it if the paper claims hand-written benchmark scripts.
# Usage (repository root): Rscript <this file>   (draws the clean figure)
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

cn_reasons <- function(df, col) {
  if (nrow(df) == 0) return("")
  t <- df |> count(.data[[col]], sort = TRUE)
  paste(sprintf("- %s: %d", t[[1]], t$n), collapse = "\n")
}

render_consort_v2 <- function(adrand, adsl, out_png, excl_reasons = TRUE,
                              disc_detail = TRUE) {
  n_screen <- nrow(adrand)
  excl     <- adrand[adrand$SFFL == "Y", ]
  n_rand   <- sum(adrand$RANDFL == "Y")
  arm <- function(a) adsl[adsl$ARM == a, ]
  recv <- function(d) sum(!is.na(d$TRTSDT) & d$TRTSDT != "")

  box <- function(id, cx, cy, w, h, label)
    data.frame(id = id, xmin = cx - w/2, xmax = cx + w/2, ymin = cy - h/2,
               ymax = cy + h/2, label = label, stringsAsFactors = FALSE)

  boxes <- list(
    box("assessed", 10, 55, 17, 12, sprintf("Assessed for\neligibility\n(n = %d)", n_screen)),
    box("excluded", 10, 20, 17, 26,
        if (excl_reasons) sprintf("Excluded\n(n = %d)\n%s", nrow(excl), cn_reasons(excl, "SFREASN"))
        else sprintf("Excluded\n(n = %d)", nrow(excl))),
    box("randomized", 31, 55, 15, 12, sprintf("Randomized\n(n = %d)", n_rand)))
  for (a in c("TREATMENT", "CONTROL")) {
    d <- arm(a); y <- if (a == "TREATMENT") 80 else 30
    disc <- d[d$EOSSTT == "DISCONTINUED", ]
    lab <- if (a == "TREATMENT") "Treatment" else "Control"
    boxes[[length(boxes) + 1]] <- box(paste0("alloc_", a), 52, y, 17, 16,
      sprintf("Allocated to\n%s (n = %d)\nReceived (n = %d)", lab, nrow(d), recv(d)))
    boxes[[length(boxes) + 1]] <- box(paste0("disc_", a), 74, y, 19, 26,
      if (disc_detail) sprintf("Discontinued\n(n = %d)\n%s", nrow(disc), cn_reasons(disc, "DCSREAS"))
      else "Discontinued\nintervention")
    boxes[[length(boxes) + 1]] <- box(paste0("analysed_", a), 92, y, 15, 16,
      sprintf("Analysed\nITT (n = %d)\nSafety (n = %d)", sum(d$ITTFL == "Y"), sum(d$SAFFL == "Y")))
  }
  bx <- bind_rows(boxes)
  g <- function(id) bx[bx$id == id, ]
  seg <- function(a, b, xa, ya, xb, yb) data.frame(x = xa, y = ya, xend = xb, yend = yb)
  arrows_df <- bind_rows(
    seg(0, 0, g("assessed")$xmax, 55, g("randomized")$xmin, 55),
    seg(0, 0, 10, g("assessed")$ymin, 10, g("excluded")$ymax),
    seg(0, 0, g("randomized")$xmax, 55, g("alloc_TREATMENT")$xmin, 80),
    seg(0, 0, g("randomized")$xmax, 55, g("alloc_CONTROL")$xmin, 30),
    seg(0, 0, g("alloc_TREATMENT")$xmax, 80, g("disc_TREATMENT")$xmin, 80),
    seg(0, 0, g("alloc_CONTROL")$xmax, 30, g("disc_CONTROL")$xmin, 30),
    seg(0, 0, g("disc_TREATMENT")$xmax, 80, g("analysed_TREATMENT")$xmin, 80),
    seg(0, 0, g("disc_CONTROL")$xmax, 30, g("analysed_CONTROL")$xmin, 30))
  heads <- data.frame(x = c(10, 31, 52, 74, 92), y = 99,
                      label = c("Enrollment", "Randomization", "Allocation",
                                "Discontinuation", "Analysis"))

  p <- ggplot() +
    geom_rect(data = bx, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
              fill = "#F7FAFC", colour = "#2D3748", linewidth = 0.5) +
    geom_text(data = bx, aes(x = (xmin + xmax)/2, y = (ymin + ymax)/2, label = label),
              size = 2.9, lineheight = 0.95) +
    geom_segment(data = arrows_df, aes(x = x, y = y, xend = xend, yend = yend),
                 arrow = arrow(length = unit(0.1, "inches"), type = "closed"), linewidth = 0.5) +
    geom_text(data = heads, aes(x = x, y = y, label = label), fontface = "bold", size = 3.6) +
    coord_cartesian(xlim = c(0, 100), ylim = c(0, 102)) +
    theme_void() +
    labs(title = "Participant flow (CONSORT 2010 style), left to right")
  ggsave(out_png, p, width = 16, height = 8, dpi = 110, device = "png")
  base <- sub("\\.png$", "", out_png)
  write.csv(bx[, c("id", "label")], paste0(base, "_boxes.csv"), row.names = FALSE)
  invisible(bx)
}

if (sys.nframe() == 0) {
  out <- "./Outputs/consort_v2"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
  render_consort_v2(read.csv("./Data/V1/ADRAND.csv", stringsAsFactors = FALSE),
                    read.csv("./Data/V1/ADSL.csv", stringsAsFactors = FALSE),
                    file.path(out, "CN_v2_clean.png"))
}
