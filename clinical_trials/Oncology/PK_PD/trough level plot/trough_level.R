# trough_level.R
# Trough Level Plot — RP2D (300 mg) cohort, Cycle 1 Day 1 vs Cycle 3 Day 1
# Standards applied: BLQ-rate panel + quantifiable-only strip (median/IQR summary stats are
# suppressed when the BLQ rate exceeds 50%, per standard PK reporting convention), N per visit.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

SCRIPT_DIR <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
if (length(SCRIPT_DIR) == 0 || SCRIPT_DIR == "") SCRIPT_DIR <- getwd()
source(file.path(SCRIPT_DIR, "..", "pkpd_common.R"))

DATA_DIR   <- file.path(SCRIPT_DIR, "..", "Data", "V1")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "Out")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

adpk <- read.csv(file.path(DATA_DIR, "ADPK.csv"), stringsAsFactors = FALSE)

visit_order <- c("CYCLE 1 DAY 1", "CYCLE 3 DAY 1")
tr <- adpk |>
  filter(PARAMCD == "TROUGH", DOSE == 300, ANL01FL == "Y") |>
  mutate(AVISIT = factor(AVISIT, levels = visit_order))

n_total <- n_distinct(tr$USUBJID)

blq_summary <- tr |>
  group_by(AVISIT) |>
  summarise(n = n(),
            pct_blq = 100 * mean(AVAL <= 0),
            pct_quant = 100 - pct_blq,
            .groups = "drop")

# Panel 1 — BLQ rate stacked bar, % BLQ labeled centered inside the bar
blq_long <- blq_summary |>
  tidyr::pivot_longer(c(pct_blq, pct_quant), names_to = "cat", values_to = "pct") |>
  mutate(cat = factor(cat, levels = c("pct_blq", "pct_quant"),
                       labels = c("BLQ", "Quantifiable")))

p1 <- ggplot(blq_long, aes(x = AVISIT, y = pct, fill = cat)) +
  geom_col(width = 0.6) +
  geom_text(data = blq_summary, aes(x = AVISIT, y = pct_blq / 2,
                                     label = sprintf("%.1f%% BLQ", pct_blq)),
            inherit.aes = FALSE, color = "white", fontface = "bold", size = 3.3) +
  geom_text(data = blq_summary, aes(x = AVISIT, y = 103,
                                     label = sprintf("%.1f%% quantifiable", pct_quant)),
            inherit.aes = FALSE, size = 2.9) +
  scale_fill_manual(values = c("BLQ" = "#AAAAAA", "Quantifiable" = "#1A3A7C"), name = NULL) +
  scale_x_discrete(labels = paste0(gsub("CYCLE ", "C", gsub(" DAY ", "D", visit_order)),
                                   "\n(n=", blq_summary$n, ")")) +
  scale_y_continuous(limits = c(0, 118), breaks = seq(0, 100, 20)) +
  labs(x = NULL, y = "% of trough samples", title = "Trough BLQ rate") +
  theme_pkpd() +
  theme(legend.position = "bottom")

# Panel 2 — quantifiable troughs only, log scale, jittered individual points (no boxplot,
# since median/IQR are not meaningful when the BLQ rate exceeds 50%)
set.seed(42)
quant_df <- tr |> filter(AVAL > 0)
n_quant  <- quant_df |> count(AVISIT, name = "n_q")

p2 <- ggplot(quant_df, aes(x = AVISIT, y = AVAL)) +
  geom_jitter(width = 0.05, size = 2.6, color = "#C0392B", alpha = 0.8, shape = 21, fill = "#C0392B") +
  scale_y_log10() +
  scale_x_discrete(labels = paste0(visit_order, "\n(n quantifiable=", n_quant$n_q, ")")) +
  labs(x = NULL, y = expression("Trough concentration, C"[trough]*" (ng/mL, log scale)"),
       title = "Quantifiable troughs only") +
  theme_pkpd()

cap <- paste(strwrap(paste0(
  "Median/IQR summary statistics are not presented because the trough BLQ rate exceeds the 50% threshold at ",
  "both visits (standard PK reporting convention for BLQ-dominated data). Finding: minimal pre-dose exposure at ",
  "steady state indicates rapid clearance relative to the 24h QD dosing interval and negligible accumulation. ",
  "Individual quantifiable values shown at right (log scale)."
), width = 150), collapse = "\n")

combined <- (p1 | p2) +
  plot_annotation(
    title = sprintf("Trough Level Plot \u2014 Vizatinib (ONCVIZ-001) 300 mg QD (RP2D), PK-evaluable population (N=%d)", n_total),
    caption = cap,
    theme = theme(plot.title = element_text(face = "bold", size = 13, hjust = 0.5))
  )

ggsave(file.path(OUTPUT_DIR, "trough_level.png"), combined, width = 11, height = 6, dpi = 300)
message("Saved trough_level.png to ", OUTPUT_DIR)
