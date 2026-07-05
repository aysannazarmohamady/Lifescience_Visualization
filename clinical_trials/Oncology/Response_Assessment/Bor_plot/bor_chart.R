library(dplyr)
library(ggplot2)
library(patchwork)
library(grid)
library(gridExtra)

DATA_DIR   <- "./Data/V1"
OUTPUT_DIR <- "./Outputs"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

adsl <- read.csv(file.path(DATA_DIR, "ADSL.csv"), stringsAsFactors = FALSE)
CUTOFF <- "05 Mar 2026"

trt <- adsl |>
  filter(ARM == "TREATMENT") |>
  mutate(COHORT = case_when(
    DOSELEVEL == 100 ~ "Cohort 1\n(100 mg)",
    DOSELEVEL == 200 ~ "Cohort 2\n(200 mg)",
    DOSELEVEL == 300 ~ "RP2D\n(300 mg)",
    DOSELEVEL == 400 ~ "Cohort 4\n(400 mg)"
  ))

COHORT_ORDER  <- c("Cohort 1\n(100 mg)","Cohort 2\n(200 mg)",
                   "RP2D\n(300 mg)","Cohort 4\n(400 mg)")
RESP_ORDER    <- c("CR","PR","SD","PD")
RESP_COLORS   <- c(CR="#1A9641", PR="#4A90C4", SD="#E8A020", PD="#C0392B")
COHORT_COLORS <- c("Cohort 1\n(100 mg)"="#E66101",
                   "Cohort 2\n(200 mg)"="#5AAE61",
                   "RP2D\n(300 mg)"    ="#3288BD",
                   "Cohort 4\n(400 mg)"="#9970AB")

cohort_pcts <- trt |>
  group_by(COHORT) |>
  mutate(n_cohort = n()) |>
  group_by(COHORT, BESTRSPC, n_cohort) |>
  summarise(count = n(), .groups = "drop") |>
  mutate(pct = round(count / n_cohort * 100)) |>
  filter(BESTRSPC %in% RESP_ORDER) |>
  mutate(COHORT    = factor(COHORT,    levels = COHORT_ORDER),
         BESTRSPC  = factor(BESTRSPC,  levels = RESP_ORDER))

n_labels <- trt |>
  group_by(COHORT) |>
  summarise(n = n(), .groups = "drop") |>
  mutate(COHORT = factor(COHORT, levels = COHORT_ORDER),
         label  = paste0("n=", n))

# ── LEFT: stacked 100% ────────────────────────────────────────────────────────
p_stacked <- ggplot(cohort_pcts,
                    aes(x = COHORT, y = pct, fill = BESTRSPC)) +
  geom_col(width = 0.55, position = "stack") +
  geom_text(aes(label = ifelse(pct >= 7, paste0(pct, "%"), "")),
            position = position_stack(vjust = 0.5),
            color = "white", fontface = "bold", size = 3.6) +
  geom_text(data = n_labels,
            aes(x = COHORT, y = 103, label = label),
            inherit.aes = FALSE,
            fontface = "bold", size = 3.6, color = "#333333") +
  scale_fill_manual(values = RESP_COLORS,
                    breaks = RESP_ORDER,
                    name   = "Response") +
  scale_y_continuous(limits = c(0, 114),
                     breaks = seq(0, 100, 20),
                     labels = function(x) paste0(x, "%")) +
  scale_x_discrete(labels = function(x) x) +
  labs(title = "Best Overall Response by Cohort",
       x     = "Dose Cohort",
       y     = "Patients (%)") +
  theme_classic(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12, hjust = 0.5,
                                    margin = margin(b = 8)),
    axis.title.x     = element_text(face = "bold", size = 11,
                                    margin = margin(t = 18)),
    axis.title.y     = element_text(face = "bold", size = 11),
    axis.text.x      = element_text(face = "bold", size = 10,
                                    lineheight = 1.3),
    axis.text.y      = element_text(size = 9),
    axis.ticks       = element_blank(),
    axis.line        = element_line(color = "#bbbbbb"),
    panel.grid.major.y = element_line(color = "#eeeeee", linewidth = 0.55),
    panel.grid.major.x = element_blank(),
    legend.position  = "none",
    plot.margin      = margin(6, 8, 6, 6)
  )

# ── RIGHT: grouped bar ────────────────────────────────────────────────────────
p_grouped <- ggplot(cohort_pcts,
                    aes(x = BESTRSPC, y = pct, fill = COHORT)) +
  geom_col(position = position_dodge(width = 0.72),
           width = 0.68, alpha = 0.92) +
  scale_fill_manual(values = COHORT_COLORS,
                    breaks = COHORT_ORDER,
                    labels = gsub("\n", " ", COHORT_ORDER),
                    name   = NULL) +
  scale_y_continuous(limits = c(0, 85),
                     breaks = seq(0, 80, 10)) +
  labs(title = "Response Distribution Across Cohorts",
       x     = "Best Overall Response",
       y     = "Patients (%)") +
  theme_classic(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 12, hjust = 0.5,
                                    margin = margin(b = 8)),
    axis.title.x     = element_text(face = "bold", size = 11,
                                    margin = margin(t = 18)),
    axis.title.y     = element_text(face = "bold", size = 11),
    axis.text.x      = element_text(face = "bold", size = 11),
    axis.text.y      = element_text(size = 9),
    axis.ticks       = element_blank(),
    axis.line        = element_line(color = "#bbbbbb"),
    panel.grid.major.y = element_line(color = "#eeeeee", linewidth = 0.55),
    panel.grid.major.x = element_blank(),
    legend.position  = "none",
    plot.margin      = margin(6, 8, 6, 6)
  )

# ── Shared legends ────────────────────────────────────────────────────────────
get_legend <- function(p) {
  g   <- ggplotGrob(p)
  idx <- which(sapply(g$grobs, function(x) x$name) == "guide-box")
  g$grobs[[idx]]
}

leg_resp <- get_legend(
  p_stacked + theme(
    legend.position    = "right",
    legend.title       = element_text(face = "bold", size = 9),
    legend.text        = element_text(size = 9),
    legend.key.size    = unit(0.42, "cm"),
    legend.margin      = margin(0, 0, 0, 4),
    legend.box.margin  = margin(0)
  )
)

leg_cohort <- get_legend(
  p_grouped + theme(
    legend.position    = "right",
    legend.text        = element_text(size = 9),
    legend.key.size    = unit(0.42, "cm"),
    legend.margin      = margin(0, 0, 0, 4),
    legend.box.margin  = margin(0)
  )
)

# ── Summary bar ───────────────────────────────────────────────────────────────
n_all   <- nrow(trt)
orr_n   <- sum(trt$BESTRSPC %in% c("CR","PR"))
orr_pct <- round(orr_n / n_all * 100)
resp_summary <- trt |>
  count(BESTRSPC) |>
  filter(BESTRSPC %in% RESP_ORDER) |>
  mutate(pct = round(n / n_all * 100),
         lbl = sprintf("%s  %d (%d%%)", BESTRSPC, n, pct))

summary_txt <- paste(
  sprintf("Total  %d", n_all),
  sprintf("ORR  %d%%  (%d)", orr_pct, orr_n),
  paste(resp_summary$lbl, collapse = "     "),
  sep = "     "
)

sum_grob <- grid::rectGrob(gp = grid::gpar(fill = "#1F4E79", col = NA))
sum_txt  <- grid::textGrob(
  summary_txt,
  gp = grid::gpar(col = "white", fontsize = 10, fontface = "bold",
                  fontfamily = "mono")
)
sum_panel <- gridExtra::arrangeGrob(sum_grob, sum_txt, nrow = 1,
                                     widths = c(1, 0))

# ── Assemble ──────────────────────────────────────────────────────────────────
main_title <- grid::textGrob(
  sprintf("Best Overall Response \u2014 ONCVIZ-001 \u00b7 Vizatinib \u00b7 Treatment Arm  (N=%d)", n_all),
  gp = grid::gpar(fontsize = 14, fontface = "bold")
)

legends_col <- gridExtra::arrangeGrob(
  leg_resp,
  grid::rectGrob(gp = grid::gpar(fill = NA, col = NA)),
  leg_cohort,
  nrow    = 3,
  heights = c(4, 0.5, 4)
)

two_panels <- gridExtra::arrangeGrob(
  ggplotGrob(p_stacked),
  ggplotGrob(p_grouped),
  nrow   = 1,
  widths = c(1, 1)
)

full_layout <- gridExtra::arrangeGrob(
  main_title,
  gridExtra::arrangeGrob(two_panels, legends_col,
                          nrow   = 1,
                          widths = c(10, 1.4)),
  gridExtra::arrangeGrob(
    grid::rectGrob(gp = grid::gpar(fill = "#1F4E79", col = NA)),
    grid::textGrob(summary_txt,
                   gp = grid::gpar(col="white", fontsize=10,
                                   fontface="bold", fontfamily="mono")),
    nrow   = 1,
    widths = c(0.01, 1)
  ),
  nrow    = 3,
  heights = c(0.8, 8.5, 0.8)
)

png(file.path(OUTPUT_DIR, "bor_plot.png"),
    width = 20, height = 9, units = "in", res = 150)
grid::grid.draw(full_layout)
dev.off()

message("Saved: bor_plot.png")
