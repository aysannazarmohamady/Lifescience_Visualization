library(dplyr)
library(ggplot2)
library(patchwork)
library(grid)
library(gridExtra)

DATA_DIR   <- "./Data/V1"
OUTPUT_DIR <- "./Outputs"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

adsl <- read.csv(file.path(DATA_DIR, "ADSL.csv"), stringsAsFactors = FALSE)
adtr <- read.csv(file.path(DATA_DIR, "ADTR.csv"), stringsAsFactors = FALSE)
adrs <- read.csv(file.path(DATA_DIR, "ADRS.csv"), stringsAsFactors = FALSE)

RESP_COLORS <- c(CR="#1A3A7C", PR="#4A90C4", SD="#E8A020", PD="#C0392B", NE="#AAAAAA")
LM_COLORS   <- c(Y="#C0392B", N="#2B6CB0")
TUMOR_ORDER <- c("NSCLC","BRCA","HCC","CRC","PDAC")
PR_TH <- -30; PD_TH <- 20

rank_map <- c(CR=1, PR=2, SD=3, PD=4, NE=5)

best_resp <- function(uid) {
  rs <- adrs[adrs$USUBJID == uid & adrs$PARAMCD == "OVRLRESP" &
               adrs$AVALC %in% names(rank_map), ]
  if (nrow(rs) == 0) return("NE")
  rs$AVALC[which.min(rank_map[rs$AVALC])]
}

best_pchg <- function(uid) {
  post <- adtr[adtr$USUBJID == uid & adtr$AVISITN > 0, ]
  vals <- suppressWarnings(as.numeric(post$PCHG))
  vals <- vals[!is.na(vals)]
  if (length(vals) == 0) return(0)
  min(vals)
}

wf <- adsl[adsl$ARM == "TREATMENT", ] |>
  rowwise() |>
  mutate(
    short = sub("^[^-]+-[^-]+-", "", USUBJID),
    resp  = best_resp(USUBJID),
    pct   = best_pchg(USUBJID),
    lm    = trimws(LIVERMETS)
  ) |>
  ungroup()

make_panel <- function(df_tumor, tumor_name, show_y = FALSE) {
  df <- df_tumor |>
    arrange(desc(pct)) |>
    mutate(x = row_number())

  n    <- nrow(df)
  n_cr <- sum(df$resp == "CR")
  n_pr <- sum(df$resp == "PR")
  orr  <- round((n_cr + n_pr) / n * 100)

  ttl <- sprintf("%s  (N=%d)\nCR:%d  PR:%d  ORR:%d%%", tumor_name, n, n_cr, n_pr, orr)

  p <- ggplot(df, aes(x = x, y = pct, fill = resp)) +
    geom_col(width = 0.78, color = NA) +
    geom_hline(yintercept = 0,     linewidth = 1.0, color = "black") +
    geom_hline(yintercept = PD_TH, linewidth = 0.7, color = "#555555",
               linetype = "dashed") +
    geom_hline(yintercept = PR_TH, linewidth = 0.7, color = "#555555",
               linetype = "dashed") +
    scale_fill_manual(values = RESP_COLORS, name = "Best Response",
                      labels = c(CR="Complete Response (CR)",
                                 PR="Partial Response (PR)",
                                 SD="Stable Disease (SD)",
                                 PD="Progressive Disease (PD)")) +
    scale_x_continuous(breaks = seq_len(n), labels = df$short,
                       expand = c(0.02, 0)) +
    scale_y_continuous(breaks = seq(-100, 140, 20),
                       labels = function(x) paste0(x, "%"),
                       limits = c(-118, 150)) +
    labs(title = ttl, x = NULL,
         y = if (show_y) "Best % Change\nfrom Baseline in SLD" else NULL) +
    theme_classic(base_size = 9) +
    theme(
      plot.title        = element_text(face = "bold", size = 10.5,
                                       hjust = 0.5, lineheight = 1.3,
                                       margin = margin(b = 6)),
      axis.text.x       = element_text(angle = 90, vjust = 0.5, hjust = 1,
                                       size = 7, family = "mono"),
      axis.text.y       = if (show_y) element_text(size = 9) else element_blank(),
      axis.ticks.y      = if (show_y) element_line() else element_blank(),
      axis.line.x       = element_blank(),
      axis.line.y       = element_line(color = "black"),
      panel.grid.major.y = element_line(color = "#e8e8e8", linewidth = 0.45),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position   = "none",
      plot.margin       = margin(4, 4, 2, 4)
    )

  if (!show_y) p <- p + theme(axis.title.y = element_blank())

  list(bar = p, df = df, n = n)
}

make_strip <- function(df_sorted, show_label = FALSE) {
  df <- df_sorted |> mutate(x = row_number())
  n  <- nrow(df)

  p <- ggplot(df, aes(x = x, y = 1, fill = lm)) +
    geom_col(width = 0.78, color = NA) +
    scale_fill_manual(values = LM_COLORS,
                      labels = c(Y = "Yes", N = "No"),
                      name = "Liver mets") +
    scale_x_continuous(breaks = seq_len(n), expand = c(0.02, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin     = margin(0, 4, 0, 4)
    )

  if (show_label) {
    p <- p + annotate("text", x = 0, y = 0.5,
                      label = "Liver\nmets", hjust = 1.1,
                      size = 3.0, fontface = "bold", color = "#333333")
  }
  p
}

panels <- lapply(seq_along(TUMOR_ORDER), function(i) {
  tumor <- TUMOR_ORDER[i]
  sub   <- wf[wf$tumor == tumor, ]
  make_panel(sub, tumor, show_y = (i == 1))
})

bar_plots   <- lapply(panels, `[[`, "bar")
strip_plots <- lapply(seq_along(panels), function(i) {
  make_strip(panels[[i]]$df, show_label = (i == 1))
})

combined_cols <- lapply(seq_along(TUMOR_ORDER), function(i) {
  bar_plots[[i]] / strip_plots[[i]] +
    plot_layout(heights = c(10, 0.8))
})

full_plot <- wrap_plots(combined_cols, nrow = 1) +
  plot_annotation(
    title = paste0(
      "Waterfall Plot \u2014 Best % Change from Baseline in SLD",
      " by Tumor Type  \u00b7  ONCVIZ-001  \u00b7  Treatment Arm"
    ),
    theme = theme(
      plot.title = element_text(face = "bold", size = 13,
                                hjust = 0.5,
                                margin = margin(b = 18))
    )
  )

legend_resp <- ggplot(data.frame(resp = names(RESP_COLORS)[1:4]),
                      aes(x = resp, fill = resp)) +
  geom_col(aes(y = 1)) +
  scale_fill_manual(
    values = RESP_COLORS[1:4],
    labels = c(CR="Complete Response (CR)", PR="Partial Response (PR)",
               SD="Stable Disease (SD)",   PD="Progressive Disease (PD)"),
    name = NULL
  ) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 4))) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text     = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"))

legend_lm <- ggplot(data.frame(lm = c("Y","N")), aes(x = lm, fill = lm)) +
  geom_col(aes(y = 1)) +
  scale_fill_manual(
    values = LM_COLORS,
    labels = c(Y="Liver mets: Yes", N="Liver mets: No"),
    name = NULL
  ) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 4))) +
  theme_void() +
  theme(legend.position = "bottom",
        legend.text     = element_text(size = 9),
        legend.key.size = unit(0.45, "cm"))

leg_grob <- gridExtra::arrangeGrob(
  get_legend <- function(p) {
    tmp <- ggplotGrob(p)
    leg <- tmp$grobs[[which(sapply(tmp$grobs, function(x) x$name) == "guide-box")]]
    leg
  },
  get_legend(legend_resp),
  get_legend(legend_lm),
  nrow = 1
)

get_legend <- function(p) {
  g <- ggplotGrob(p)
  g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
}
leg_combined <- gridExtra::arrangeGrob(
  get_legend(legend_resp), get_legend(legend_lm),
  nrow = 1, widths = c(4, 1.5)
)

png(file.path(OUTPUT_DIR, "waterfall_5panel.png"),
    width = 28, height = 11, units = "in", res = 150)
gridExtra::grid.arrange(
  ggplotGrob(full_plot),
  leg_combined,
  nrow   = 2,
  heights = c(10.3, 0.7)
)
dev.off()

make_waterfall_single <- function(df_in, title_str, fname) {
  df <- df_in |>
    arrange(desc(pct)) |>
    mutate(x = row_number())

  n    <- nrow(df)
  n_cr <- sum(df$resp == "CR"); n_pr <- sum(df$resp == "PR")
  n_sd <- sum(df$resp == "SD"); n_pd <- sum(df$resp == "PD")
  orr  <- round((n_cr + n_pr) / n * 100)
  stats_lbl <- sprintf("N=%d   CR:%d  PR:%d  SD:%d  PD:%d   ORR=%d%%",
                       n, n_cr, n_pr, n_sd, n_pd, orr)

  bar_p <- ggplot(df, aes(x = x, y = pct, fill = resp)) +
    geom_col(width = 0.78, color = NA) +
    geom_hline(yintercept = 0,     linewidth = 1.2, color = "black") +
    geom_hline(yintercept = PD_TH, linewidth = 0.8, color = "#444444",
               linetype = "dashed") +
    geom_hline(yintercept = PR_TH, linewidth = 0.8, color = "#444444",
               linetype = "dashed") +
    annotate("text", x = n + 0.7, y = PD_TH + 2,
             label = "+20%: PD threshold", hjust = 0, vjust = 0,
             size = 3.2, fontface = "italic", color = "#444444") +
    annotate("text", x = n + 0.7, y = PR_TH - 2,
             label = "-30%: PR threshold", hjust = 0, vjust = 1,
             size = 3.2, fontface = "italic", color = "#444444") +
    scale_fill_manual(
      values = RESP_COLORS,
      labels = c(CR="Complete Response (CR)", PR="Partial Response (PR)",
                 SD="Stable Disease (SD)",    PD="Progressive Disease (PD)"),
      name = "Best Response"
    ) +
    scale_x_continuous(breaks = seq_len(n), labels = df$short,
                       expand = c(0.01, 0)) +
    scale_y_continuous(breaks = seq(-100, 140, 20),
                       labels = function(x) paste0(x, "%"),
                       limits = c(-118, 155)) +
    labs(title = paste0(title_str, "\n", stats_lbl),
         x = NULL,
         y = "Best % Change from Baseline in SLD") +
    theme_classic(base_size = 10) +
    theme(
      plot.title         = element_text(face="bold", size=12,
                                        hjust=0.5, lineheight=1.3,
                                        margin=margin(b=8)),
      axis.text.x        = element_text(angle=90, vjust=0.5, hjust=1,
                                        size=8, family="mono"),
      axis.text.y        = element_text(size=10),
      panel.grid.major.y = element_line(color="#e8e8e8", linewidth=0.5),
      panel.grid.minor   = element_blank(),
      panel.grid.major.x = element_blank(),
      legend.position    = "inside",
      legend.position.inside = c(0.98, 0.98),
      legend.justification   = c(1, 1),
      legend.background  = element_rect(fill="white", color="#aaaaaa",
                                        linewidth=0.5),
      legend.text        = element_text(size=9),
      legend.title       = element_text(size=9, face="bold"),
      legend.key.size    = unit(0.45, "cm"),
      plot.margin        = margin(6, 6, 2, 6)
    )

  strip_p <- ggplot(df |> mutate(x = row_number()),
                    aes(x = x, y = 1, fill = lm)) +
    geom_col(width = 0.78, color = NA) +
    scale_fill_manual(values  = LM_COLORS,
                      labels  = c(Y="Yes", N="No"),
                      name    = "Liver mets") +
    scale_x_continuous(expand = c(0.01, 0)) +
    scale_y_continuous(expand = c(0, 0)) +
    annotate("text", x = 0, y = 0.5,
             label = "Liver mets", hjust = 1.05, size = 3.2,
             fontface = "bold", color = "#333333") +
    theme_void() +
    theme(
      legend.position = "none",
      plot.margin     = margin(0, 6, 0, 6)
    )

  lm_leg_p <- ggplot(data.frame(lm = c("Y","N")), aes(x=lm, fill=lm)) +
    geom_col(aes(y=1)) +
    scale_fill_manual(values=LM_COLORS,
                      labels=c(Y="Yes", N="No"), name="Liver mets") +
    guides(fill=guide_legend(nrow=1, override.aes=list(size=3.5))) +
    theme_void() +
    theme(legend.position="bottom",
          legend.text=element_text(size=9),
          legend.key.size=unit(0.42,"cm"))

  lm_leg_grob <- get_legend(lm_leg_p)

  full <- bar_p / strip_p + plot_layout(heights = c(10, 0.75))

  png(fname, width=16, height=9, units="in", res=150)
  gridExtra::grid.arrange(
    ggplotGrob(full),
    lm_leg_grob,
    nrow=2, heights=c(9.3, 0.7)
  )
  dev.off()
  message("  Saved: ", basename(fname))
}

make_waterfall_single(wf,
  "Waterfall Plot \u2014 Best % Change from Baseline  \u00b7  ONCVIZ-001  \u00b7  Treatment Arm",
  file.path(OUTPUT_DIR, "waterfall_all_treatment.png"))

for (tumor in TUMOR_ORDER) {
  sub <- wf[wf$tumor == tumor, ]
  make_waterfall_single(sub,
    sprintf("Waterfall Plot \u2014 Best %% Change from Baseline  \u00b7  ONCVIZ-001  \u00b7  %s  \u00b7  Treatment Arm", tumor),
    file.path(OUTPUT_DIR, sprintf("waterfall_%s.png", tolower(tumor))))
}

message("All waterfall plots saved to: ", OUTPUT_DIR)
