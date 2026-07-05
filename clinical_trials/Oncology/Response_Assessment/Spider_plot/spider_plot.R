library(dplyr)
library(ggplot2)
library(patchwork)

DATA_DIR   <- "./Data/V1"
OUTPUT_DIR <- "./Outputs"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

adsl <- read.csv(file.path(DATA_DIR, "ADSL.csv"), stringsAsFactors = FALSE)
adtr <- read.csv(file.path(DATA_DIR, "ADTR.csv"), stringsAsFactors = FALSE)
adrs <- read.csv(file.path(DATA_DIR, "ADRS.csv"), stringsAsFactors = FALSE)

DAY2MO      <- 30.4375
PR_TH       <- -30
PD_TH       <- 20
CLIP_TOP    <- 100
TUMOR_ORDER <- c("NSCLC","CRC","HCC","PDAC","BRCA")
TUMOR_COLORS <- c(NSCLC="#2166AC", CRC="#1A9641", HCC="#E6007E",
                   PDAC="#D95F02", BRCA="#762A83")
rank_map <- c(CR=1, PR=2, SD=3, PD=4, NE=5)

best_resp_from_adrs <- function(uid) {
  rs <- adrs[adrs$USUBJID == uid & adrs$PARAMCD == "OVRLRESP" &
               adrs$AVALC %in% names(rank_map), ]
  if (nrow(rs) == 0) return("NE")
  rs$AVALC[which.min(rank_map[rs$AVALC])]
}

trt_ids <- adsl$USUBJID[adsl$ARM == "TREATMENT"]
adtr_trt <- adtr[adtr$USUBJID %in% trt_ids, ]
adtr_trt$ADTN <- suppressWarnings(as.numeric(adtr_trt$ADTN))
adtr_trt$PCHG <- suppressWarnings(as.numeric(adtr_trt$PCHG))
adtr_trt$months <- adtr_trt$ADTN / DAY2MO

post    <- adtr_trt[adtr_trt$AVISITN > 0, ]
bl      <- adtr_trt[adtr_trt$AVISITN == 0, ]
multi   <- names(which(table(post$USUBJID) >= 2))
valid   <- intersect(
  bl$USUBJID[!is.na(bl$AVAL) & bl$AVAL > 0],
  multi
)

post_valid <- post[post$USUBJID %in% valid, c("USUBJID","months","PCHG")]
post_valid$PCHG <- pmin(post_valid$PCHG, CLIP_TOP, na.rm=TRUE)
bl_rows <- data.frame(USUBJID=valid, months=0.0, PCHG=0.0,
                      stringsAsFactors=FALSE)
spider_base <- bind_rows(bl_rows, post_valid)

meta <- adsl[adsl$ARM == "TREATMENT", c("USUBJID","TUMORTYPE","TMBHIGH")]
meta$BESTRESP <- sapply(meta$USUBJID, best_resp_from_adrs)

spider <- spider_base |>
  left_join(meta, by="USUBJID") |>
  arrange(USUBJID, months)

xmax_data <- max(spider$months, na.rm=TRUE)
XMAX      <- ceiling(xmax_data / 3) * 3 + 3

make_spider_panel <- function(sub_df, panel_title, color_col,
                               color_map, show_y=FALSE) {
  n     <- n_distinct(sub_df$USUBJID)
  n_cr  <- sum(meta$BESTRESP[meta$USUBJID %in% unique(sub_df$USUBJID)] == "CR")
  n_pr  <- sum(meta$BESTRESP[meta$USUBJID %in% unique(sub_df$USUBJID)] == "PR")
  orr   <- round((n_cr + n_pr) / n * 100)
  stats <- sprintf("N=%d   CR:%d  PR:%d  ORR=%d%%", n, n_cr, n_pr, orr)

  sub_df[[color_col]] <- as.character(sub_df[[color_col]])

  p <- ggplot(sub_df, aes(x=months, y=PCHG,
                           group=USUBJID,
                           color=.data[[color_col]])) +
    geom_hline(yintercept=0,     linewidth=1.1, color="#222222") +
    geom_hline(yintercept=PR_TH, linewidth=0.65, color="#555555",
               linetype="dashed") +
    geom_hline(yintercept=PD_TH, linewidth=0.65, color="#555555",
               linetype="dashed") +
    annotate("text", x=XMAX+0.4, y=PD_TH+2.5,
             label="+20% PD", hjust=0, vjust=0,
             size=2.6, fontface="italic", color="#555555") +
    annotate("text", x=XMAX+0.4, y=PR_TH-2.5,
             label="-30% PR", hjust=0, vjust=1,
             size=2.6, fontface="italic", color="#555555") +
    geom_line(linewidth=0.85, alpha=0.72,
              lineend="butt", linejoin="miter") +
    geom_point(data=sub_df |>
                 group_by(USUBJID) |>
                 slice_max(months, n=1, with_ties=FALSE) |>
                 ungroup(),
               aes(x=months, y=PCHG),
               size=1.6, alpha=0.88) +
    scale_color_manual(values=color_map) +
    scale_x_continuous(
      breaks = seq(0, XMAX, by=6),
      labels = as.character(seq(0, XMAX, by=6)),
      limits = c(-0.6, XMAX+3.5),
      expand = c(0,0)
    ) +
    scale_y_continuous(
      breaks = seq(-100, 120, by=20),
      labels = function(x) paste0(x, "%"),
      limits = c(-105, 130)
    ) +
    labs(
      title    = panel_title,
      subtitle = stats,
      x        = "Time from Treatment Start (months)",
      y        = if (show_y) "Change from Baseline in SLD (%)" else NULL
    ) +
    theme_classic(base_size=9) +
    theme(
      plot.title       = element_text(face="bold", size=11,
                                      hjust=0.5, color="black",
                                      margin=margin(b=2)),
      plot.subtitle    = element_text(size=6.8, hjust=0.5,
                                      color="#555555", margin=margin(b=4)),
      axis.title.x     = element_text(face="bold", size=8.5),
      axis.title.y     = if (show_y) element_text(face="bold", size=8.5)
                         else element_blank(),
      axis.text        = element_text(size=7.5),
      axis.text.y      = if (show_y) element_text(size=7.5)
                         else element_blank(),
      axis.ticks.y     = if (show_y) element_line() else element_blank(),
      axis.line        = element_line(color="#bbbbbb", linewidth=0.5),
      panel.background = element_rect(fill="#fafafa", color=NA),
      panel.grid.major.y = element_line(color="#eeeeee", linewidth=0.4),
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      legend.position    = "none",
      plot.margin        = margin(4, 4, 4, 4),
      plot.clip          = "off"
    )

  p
}

panels_A <- lapply(seq_along(TUMOR_ORDER), function(i) {
  tumor <- TUMOR_ORDER[i]
  sub   <- spider[spider$TUMORTYPE == tumor, ]
  title_col <- names(TUMOR_COLORS)[names(TUMOR_COLORS) == tumor]
  make_spider_panel(
    sub_df      = sub,
    panel_title = tumor,
    color_col   = "BESTRESP",
    color_map   = c(CR="#1A6B3C", PR="#2B83BA", SD="#F0A500",
                    PD="#C0392B", NE="#999999"),
    show_y      = (i == 1)
  ) + theme(plot.title = element_text(color=TUMOR_COLORS[tumor]))
})

shared_x_label <- ggplot() +
  annotate("text", x=0.5, y=0.5,
           label="Time from Treatment Start (months)",
           hjust=0.5, vjust=0.5, size=3.5, fontface="bold") +
  theme_void() + theme(plot.margin=margin(0,0,2,0))

legend_A <- ggplot(
  data.frame(
    resp   = c("CR","PR","SD","PD"),
    label  = c("Complete Response (CR)","Partial Response (PR)",
               "Stable Disease (SD)",  "Progressive Disease (PD)"),
    color  = c("#1A6B3C","#2B83BA","#F0A500","#C0392B")
  )) +
  geom_point(aes(x=resp, y=1, color=resp), shape=15, size=5) +
  scale_color_manual(
    values = c(CR="#1A6B3C", PR="#2B83BA", SD="#F0A500", PD="#C0392B"),
    labels = c(CR="Complete Response (CR)", PR="Partial Response (PR)",
               SD="Stable Disease (SD)",    PD="Progressive Disease (PD)"),
    name   = "Best Overall Response"
  ) +
  guides(color=guide_legend(nrow=1, override.aes=list(size=4))) +
  theme_void() +
  theme(
    legend.position   = "bottom",
    legend.text       = element_text(size=8),
    legend.title      = element_text(size=8, face="bold"),
    legend.key.size   = unit(0.4,"cm"),
    legend.margin     = margin(t=4)
  )

for (i in seq_along(panels_A))
  panels_A[[i]] <- panels_A[[i]] + labs(x=NULL)

fig_A <- wrap_plots(panels_A, nrow=1) +
  plot_annotation(
    title = "Spider Plot \u2014 Percent Change from Baseline in SLD by Tumor Type  \u00b7  ONCVIZ-001  \u00b7  Treatment Arm",
    caption = "RECIST 1.1  \u00b7  Lines clipped at +100%  \u00b7  Data cutoff: 05 Mar 2026",
    theme = theme(
      plot.title   = element_text(face="bold", size=11, hjust=0.5,
                                  margin=margin(b=4)),
      plot.caption = element_text(size=8, hjust=0.5, color="#555555",
                                  face="italic", margin=margin(t=4))
    )
  )

get_legend <- function(p) {
  g <- ggplotGrob(p)
  g$grobs[[which(sapply(g$grobs, function(x) x$name) == "guide-box")]]
}

png(file.path(OUTPUT_DIR, "spider_plot_A_tumor_type.png"),
    width=22, height=8.2, units="in", res=300)
gridExtra::grid.arrange(
  ggplotGrob(fig_A),
  grid::textGrob("Time from Treatment Start (months)",
                 gp=grid::gpar(fontsize=9, fontface="bold")),
  get_legend(legend_A),
  nrow=3, heights=c(7.2, 0.3, 0.7)
)
dev.off()
message("  Saved: spider_plot_A_tumor_type.png")

panels_B <- lapply(seq_along(c("Y","N")), function(i) {
  vals <- c("Y","N")
  lbls <- c("TMB-High","TMB-Low")
  sub  <- spider[spider$TMBHIGH == vals[i], ]
  make_spider_panel(
    sub_df      = sub,
    panel_title = lbls[i],
    color_col   = "TUMORTYPE",
    color_map   = TUMOR_COLORS,
    show_y      = (i == 1)
  ) + theme(
    plot.title = element_text(
      color = c(Y="#C0392B", N="#2166AC")[vals[i]])
  )
})

legend_B <- ggplot(
  data.frame(tumor=TUMOR_ORDER,
             color=TUMOR_COLORS[TUMOR_ORDER])) +
  geom_point(aes(x=tumor, y=1, color=tumor), shape=15, size=5) +
  scale_color_manual(values=TUMOR_COLORS, name="Tumor Type") +
  guides(color=guide_legend(nrow=1, override.aes=list(size=4))) +
  theme_void() +
  theme(
    legend.position = "bottom",
    legend.text     = element_text(size=8),
    legend.title    = element_text(size=8, face="bold"),
    legend.key.size = unit(0.4,"cm"),
    legend.margin   = margin(t=4)
  )

for (i in seq_along(panels_B))
  panels_B[[i]] <- panels_B[[i]] + labs(x=NULL)

fig_B <- wrap_plots(panels_B, nrow=1) +
  plot_annotation(
    title = "Spider Plot \u2014 Percent Change from Baseline in SLD by TMB Status  \u00b7  ONCVIZ-001  \u00b7  Treatment Arm",
    caption = "RECIST 1.1  \u00b7  Lines clipped at +100%  \u00b7  Data cutoff: 05 Mar 2026",
    theme = theme(
      plot.title   = element_text(face="bold", size=11, hjust=0.5,
                                  margin=margin(b=4)),
      plot.caption = element_text(size=8, hjust=0.5, color="#555555",
                                  face="italic", margin=margin(t=4))
    )
  )

png(file.path(OUTPUT_DIR, "spider_plot_B_tmb.png"),
    width=15, height=8.2, units="in", res=300)
gridExtra::grid.arrange(
  ggplotGrob(fig_B),
  grid::textGrob("Time from Treatment Start (months)",
                 gp=grid::gpar(fontsize=9, fontface="bold")),
  get_legend(legend_B),
  nrow=3, heights=c(7.2, 0.3, 0.7)
)
dev.off()
message("  Saved: spider_plot_B_tmb.png")

message("All spider plots saved to: ", OUTPUT_DIR)
