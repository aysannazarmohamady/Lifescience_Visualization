library(tidyverse)
library(patchwork)

OUTPUT_A   <- "spider_plot_A_tumor_type.png"
OUTPUT_B   <- "spider_plot_B_tmb.png"
DPI        <- 300
DAY2MO     <- 30.4375
PR_TH      <- -30
PD_TH      <- 20
YMIN       <- -105
YMAX       <- 125
CLIP_TOP   <- 120
XMAX       <- 48
CUTOFF_LBL <- "Data cutoff: 05 Mar 2026"

TUMOR_ORDER  <- c("NSCLC","CRC","HCC","PDAC","BRCA")
TUMOR_COLORS <- c(NSCLC="#2166AC", CRC="#1A9641", HCC="#E6007E", PDAC="#D95F02", BRCA="#762A83")
RESP_COLORS  <- c(CR="#1A6B3C", PR="#2B83BA", SD="#F0A500", PD="#C0392B", NE="#999999")
RESP_LABELS  <- c(CR="Complete Response (CR)", PR="Partial Response (PR)",
                  SD="Stable Disease (SD)",    PD="Progressive Disease (PD)")
TMB_COLORS   <- c(Y="#C0392B", N="#2166AC")

adsl <- read_csv("ADSL.csv", show_col_types=FALSE)
adtr <- read_csv("ADTR.csv", show_col_types=FALSE)
adrs <- read_csv("ADRS.csv", show_col_types=FALSE)

rank_map <- c(CR=0, PR=1, SD=2, PD=3, NE=4)
best <- adrs %>%
  filter(PARAMCD == "OVRLRESP") %>%
  mutate(rank = rank_map[AVALC]) %>%
  mutate(rank = replace_na(rank, 4)) %>%
  arrange(rank) %>%
  group_by(USUBJID) %>%
  slice(1) %>%
  ungroup() %>%
  select(USUBJID, BESTRESP=AVALC)

meta <- adsl %>%
  filter(ARM == "TREATMENT") %>%
  select(USUBJID, TUMORTYPE, TMBHIGH, BESTRSPC) %>%
  left_join(best, by="USUBJID") %>%
  mutate(BESTRESP = replace_na(BESTRESP, "NE"))

trt_ids <- meta$USUBJID
adtr_trt <- adtr %>%
  filter(USUBJID %in% trt_ids) %>%
  mutate(across(c(ADTN, PCHG, AVAL), ~suppressWarnings(as.numeric(.))))

bl   <- adtr_trt %>% filter(AVISIT == "BASELINE")
post <- adtr_trt %>% filter(AVISIT != "BASELINE")

multi_visit_ids <- post %>%
  group_by(USUBJID) %>%
  summarise(n=n()) %>%
  filter(n >= 2) %>%
  pull(USUBJID)

valid_ids <- bl %>%
  filter(!is.na(AVAL), AVAL > 0, USUBJID %in% multi_visit_ids) %>%
  pull(USUBJID) %>%
  unique()

post_valid <- post %>%
  filter(USUBJID %in% valid_ids) %>%
  left_join(meta %>% select(USUBJID, BESTRSPC), by="USUBJID")

bl_rows <- tibble(USUBJID=valid_ids, months=0, PCHG=0)
spider <- post_valid %>%
  mutate(months = ADTN / DAY2MO) %>%
  filter(!is.na(PCHG), !is.na(months)) %>%
  select(USUBJID, months, PCHG) %>%
  bind_rows(bl_rows) %>%
  left_join(meta, by="USUBJID") %>%
  mutate(PCHG = pmin(PCHG, CLIP_TOP))

cat("Patients:", length(valid_ids), "\n")
cat("Spider rows:", nrow(spider), "\n")

ref_lines <- list(
  geom_hline(yintercept=0,     color="#222222", linewidth=1.0),
  geom_hline(yintercept=PR_TH, color="#666666", linewidth=0.6, linetype="dashed"),
  geom_hline(yintercept=PD_TH, color="#666666", linewidth=0.6, linetype="dashed"),
  annotate("text", x=XMAX+0.5, y=PD_TH+3, label="+20% PD",
           hjust=0, vjust=0, size=2.3, color="#666666", fontface="italic"),
  annotate("text", x=XMAX+0.5, y=PR_TH-3, label="\u221230% PR",
           hjust=0, vjust=1, size=2.3, color="#666666", fontface="italic")
)

base_theme <- function(show_y=FALSE) {
  t <- theme_minimal(base_size=9) +
    theme(
      panel.grid.major.x = element_blank(),
      panel.grid.minor   = element_blank(),
      panel.grid.major.y = element_line(color="#eeeeee", linewidth=0.4),
      axis.line          = element_blank(),
      panel.border       = element_rect(color="#bbbbbb", fill=NA, linewidth=0.5),
      plot.title         = element_text(size=7, color="#444444", hjust=0.5,
                                        margin=margin(b=2)),
      axis.title.x       = element_text(size=8.5, face="bold"),
      axis.title.y       = element_text(size=8.5, face="bold"),
      legend.position    = "none",
      plot.margin        = margin(5,5,5,5)
    )
  if (!show_y)
    t <- t + theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
                   axis.title.y=element_blank())
  t
}

make_spider_panel <- function(sub, color_col, color_map, show_y=FALSE,
                               panel_title="", title_color="black") {
  n   <- n_distinct(sub$USUBJID)
  rc  <- sub %>% distinct(USUBJID, .keep_all=TRUE) %>% count(BESTRESP)
  cr_ <- rc %>% filter(BESTRESP=="CR") %>% pull(n) %>% sum()
  pr_ <- rc %>% filter(BESTRESP=="PR") %>% pull(n) %>% sum()
  orr <- if (n>0) round((cr_+pr_)/n*100) else 0
  stats_lbl <- sprintf("N=%d   CR:%d  PR:%d  ORR=%d%%", n, cr_, pr_, orr)

  sub <- sub %>% arrange(USUBJID, months)
  sub$color <- color_map[as.character(sub[[color_col]])]
  sub$color[is.na(sub$color)] <- "#aaaaaa"

  ends <- sub %>% group_by(USUBJID) %>% slice_tail(n=1) %>% ungroup()
  ends$color <- color_map[as.character(ends[[color_col]])]
  ends$color[is.na(ends$color)] <- "#aaaaaa"

  p <- ggplot() +
    ref_lines +
    geom_line(data=sub, aes(x=months, y=PCHG, group=USUBJID, color=color),
              linewidth=0.75, alpha=0.70, lineend="butt", linejoin="miter") +
    geom_point(data=ends, aes(x=months, y=PCHG, color=color),
               size=1.5, alpha=0.85) +
    scale_color_identity() +
    scale_x_continuous(limits=c(-0.5, XMAX+3),
                       breaks=seq(0, XMAX, 6),
                       labels=seq(0, XMAX, 6),
                       expand=c(0,0)) +
    scale_y_continuous(limits=c(YMIN, YMAX),
                       breaks=seq(-100, 120, 20),
                       labels=paste0(seq(-100,120,20),"%"),
                       expand=c(0,0)) +
    labs(x="Time from Treatment Start (months)",
         y="Change from Baseline in SLD (%)") +
    ggtitle(stats_lbl) +
    base_theme(show_y=show_y) +
    annotate("text", x=0.5*(XMAX+3)/2, y=YMAX*0.93,
             label=panel_title, hjust=0.5, vjust=1,
             size=4.5, fontface="bold", color=title_color)
  p
}

panels_A <- lapply(seq_along(TUMOR_ORDER), function(i) {
  tumor <- TUMOR_ORDER[i]
  sub   <- spider %>% filter(TUMORTYPE == tumor)
  make_spider_panel(sub, color_col="BESTRESP", color_map=RESP_COLORS,
                    show_y=(i==1), panel_title=tumor,
                    title_color=TUMOR_COLORS[tumor])
})

legend_patches_A <- tibble(
  resp  = c("CR","PR","SD","PD"),
  label = RESP_LABELS[c("CR","PR","SD","PD")],
  color = RESP_COLORS[c("CR","PR","SD","PD")]
)
legend_A <- ggplot(legend_patches_A, aes(xmin=0,xmax=1,ymin=0,ymax=1,fill=label)) +
  geom_rect() +
  scale_fill_manual(values=setNames(legend_patches_A$color, legend_patches_A$label),
                    name="Best Overall Response") +
  theme_void() +
  theme(legend.position="bottom",
        legend.title=element_text(size=8, face="bold"),
        legend.text=element_text(size=8),
        legend.key.size=unit(0.45,"cm"))
leg_A <- cowplot::get_legend(legend_A)

combined_A <- wrap_plots(panels_A, nrow=1) +
  plot_annotation(
    title    = "Spider Plot \u2014 Percent Change from Baseline in SLD by Tumor Type   \u00b7   ONCVIZ-001  \u00b7  Treatment Arm",
    subtitle = paste0("RECIST 1.1  \u00b7  Lines clipped at +", CLIP_TOP, "%  \u00b7  ", CUTOFF_LBL),
    theme    = theme(
      plot.title    = element_text(size=12, face="bold",  hjust=0.5),
      plot.subtitle = element_text(size=8.5, hjust=0.5, face="italic", color="#555555")
    )
  )

png(OUTPUT_A, width=22, height=7.8, units="in", res=DPI, bg="white")
print(combined_A)
grid::grid.draw(leg_A)
dev.off()
cat("\u2713", OUTPUT_A, "\n")

panels_B <- lapply(seq_along(c("Y","N")), function(i) {
  tval <- c("Y","N")[i]
  tlbl <- c("TMB-High","TMB-Low")[i]
  sub  <- spider %>% filter(TMBHIGH == tval)
  make_spider_panel(sub, color_col="TUMORTYPE", color_map=TUMOR_COLORS,
                    show_y=(i==1), panel_title=tlbl,
                    title_color=TMB_COLORS[tval])
})

legend_patches_B <- tibble(
  tumor = TUMOR_ORDER,
  color = TUMOR_COLORS[TUMOR_ORDER]
)
legend_B <- ggplot(legend_patches_B, aes(xmin=0,xmax=1,ymin=0,ymax=1,fill=tumor)) +
  geom_rect() +
  scale_fill_manual(values=setNames(legend_patches_B$color, legend_patches_B$tumor),
                    name="Tumor Type") +
  theme_void() +
  theme(legend.position="bottom",
        legend.title=element_text(size=8, face="bold"),
        legend.text=element_text(size=8),
        legend.key.size=unit(0.45,"cm"))
leg_B <- cowplot::get_legend(legend_B)

combined_B <- wrap_plots(panels_B, nrow=1) +
  plot_annotation(
    title    = "Spider Plot \u2014 Percent Change from Baseline in SLD by TMB Status   \u00b7   ONCVIZ-001  \u00b7  Treatment Arm",
    subtitle = paste0("RECIST 1.1  \u00b7  Lines clipped at +", CLIP_TOP, "%  \u00b7  ", CUTOFF_LBL),
    theme    = theme(
      plot.title    = element_text(size=12, face="bold",  hjust=0.5),
      plot.subtitle = element_text(size=8.5, hjust=0.5, face="italic", color="#555555")
    )
  )

png(OUTPUT_B, width=15, height=7.8, units="in", res=DPI, bg="white")
print(combined_B)
grid::grid.draw(leg_B)
dev.off()
cat("\u2713", OUTPUT_B, "\n")
cat("Done.\n")
