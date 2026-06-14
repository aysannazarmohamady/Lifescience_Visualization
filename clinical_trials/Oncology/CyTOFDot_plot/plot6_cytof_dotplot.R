library(ggplot2)
library(scales)

set.seed(42)

OUT_DIR <- "plots"
dir.create(OUT_DIR, showWarnings = FALSE)

CELL_COLS <- c(
  "CD4+ T"  = "#2166AC",
  "CD8+ T"  = "#D6604D",
  "NK"      = "#4DAC26",
  "B cell"  = "#F4A582",
  "Monocyte"= "#762A83",
  "T-reg"   = "#1B7837",
  "pDC"     = "#FEE08B",
  "Neutro"  = "#8C510A"
)

# ── PLOT 6: CyTOF Dot Plot ────────────────────────────────────────────────────
markers_dot <- c("CD3","CD4","CD8a","CD56","CD19","CD14","CD16",
                 "CD25","CD127","FoxP3","GZMB","PD-1","TIM-3",
                 "LAG-3","IFN-\u03b3","TNF-\u03b1","IL-2","Ki-67")

cell_types_dot <- names(CELL_COLS)

expr_profiles <- list(
  "CD3"    = c(.95,.95,.95,.05,.05,.05,.05,.05),
  "CD4"    = c(.85,.05,.80,.05,.05,.05,.05,.05),
  "CD8a"   = c(.05,.90,.05,.05,.05,.05,.05,.05),
  "CD56"   = c(.05,.05,.05,.90,.05,.05,.05,.10),
  "CD19"   = c(.05,.05,.05,.05,.95,.10,.05,.05),
  "CD14"   = c(.05,.05,.05,.05,.05,.90,.05,.10),
  "CD16"   = c(.05,.05,.05,.55,.05,.50,.05,.80),
  "CD25"   = c(.10,.10,.95,.05,.05,.05,.05,.05),
  "CD127"  = c(.65,.55,.05,.50,.55,.50,.50,.45),
  "FoxP3"  = c(.05,.05,.92,.05,.05,.05,.05,.05),
  "GZMB"   = c(.05,.65,.05,.55,.05,.10,.05,.10),
  "PD-1"   = c(.35,.45,.25,.20,.05,.05,.05,.05),
  "TIM-3"  = c(.10,.35,.15,.10,.05,.05,.05,.05),
  "LAG-3"  = c(.10,.30,.10,.05,.05,.05,.05,.05),
  "IFN-\u03b3" = c(.30,.60,.10,.55,.05,.45,.10,.10),
  "TNF-\u03b1"  = c(.30,.50,.10,.45,.10,.50,.15,.15),
  "IL-2"   = c(.30,.40,.10,.10,.05,.10,.05,.05),
  "Ki-67"  = c(.10,.15,.10,.10,.10,.15,.10,.10)
)

dot_df <- do.call(rbind, lapply(names(expr_profiles), function(mk) {
  vals <- expr_profiles[[mk]]
  data.frame(
    Marker   = mk,
    CellType = cell_types_dot,
    MeanExpr = vals,
    PctExpr  = pmin(vals * 110 + runif(length(vals), 0, 8), 100)
  )
}))

dot_df$Marker   <- factor(dot_df$Marker,   levels=rev(markers_dot))
dot_df$CellType <- factor(dot_df$CellType, levels=cell_types_dot)

p6 <- ggplot(dot_df, aes(CellType, Marker, size=PctExpr, fill=MeanExpr)) +
  geom_point(shape=21, stroke=0.3, color="grey60") +
  scale_fill_gradientn(
    colors = c("white","#FEE08B","#FC8D59","#D73027","#67001F"),
    values = scales::rescale(c(0, 0.25, 0.5, 0.75, 1)),
    limits = c(0, 1),
    name   = "Mean Expression\n(normalized 0-1)"
  ) +
  scale_size_continuous(
    range  = c(1, 12),
    limits = c(0, 100),
    breaks = c(25, 50, 75, 100),
    name   = "% Cells\nExpressing"
  ) +
  labs(
    title    = "CyTOF Dot Plot - Marker Expression Profile",
    subtitle = "Across Immune Cell Populations",
    x="Cell Type", y="Marker"
  ) +
  theme_bw(base_size=10) +
  theme(
    plot.title    = element_text(hjust=0.5, face="bold"),
    plot.subtitle = element_text(hjust=0.5),
    axis.text.x   = element_text(angle=30, hjust=1),
    panel.grid    = element_line(color="grey92"),
    legend.title  = element_text(size=8, face="bold"),
    legend.text   = element_text(size=7)
  )

ggsave(file.path(OUT_DIR,"plot6_cytof_dotplot.png"), p6, width=11, height=9, dpi=150)
cat("Plot 6 done\n")
