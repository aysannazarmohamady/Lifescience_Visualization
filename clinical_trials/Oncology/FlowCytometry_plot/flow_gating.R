library(ggplot2)
library(dplyr)
library(tidyr)
library(gridExtra)
library(ggpubr)

set.seed(42)

DATA_DIR <- "data/v1"
OUT_DIR  <- "plots"
dir.create(OUT_DIR, showWarnings = FALSE)

adsl <- read.csv(file.path(DATA_DIR, "ADSL.csv"), stringsAsFactors = FALSE)
adbm <- read.csv(file.path(DATA_DIR, "ADBM.csv"), stringsAsFactors = FALSE)

# ── PLOT 1: Flow Cytometry Hierarchical Gating ────────────────────────────────
bm_bl <- adbm %>%
  filter(ANL01FL == "Y", AVISIT == "BASELINE") %>%
  select(USUBJID, PARAMCD, AVAL) %>%
  pivot_wider(names_from = PARAMCD, values_from = AVAL, values_fn = mean) %>%
  inner_join(adsl %>% select(USUBJID, ARM, TUMORTYPE), by = "USUBJID") %>%
  filter(!is.na(CD4), !is.na(CD8))

n_sub  <- nrow(bm_bl)

sim_fsc_all <- c(
  rnorm(n_sub * 4, 5.2, 0.9),
  rnorm(n_sub * 3, 2.0, 0.7),
  rnorm(n_sub * 3, 3.5, 1.2),
  rnorm(n_sub * 4, 6.5, 0.8),
  rnorm(n_sub * 4, 1.5, 0.5)
)
sim_ssc_all <- c(
  rnorm(n_sub * 4, 5.5, 0.7),
  rnorm(n_sub * 3, 2.5, 0.6),
  rnorm(n_sub * 3, 4.0, 1.0),
  rnorm(n_sub * 4, 3.5, 1.0),
  rnorm(n_sub * 4, 1.5, 0.4)
)
cell_label_all <- rep(
  c("CD4+ T","CD8+ T","B cell","Monocyte","Debris"),
  c(n_sub*4, n_sub*3, n_sub*3, n_sub*4, n_sub*4)
)

df_fsc <- data.frame(
  FSC = pmin(pmax(sim_fsc_all, 0.2), 9),
  SSC = pmin(pmax(sim_ssc_all, 0.2), 9),
  CellType = cell_label_all
)

lymph_idx <- df_fsc$FSC >= 3 & df_fsc$FSC <= 7.5 &
             df_fsc$SSC >= 3.5 & df_fsc$SSC <= 7.5

df_lymph <- df_fsc[lymph_idx, ]
df_lymph$CD3  <- ifelse(df_lymph$CellType %in% c("CD4+ T","CD8+ T"),
                        rnorm(sum(lymph_idx), 7, 0.4),
                        rnorm(sum(lymph_idx), 1.2, 0.3))
df_lymph$CD19 <- ifelse(df_lymph$CellType == "B cell",
                        rnorm(sum(lymph_idx), 7, 0.4),
                        rnorm(sum(lymph_idx), 1.0, 0.3))

t_idx  <- df_lymph$CellType %in% c("CD4+ T","CD8+ T")
df_t   <- df_lymph[t_idx, ]
df_t$CD4val <- ifelse(df_t$CellType == "CD4+ T",
                      rnorm(sum(t_idx), 7, 0.45),
                      rnorm(sum(t_idx), 1.0, 0.3))
df_t$CD8val <- ifelse(df_t$CellType == "CD8+ T",
                      rnorm(sum(t_idx), 7, 0.45),
                      rnorm(sum(t_idx), 1.0, 0.3))

cd4_pct  <- round(mean(bm_bl$CD4), 1)
cd8_pct  <- round(mean(bm_bl$CD8), 1)
treg_pct <- round(mean(bm_bl$TREG, na.rm = TRUE), 1)

col_map <- c("CD4+ T"="#2166AC","CD8+ T"="#D6604D","B cell"="#F4A582",
             "Monocyte"="#762A83","T-reg"="#1B7837","Debris"="#BBBBBB","NK"="#4DAC26")

p1a <- ggplot(df_fsc, aes(FSC, SSC, color = CellType)) +
  geom_point(size = 0.3, alpha = 0.4) +
  annotate("rect", xmin=3, xmax=7.5, ymin=3.5, ymax=7.5,
           fill="red", alpha=0.08, color="red", linewidth=0.8) +
  annotate("label", x=5.8, y=7.1, label="Lymphocyte\nGate",
           color="red", size=3, fontface="bold", fill="white") +
  scale_color_manual(values=col_map) +
  labs(title="A  FSC vs SSC — All Events",
       x="FSC-A (Forward Scatter)", y="SSC-A (Side Scatter)") +
  xlim(0,9) + ylim(0,9) +
  theme_bw(base_size=9) +
  theme(legend.position="none", plot.title=element_text(size=8,face="bold"))

p1b <- ggplot(df_lymph, aes(CD3, CD19, color=CellType)) +
  geom_point(size=0.4, alpha=0.5) +
  geom_vline(xintercept=4.2, linetype="dashed", color="grey40") +
  geom_hline(yintercept=3.5, linetype="dashed", color="grey40") +
  annotate("text", x=2.0, y=7.0, label="CD3\u207BCD19\u207A\n(B cells)",
           color="#F4A582", size=2.8, fontface="bold") +
  annotate("text", x=7.0, y=7.0, label="CD3\u207ACD19\u207B\n(T cells)",
           color="#2166AC", size=2.8, fontface="bold") +
  annotate("text", x=1.5, y=1.2, label="DN", color="grey40", size=2.5) +
  scale_color_manual(values=col_map) +
  labs(title="B  CD3 vs CD19 — Lymphocyte Gate",
       x="CD3 (FITC)", y="CD19 (PE)") +
  xlim(0,9) + ylim(0,9) +
  theme_bw(base_size=9) +
  theme(legend.position="none", plot.title=element_text(size=8,face="bold"))

p1c <- ggplot(df_t, aes(CD4val, CD8val, color=CellType)) +
  geom_point(size=0.4, alpha=0.5) +
  geom_vline(xintercept=3.5, linetype="dashed", color="grey40") +
  geom_hline(yintercept=3.5, linetype="dashed", color="grey40") +
  annotate("text", x=7, y=8.2, label=paste0("CD4\u207A: ", cd4_pct, "%"),
           color="#2166AC", size=3, fontface="bold") +
  annotate("text", x=7, y=3.5, label=paste0("CD8\u207A: ", cd8_pct, "%"),
           color="#D6604D", size=3, fontface="bold") +
  annotate("text", x=1.5, y=1.5, label=paste0("T-reg: ", treg_pct, "%"),
           color="#1B7837", size=3, fontface="bold") +
  scale_color_manual(values=col_map) +
  labs(title="C  CD4 vs CD8 — T-cell Gate",
       x="CD4 (BV421)", y="CD8 (APC)") +
  xlim(0,9) + ylim(0,9) +
  theme_bw(base_size=9) +
  theme(legend.position="none", plot.title=element_text(size=8,face="bold"))

legend_df <- data.frame(
  CellType = names(col_map)[names(col_map) != "Debris"],
  x=1, y=seq_along(names(col_map)[names(col_map)!="Debris"])
)
p_leg <- ggplot(legend_df, aes(x, y, color=CellType)) +
  geom_point(size=3) +
  scale_color_manual(values=col_map, name="Cell Type") +
  theme_void() +
  theme(legend.position="bottom",
        legend.text=element_text(size=7),
        legend.title=element_text(size=7,face="bold"))
leg <- ggpubr::get_legend(p_leg)

p1 <- grid.arrange(
  arrangeGrob(p1a, p1b, p1c, nrow=1),
  leg, nrow=2, heights=c(10,1),
  top=grid::textGrob("Flow Cytometry — Hierarchical Gating Strategy",
                     gp=grid::gpar(fontsize=11, fontface="bold"))
)
ggsave(file.path(OUT_DIR,"plot1_flow_gating.png"), p1, width=14, height=5, dpi=150)
cat("Plot 1 done\n")
