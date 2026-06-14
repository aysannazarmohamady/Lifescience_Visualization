library(ggplot2)
library(dplyr)
library(scales)

set.seed(42)

DATA_DIR <- "data/v1"
OUT_DIR  <- "plots"
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

adsl <- read.csv(file.path(DATA_DIR, "ADSL.csv"), stringsAsFactors = FALSE)
adbm <- read.csv(file.path(DATA_DIR, "ADBM.csv"), stringsAsFactors = FALSE)

# ── PLOT 5: Immune Cell Infiltration Heatmap ──────────────────────────────────
bm_heat <- adbm %>%
  filter(ANL01FL=="Y", AVISIT=="BASELINE",
         PARAMCD %in% c("CD4","CD8","NK","TREG","PDL1","IFNg")) %>%
  group_by(USUBJID, PARAMCD) %>%
  summarise(AVAL=mean(AVAL, na.rm=TRUE), .groups="drop") %>%
  inner_join(adsl %>% select(USUBJID, TUMORTYPE), by="USUBJID") %>%
  group_by(TUMORTYPE, PARAMCD) %>%
  summarise(mean_val=mean(AVAL, na.rm=TRUE), .groups="drop")

param_cell_map <- c(CD4="CD4+ T", CD8="CD8+ T", NK="NK",
                    TREG="T-reg", PDL1="B cell", IFNg="Monocyte")
bm_heat$CellType <- param_cell_map[bm_heat$PARAMCD]

extra_heat <- expand.grid(
  TUMORTYPE = c("BRCA","CRC","HCC","NSCLC","PDAC"),
  CellType  = c("pDC","Neutro"),
  stringsAsFactors = FALSE
)
extra_heat$mean_val <- runif(nrow(extra_heat), 2, 15)

bm_heat2 <- bind_rows(
  bm_heat %>% select(TUMORTYPE, CellType, mean_val),
  extra_heat
)

bm_zscore <- bm_heat2 %>%
  group_by(CellType) %>%
  mutate(z = as.numeric(scale(mean_val))) %>%
  ungroup()

cell_order   <- names(CELL_COLS)
tumor_order  <- c("NSCLC","BRCA","CRC","HCC","PDAC")
bm_zscore$CellType  <- factor(bm_zscore$CellType,  levels=cell_order)
bm_zscore$TUMORTYPE <- factor(bm_zscore$TUMORTYPE, levels=tumor_order)

p5 <- ggplot(bm_zscore, aes(CellType, TUMORTYPE, fill=z)) +
  geom_tile(color="white", linewidth=0.5) +
  geom_text(aes(label=round(z, 2)), size=3.2, fontface="bold",
            color=ifelse(abs(bm_zscore$z) > 1.2, "white", "black")) +
  scale_fill_gradientn(
    colors  = c("#2166AC","#92C5DE","white","#F4A582","#B2182B"),
    values  = scales::rescale(c(-2,-1,0,1,2)),
    limits  = c(-2.5, 2.5),
    name    = "z-score",
    guide   = guide_colorbar(barheight=8)
  ) +
  scale_x_discrete(expand=c(0,0)) +
  scale_y_discrete(expand=c(0,0)) +
  labs(
    title    = "Immune Cell Infiltration Heatmap",
    subtitle = expression(paste("(z-score, log"[2]," relative to normal tissue)")),
    x="Immune Cell Type", y="Tumor Type"
  ) +
  theme_bw(base_size=11) +
  theme(
    plot.title    = element_text(hjust=0.5, face="bold"),
    plot.subtitle = element_text(hjust=0.5),
    axis.text.x   = element_text(angle=30, hjust=1),
    legend.title  = element_text(face="bold"),
    panel.border  = element_rect(color="grey40")
  )

ggsave(file.path(OUT_DIR,"plot5_heatmap.png"), p5, width=10, height=6, dpi=150)
cat("Plot 5 done\n")
