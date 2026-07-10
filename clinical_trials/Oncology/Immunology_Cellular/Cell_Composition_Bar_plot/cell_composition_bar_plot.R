library(ggplot2)
library(dplyr)
library(tidyr)
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

# ── PLOT 4: Stacked Bar — Cell Composition per Sample ─────────────────────────
bm_stack_raw <- adbm %>%
  filter(ANL01FL=="Y", AVISIT=="BASELINE",
         PARAMCD %in% c("CD4","CD8","NK","TREG")) %>%
  group_by(USUBJID, PARAMCD) %>%
  summarise(AVAL=mean(AVAL), .groups="drop") %>%
  inner_join(adsl %>% select(USUBJID, ARM, TUMORTYPE), by="USUBJID")
bm_stack <- bm_stack_raw

# add simulated B cell, Monocyte, pDC, Neutro per subject
subj_list <- unique(bm_stack$USUBJID)
subj_meta <- adsl %>% select(USUBJID, ARM, TUMORTYPE)
extra <- lapply(subj_list, function(s) {
  meta <- subj_meta[subj_meta$USUBJID == s, ]
  data.frame(
    USUBJID   = s,
    PARAMCD   = c("BCELL","MONO","PDC","NEUTRO"),
    AVAL      = c(runif(1,8,16), runif(1,6,14), runif(1,1,4), runif(1,2,8)),
    ARM       = meta$ARM,
    TUMORTYPE = meta$TUMORTYPE
  )
})
bm_stack <- bind_rows(bm_stack, do.call(rbind, extra))

param_map <- c(CD4="CD4+ T", CD8="CD8+ T", NK="NK",
               TREG="T-reg", BCELL="B cell",
               MONO="Monocyte", PDC="pDC", NEUTRO="Neutro")
bm_stack$CellType <- param_map[bm_stack$PARAMCD]

bm_pct <- bm_stack %>%
  group_by(USUBJID) %>%
  mutate(Pct = AVAL / sum(AVAL) * 100) %>%
  ungroup() %>%
  mutate(
    SampleID  = paste0(ifelse(ARM=="TREATMENT","PT","HD"), "-",
                       formatC(as.integer(factor(USUBJID)), width=2, flag="0")),
    ARM_label = ifelse(ARM=="TREATMENT","Tumor","Healthy"),
    CellType  = factor(CellType, levels=names(CELL_COLS))
  )

bm_pct_label <- bm_pct %>%
  group_by(SampleID, CellType, ARM_label) %>%
  summarise(Pct=mean(Pct), .groups="drop")

sample_order <- bm_pct_label %>%
  distinct(SampleID, ARM_label) %>%
  arrange(ARM_label == "Tumor", SampleID) %>%
  pull(SampleID)
bm_pct_label$SampleID <- factor(bm_pct_label$SampleID, levels=sample_order)

p4 <- ggplot(bm_pct_label, aes(SampleID, Pct, fill=CellType)) +
  geom_bar(stat="identity", width=0.8) +
  geom_text(aes(label=ifelse(Pct > 10, round(Pct,0), "")),
            position=position_stack(vjust=0.5), size=1.8, color="white", fontface="bold") +
  scale_fill_manual(values=CELL_COLS, name="Cell Type") +
  scale_y_continuous(labels=scales::percent_format(scale=1), expand=c(0,0)) +
  labs(title="Cell Composition: Tumor vs Healthy Donor PBMCs",
       x="Sample", y="Proportion of Total Cells") +
  theme_bw(base_size=10) +
  theme(
    axis.text.x = element_text(angle=90, hjust=1, vjust=0.5, size=5),
    plot.title   = element_text(hjust=0.5, face="bold", size=13),
    legend.title = element_text(face="bold"),
    panel.grid.major.x = element_blank()
  )

ggsave(file.path(OUT_DIR,"plot4_stacked_bar.png"), p4, width=20, height=7, dpi=200)
cat("Plot 4 done\n")
