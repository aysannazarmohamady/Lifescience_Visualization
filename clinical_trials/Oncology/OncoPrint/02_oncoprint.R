library(ggplot2); library(dplyr); library(tidyr); library(patchwork)

set.seed(2)
GENES <- c("TP53","KRAS","CDKN2A","PTEN","SMAD4","PIK3CA",
           "ARID1A","RB1","NF1","MET","BRAF","KEAP1")
VC_COLORS <- c(Missense_Mutation="#2c5f8a", Nonsense_Mutation="#c0392b",
               Frame_Shift_Del="#e67e22", Frame_Shift_Ins="#f1c40f",
               Splice_Site="#8e44ad", In_Frame_Del="#27ae60",
               Amplification="#7b0000")
BOR_COLORS <- c(CR="#1a3f8f", PR="#4f8fd4", SD="#e8a020", PD="#c0392b")
TT_COLORS  <- c(NSCLC="#2c5f8a", BRCA="#c0392b", HCC="#27ae60",
                CRC="#e67e22",   PDAC="#8e44ad")

pts <- paste0("PT",sprintf("%03d",1:80))
adsl <- data.frame(
  USUBJID  = pts,
  TMB      = round(runif(80, 1, 45)),
  BESTRSPC = sample(c("CR","PR","SD","PD"), 80, TRUE, prob=c(.1,.35,.35,.2)),
  TUMORTYPE= sample(names(TT_COLORS), 80, TRUE)
)

mut_list <- lapply(pts, function(p) {
  g <- sample(GENES, sample(1:6,1))
  data.frame(USUBJID=p, HUGO_SYMBOL=g,
             VARIANT_CLASS=sample(names(VC_COLORS), length(g), TRUE))
})
mat <- do.call(rbind, mut_list) %>%
  group_by(USUBJID, HUGO_SYMBOL) %>% slice(1) %>% ungroup()

pt_order <- mat %>% count(USUBJID, sort=TRUE) %>% pull(USUBJID)
pt_order_full <- c(pt_order, setdiff(pts, pt_order))

gene_pct <- mat %>% group_by(HUGO_SYMBOL) %>%
  summarise(pct=n_distinct(USUBJID)/80*100) %>%
  mutate(HUGO_SYMBOL=factor(HUGO_SYMBOL, levels=rev(GENES)))

tmb_df <- adsl %>%
  mutate(USUBJID=factor(USUBJID, levels=pt_order_full))
mat2 <- mat %>%
  mutate(USUBJID=factor(USUBJID, levels=pt_order_full),
         HUGO_SYMBOL=factor(HUGO_SYMBOL, levels=rev(GENES)))

th <- theme_classic(base_size=8) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.x=element_blank())

p_tmb <- ggplot(tmb_df, aes(USUBJID, TMB)) +
  geom_col(fill="#555555") + th + labs(y="TMB")

p_bor <- ggplot(tmb_df, aes(USUBJID, 1, fill=BESTRSPC)) +
  geom_tile(color="white") +
  scale_fill_manual(values=BOR_COLORS) + theme_void()

p_tt <- ggplot(tmb_df, aes(USUBJID, 1, fill=TUMORTYPE)) +
  geom_tile(color="white") +
  scale_fill_manual(values=TT_COLORS) + theme_void()

p_heat <- ggplot(mat2, aes(USUBJID, HUGO_SYMBOL, fill=VARIANT_CLASS)) +
  geom_tile(color="white", height=0.85, width=0.95) +
  scale_fill_manual(values=VC_COLORS, na.value="grey95") +
  theme_classic(base_size=8) + th + labs(y=NULL)

p_pct <- ggplot(gene_pct, aes(pct, HUGO_SYMBOL)) +
  geom_col(fill="#2c5f8a", width=0.6) +
  geom_text(aes(label=paste0(round(pct),"%")), hjust=-0.1, size=2.8) +
  scale_x_continuous(limits=c(0,70), name="Altered (%)") +
  theme_classic(base_size=8) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank())

p <- (p_tmb / p_bor / p_tt / p_heat) | p_pct +
  plot_layout(widths=c(19,1)) +
  plot_annotation(title="OncoPrint — Somatic Mutation Landscape, Top 12 Genes",
                  subtitle="n = 80 patients")

ggsave("plots/02_oncoprint.png", p, width=20, height=10, dpi=180, bg="white")
cat("02 done\n")
