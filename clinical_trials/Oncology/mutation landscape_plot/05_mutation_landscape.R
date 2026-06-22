library(ggplot2); library(dplyr); library(tidyr); library(patchwork)

set.seed(5)
GENES15 <- c("TP53","KRAS","PTEN","ARID1A","CDKN2A","PIK3CA",
             "KEAP1","SMAD4","RB1","NF1","MET","EGFR","BRAF","RET","STK11")
TT <- c(NSCLC=25,BRCA=13,HCC=12,CRC=16,PDAC=14)
TT_COLORS <- c(NSCLC="#2c5f8a",BRCA="#c0392b",HCC="#27ae60",CRC="#e67e22",PDAC="#8e44ad")

pts <- data.frame(
  USUBJID   = paste0("PT",sprintf("%03d",1:80)),
  TUMORTYPE = rep(names(TT), TT)
)
mut_list <- lapply(1:80, function(i)
  data.frame(USUBJID=pts$USUBJID[i], TUMORTYPE=pts$TUMORTYPE[i],
             HUGO_SYMBOL=sample(GENES15, sample(2:7,1))))
point <- do.call(rbind, mut_list) %>% distinct(USUBJID, HUGO_SYMBOL, .keep_all=TRUE)

mut_tt <- point %>% group_by(HUGO_SYMBOL, TUMORTYPE) %>%
  summarise(n=n_distinct(USUBJID), .groups="drop") %>%
  mutate(pct=n/TT[TUMORTYPE]*100)

pan_pct <- point %>% group_by(HUGO_SYMBOL) %>%
  summarise(pan=n_distinct(USUBJID)/80*100, .groups="drop")
gene_ord <- pan_pct %>% arrange(pan) %>% pull(HUGO_SYMBOL)

mut_tt  <- mut_tt  %>% mutate(HUGO_SYMBOL=factor(HUGO_SYMBOL,levels=gene_ord),
                               TUMORTYPE=factor(TUMORTYPE,levels=names(TT_COLORS)))
pan_pct <- pan_pct %>% mutate(HUGO_SYMBOL=factor(HUGO_SYMBOL,levels=gene_ord))

p5a <- ggplot(mut_tt, aes(pct, HUGO_SYMBOL, fill=TUMORTYPE)) +
  geom_col(position=position_dodge(.8), width=.7, alpha=.88) +
  scale_fill_manual(values=TT_COLORS, name="Tumor Type") +
  scale_x_continuous(name="Mutation Prevalence (%)", limits=c(0,110)) +
  theme_classic(base_size=10) + theme(axis.title.y=element_blank()) +
  labs(title="By Tumor Type")

p5b <- ggplot(pan_pct, aes(pan, HUGO_SYMBOL)) +
  geom_col(fill="#444444", width=.6) +
  geom_text(aes(label=paste0(round(pan),"%")), hjust=-0.1, size=3) +
  scale_x_continuous(name="Pan-Cohort (%)", limits=c(0,65)) +
  theme_classic(base_size=10) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.title.y=element_blank()) +
  labs(title="Overall")

p <- p5a + p5b + plot_layout(widths=c(4,1)) +
  plot_annotation(title="Mutation Landscape — Somatic Alteration Frequencies",
                  subtitle="n = 80 patients · 15 genes")

ggsave("plots/05_mutation_landscape.png", p, width=15, height=10, dpi=180, bg="white")
cat("05 done\n")
