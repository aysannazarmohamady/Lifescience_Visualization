library(ggplot2); library(ggrepel); library(dplyr)

set.seed(9)
GENES <- c("ARID1A","PIK3CA","EGFR","MET","BRAF","CDKN2A","RET","PTEN",
           "KRAS","RB1","TP53","NF1","SMAD4","KEAP1","STK11")
CHR_MAP   <- c(ARID1A="chr1",PIK3CA="chr3",EGFR="chr7",MET="chr7",BRAF="chr7",
               CDKN2A="chr9",RET="chr10",PTEN="chr10",KRAS="chr12",RB1="chr13",
               TP53="chr17",NF1="chr17",SMAD4="chr18",KEAP1="chr19",STK11="chr19")
CHR_ORDER <- c("chr1","chr3","chr7","chr9","chr10","chr12","chr13","chr17","chr18","chr19")

vars <- do.call(rbind, lapply(1:138, function(i) {
  g <- sample(GENES,1); data.frame(HUGO_SYMBOL=g, CHR=CHR_MAP[g])
}))

gene_stats <- data.frame(
  HUGO_SYMBOL = GENES,
  pval_tmb = runif(15,.001,.5),
  pval_msi = runif(15,.001,.5)
) %>% mutate(
  pval_tmb = ifelse(HUGO_SYMBOL %in% c("TP53","KRAS"), pval_tmb*.05, pval_tmb),
  pval_msi = ifelse(HUGO_SYMBOL %in% c("ARID1A","PTEN"), pval_msi*.05, pval_msi))

miami_tmb <- vars %>% left_join(gene_stats, by="HUGO_SYMBOL") %>%
  mutate(CHR=factor(CHR,levels=CHR_ORDER),
         xi=as.integer(CHR)+runif(n(),-.28,.28),
         log10p_dir=-log10(pval_tmb), panel="TMB-High")

miami_msi <- vars %>% left_join(gene_stats, by="HUGO_SYMBOL") %>%
  mutate(CHR=factor(CHR,levels=CHR_ORDER),
         xi=as.integer(CHR)+runif(n(),-.28,.28),
         log10p_dir=log10(pval_msi), panel="MSI-High")

all_df <- bind_rows(miami_tmb, miami_msi)

lbl_df <- all_df %>% group_by(HUGO_SYMBOL,panel) %>%
  slice_max(abs(log10p_dir),n=1) %>% filter(abs(log10p_dir)>.5)

p <- ggplot(all_df, aes(xi, log10p_dir)) +
  geom_hline(data=data.frame(y=-log10(.05),panel="TMB-High"),
             aes(yintercept=y), color="#e8a020",linetype="dashed") +
  geom_hline(data=data.frame(y=log10(.05),panel="MSI-High"),
             aes(yintercept=y), color="#e8a020",linetype="dashed") +
  geom_point(aes(color=panel), size=2, alpha=.75) +
  geom_text_repel(data=lbl_df, aes(label=HUGO_SYMBOL),
                  size=2.8, segment.color="#cccccc", box.padding=.4) +
  scale_color_manual(values=c("TMB-High"="#2c5f8a","MSI-High"="#c0392b"),
                     name="Association") +
  facet_wrap(~panel, ncol=1, scales="free_y", strip.position="left") +
  scale_x_continuous(breaks=1:length(CHR_ORDER), labels=CHR_ORDER, name="Chromosome") +
  scale_y_continuous(name="−log₁₀(p)", labels=abs) +
  theme_classic(base_size=10) +
  theme(strip.placement="outside") +
  labs(title="Miami Plot — TMB-High (top) vs MSI-High (bottom)",
       subtitle="n = 138 mutations · 15 genes · 80 patients")

ggsave("plots/09_miami.png", p, width=16, height=10, dpi=180, bg="white")
cat("09 done\n")
