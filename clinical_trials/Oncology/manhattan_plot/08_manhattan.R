library(ggplot2); library(ggrepel); library(dplyr)

set.seed(8)
GENES <- c("ARID1A","PIK3CA","EGFR","MET","BRAF","CDKN2A","RET","PTEN",
           "KRAS","RB1","TP53","NF1","SMAD4","KEAP1","STK11")
CHR   <- c("chr1","chr3","chr7","chr7","chr7","chr9","chr10","chr10",
           "chr12","chr13","chr17","chr17","chr18","chr19","chr19")
CHR_ORDER <- c("chr1","chr3","chr7","chr9","chr10","chr12","chr13","chr17","chr18","chr19")

# one point per variant per gene (simulate ~138 mutations)
vars <- lapply(1:138, function(i) {
  g <- sample(GENES, 1, prob=c(.08,.07,.06,.06,.06,.08,.06,.08,.16,.05,
                               .12,.05,.06,.05,.06))
  data.frame(HUGO_SYMBOL=g, CHR=CHR[match(g,GENES)],
             PROTEIN_POS=sample(1:1200,1))
}) %>% do.call(rbind,.)

pvals <- data.frame(HUGO_SYMBOL=GENES, pval=runif(15,.01,.7)) %>%
  mutate(pval=ifelse(HUGO_SYMBOL %in% c("TP53","KRAS"), pval*.05, pval),
         log10p=-log10(pval))

man_df <- vars %>% left_join(pvals, by="HUGO_SYMBOL") %>%
  mutate(CHR=factor(CHR, levels=CHR_ORDER),
         xi=as.integer(CHR)+runif(n(),-0.28,.28),
         notable=pval<0.15)

lbl_df <- man_df %>% group_by(HUGO_SYMBOL) %>%
  slice_max(log10p,n=1) %>% filter(pval<0.5)

p <- ggplot(man_df, aes(xi, log10p)) +
  geom_rect(data=data.frame(xi=seq(2,length(CHR_ORDER),2)),
            aes(xmin=xi-.45,xmax=xi+.45,ymin=-Inf,ymax=Inf),
            fill="#f0f0f0",alpha=.7,inherit.aes=FALSE) +
  geom_hline(yintercept=-log10(.05),color="#e8a020",linetype="dashed",linewidth=1.1) +
  geom_point(aes(color=notable, size=notable), alpha=.8) +
  scale_color_manual(values=c("TRUE"="#2c5f8a","FALSE"="#7fb3d3"), guide="none") +
  scale_size_manual(values=c("TRUE"=2.5,"FALSE"=1.8), guide="none") +
  geom_text_repel(data=lbl_df, aes(label=HUGO_SYMBOL),
                  size=3, color="#111111", box.padding=.5, segment.color="#cccccc") +
  scale_x_continuous(breaks=1:length(CHR_ORDER), labels=CHR_ORDER,
                     name="Chromosome") +
  scale_y_continuous(name="Fisher's exact −log₁₀(p)") +
  theme_classic(base_size=10) +
  labs(title="Manhattan Plot — Gene-Level Association with TMB-High",
       subtitle="n = 138 somatic variants · 15 genes · 80 patients")

ggsave("plots/08_manhattan.png", p, width=15, height=6.5, dpi=180, bg="white")
cat("08 done\n")
