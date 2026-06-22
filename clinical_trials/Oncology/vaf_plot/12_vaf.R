library(ggplot2); library(dplyr); library(patchwork)

set.seed(12)
TT_COLORS <- c(NSCLC="#2c5f8a",BRCA="#c0392b",HCC="#27ae60",CRC="#e67e22",PDAC="#8e44ad")
CLONAL_C  <- c(Y="#1a3f8f",N="#aec6ef")
GENES8    <- c("PTEN","KEAP1","PIK3CA","SMAD4","CDKN2A","KRAS","ARID1A","TP53")

TT <- rep(names(TT_COLORS), c(25,13,12,16,14))
point <- do.call(rbind, lapply(1:138, function(i) {
  g <- sample(c(GENES8,"OTHER"), 1, prob=c(rep(.1,8),.2))
  tt <- sample(TT, 1)
  data.frame(HUGO_SYMBOL=g, TUMORTYPE=tt,
             VAF=round(runif(1,.05,.80),3),
             CLONAL=sample(c("Y","N"),1,prob=c(.65,.35)))
})) %>% mutate(TUMORTYPE=factor(TUMORTYPE,levels=names(TT_COLORS)))

p12a <- ggplot(point, aes(VAF, fill=CLONAL)) +
  geom_histogram(binwidth=.05, alpha=.85, position="stack", color="white") +
  geom_vline(xintercept=.4, color="#c0392b", linetype="dashed", linewidth=1.4) +
  scale_fill_manual(values=CLONAL_C, labels=c(Y="Clonal",N="Subclonal"), name="Clonality") +
  theme_classic(base_size=10) + labs(x="VAF", y="Count", title="VAF Distribution")

p12b <- ggplot(point, aes(TUMORTYPE, VAF, fill=TUMORTYPE)) +
  geom_violin(alpha=.38, color=NA) +
  geom_boxplot(width=.12, outlier.shape=NA, color="#222222") +
  geom_jitter(aes(color=CLONAL), width=.08, size=1.8, alpha=.72) +
  geom_hline(yintercept=.4, color="#c0392b", linetype="dashed") +
  scale_fill_manual(values=TT_COLORS, guide="none") +
  scale_color_manual(values=CLONAL_C, labels=c(Y="Clonal",N="Subclonal"), name="Clonality") +
  theme_classic(base_size=10) + theme(axis.title.x=element_blank()) +
  labs(title="VAF by Tumor Type")

p12c <- ggplot(point %>% filter(HUGO_SYMBOL %in% GENES8) %>%
                 mutate(HUGO_SYMBOL=factor(HUGO_SYMBOL,levels=GENES8)),
               aes(HUGO_SYMBOL, VAF)) +
  geom_violin(fill="#7f8fa6", alpha=.38, color=NA) +
  geom_boxplot(width=.12, outlier.shape=NA, color="#222222") +
  geom_jitter(aes(color=CLONAL), width=.08, size=1.8, alpha=.72) +
  geom_hline(yintercept=.4, color="#c0392b", linetype="dashed") +
  scale_color_manual(values=CLONAL_C, name="Clonality") +
  theme_classic(base_size=10) +
  theme(axis.text.x=element_text(angle=22,hjust=1), axis.title.x=element_blank()) +
  labs(title="VAF by Gene (Top 8)")

p <- p12a / (p12b | p12c) + plot_layout(heights=c(1,1.2)) +
  plot_annotation(title="VAF Distribution — Clonal Architecture",
                  subtitle="n = 138 mutations · 66 patients")

ggsave("plots/12_vaf.png", p, width=17, height=12, dpi=180, bg="white")
cat("12 done\n")
