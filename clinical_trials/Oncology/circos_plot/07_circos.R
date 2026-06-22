library(ggplot2)

set.seed(7)
GENES8 <- c("ARID1A","PIK3CA","CDKN2A","PTEN","KRAS","SMAD4","TP53","KEAP1")
df <- data.frame(gene=factor(GENES8,levels=GENES8),
                 pct=c(12,11,11,12,30,9,48,10))
COLS <- c("#2980b9","#e67e22","#27ae60","#8e44ad",
          "#c0392b","#16a085","#6c3483","#1abc9c")

p <- ggplot(df, aes(gene, pct, fill=gene)) +
  geom_col(width=.7, alpha=.85) +
  geom_text(aes(y=pct+4, label=paste0(gene,"\n",pct,"%")), size=3, fontface="bold") +
  scale_fill_manual(values=setNames(COLS,GENES8)) +
  coord_polar() +
  scale_y_continuous(limits=c(0,65)) +
  theme_void() +
  theme(legend.position="none",
        plot.title=element_text(face="bold",size=13,hjust=.5),
        plot.subtitle=element_text(size=10,hjust=.5,color="#444444")) +
  labs(title="Circos Plot — Top 8 Mutated Genes",
       subtitle="n = 80 patients · Arc height = % altered")

ggsave("plots/07_circos.png", p, width=13, height=13, dpi=180, bg="white")
cat("07 done\n")
