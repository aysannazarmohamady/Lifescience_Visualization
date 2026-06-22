library(ggplot2); library(dplyr); library(tidyr); library(patchwork)

set.seed(6)
GENES_CNV <- c("PIK3CA","EGFR","MET","BRAF","CDKN2A","PTEN","KRAS","RB1","NF1","SMAD4")
CYTO <- c(PIK3CA="3q26.3",EGFR="7p11.2",MET="7q31.2",BRAF="7q34",CDKN2A="9p21.3",
           PTEN="10q23.3",KRAS="12p12.1",RB1="13q14.2",NF1="17q11.2",SMAD4="18q21.2")
TT_COLORS <- c(NSCLC="#2c5f8a",BRCA="#c0392b",HCC="#27ae60",CRC="#e67e22",PDAC="#8e44ad")
TT <- c(NSCLC=25,BRCA=13,HCC=12,CRC=16,PDAC=14)

pts <- data.frame(USUBJID=paste0("PT",sprintf("%03d",1:49)),
                  TUMORTYPE=sample(rep(names(TT), c(15,8,8,10,8))))

SCORES <- c(-4,-2,-1,0,1,3,4)
cnv <- expand.grid(USUBJID=pts$USUBJID, HUGO_SYMBOL=GENES_CNV) %>%
  filter(runif(n()) > 0.35) %>%
  mutate(SCORE=sample(SCORES, n(), TRUE, prob=c(.05,.1,.1,.35,.15,.1,.15)),
         GENE_LBL=paste0(HUGO_SYMBOL,"\n(",CYTO[HUGO_SYMBOL],")"),
         GENE_LBL=factor(GENE_LBL, levels=paste0(GENES_CNV,"\n(",CYTO[GENES_CNV],")"))) %>%
  left_join(pts, by="USUBJID") %>%
  mutate(USUBJID=factor(USUBJID, levels=pts$USUBJID))

tt_track <- pts %>% mutate(USUBJID=factor(USUBJID,levels=pts$USUBJID),
                            TUMORTYPE=factor(TUMORTYPE,levels=names(TT_COLORS)))

sum_cnv <- cnv %>% group_by(HUGO_SYMBOL,GENE_LBL) %>%
  summarise(n_amp=sum(SCORE>0), n_del=sum(SCORE<0), .groups="drop") %>%
  pivot_longer(c(n_amp,n_del), names_to="type", values_to="n") %>%
  mutate(n_signed=ifelse(type=="n_del",-n,n),
         GENE_LBL=factor(GENE_LBL,levels=paste0(GENES_CNV,"\n(",CYTO[GENES_CNV],")")))

p_tt <- ggplot(tt_track, aes(USUBJID,1,fill=TUMORTYPE)) +
  geom_tile(color="white") + scale_fill_manual(values=TT_COLORS) + theme_void()

p_heat <- ggplot(cnv, aes(USUBJID,GENE_LBL,fill=SCORE)) +
  geom_tile(color="white") +
  scale_fill_gradient2(low="#1a5276",mid="#f8f9fa",high="#922b21",midpoint=0,
                       limits=c(-4,4), name="CNV Score") +
  theme_classic(base_size=8) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.x=element_blank()) + labs(y=NULL)

p_bar <- ggplot(sum_cnv, aes(n_signed, GENE_LBL,
                              fill=ifelse(n_signed>0,"#c0392b","#1a5276"))) +
  geom_col(width=.6) + geom_vline(xintercept=0,color="#555555") +
  scale_fill_identity() +
  scale_x_continuous(name="n ← Del  Amp →", labels=abs) +
  theme_classic(base_size=8) +
  theme(axis.text.y=element_blank(), axis.ticks.y=element_blank(),
        axis.title.y=element_blank())

p <- (p_tt / p_heat + plot_layout(heights=c(1,12))) | p_bar +
  plot_layout(widths=c(13,2)) +
  plot_annotation(title="CNV Plot — Gene-Level Copy Number Alterations",
                  subtitle="n = 49 patients")

ggsave("plots/06_cnv.png", p, width=17, height=9, dpi=180, bg="white")
cat("06 done\n")
