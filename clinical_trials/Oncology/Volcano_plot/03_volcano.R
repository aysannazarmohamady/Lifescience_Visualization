library(ggplot2); library(ggrepel); library(dplyr)

set.seed(3)
PARAMS <- c("TMB","CTDNA","PDL1","MSI","CEA","CA125","CD8","CD4","NK","TREG","IFNg")
VISITS <- c("C1","C2","C3","C4","C6","C8")
FAM <- c(TMB="Tumor genomics",CTDNA="Tumor genomics",PDL1="Tumor genomics",MSI="Tumor genomics",
         CEA="Serum tumor marker",CA125="Serum tumor marker",CD8="Immune cell subset",
         CD4="Immune cell subset",NK="Immune cell subset",TREG="Immune cell subset",IFNg="Cytokine")
FAM_COLORS <- c("Tumor genomics"="#2c5f8a","Serum tumor marker"="#c0392b",
                "Immune cell subset"="#27ae60","Cytokine"="#8e44ad")

vol <- expand.grid(PARAMCD=PARAMS, AVISIT=VISITS, stringsAsFactors=FALSE) %>%
  mutate(mean_pchg = rnorm(n(), mean=c(TMB=30,CTDNA=-20,PDL1=50,MSI=10,CEA=80,
                                        CA125=60,CD8=-15,CD4=-5,NK=-10,TREG=20,IFNg=25)[PARAMCD], sd=15),
         pval = runif(n(), 0.001, 0.5),
         pval = ifelse(PARAMCD %in% c("CEA","PDL1","CTDNA"), pval*0.05, pval),
         log10p = -log10(pmax(pval,1e-10)),
         family = FAM[PARAMCD],
         label = paste0(PARAMCD," (",AVISIT,")"),
         sig = pval < 0.05)

p <- ggplot(vol, aes(mean_pchg, log10p, color=family)) +
  annotate("rect", xmin=-Inf, xmax=Inf, ymin=-log10(0.05), ymax=Inf,
           fill="#fff3e0", alpha=0.5) +
  geom_hline(yintercept=-log10(0.05), linetype="dashed", color="#aaaaaa") +
  geom_vline(xintercept=0, color="#cccccc") +
  geom_point(aes(size=sig), alpha=0.85) +
  scale_size_manual(values=c("TRUE"=3.5,"FALSE"=2.2), guide="none") +
  geom_text_repel(data=vol %>% filter(log10p>1),
                  aes(label=label), size=2.8, color="#1a1a1a",
                  box.padding=0.5, segment.color="#bbbbbb", max.overlaps=25) +
  scale_color_manual(values=FAM_COLORS, name="Biomarker Family") +
  labs(x="Mean % Change from Baseline", y="−log₁₀(p-value)",
       title="Volcano Plot — On-Treatment Biomarker Change",
       subtitle="Treatment Arm · 11 biomarkers × visit") +
  theme_classic(base_size=10) + theme(legend.position=c(0.85,0.85))

ggsave("plots/03_volcano.png", p, width=13, height=9, dpi=180, bg="white")
cat("03 done\n")
