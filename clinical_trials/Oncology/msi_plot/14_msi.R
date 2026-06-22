library(ggplot2); library(dplyr); library(tidyr); library(patchwork)

set.seed(14)
TT_COLORS <- c(NSCLC="#2c5f8a",BRCA="#c0392b",HCC="#27ae60",CRC="#e67e22",PDAC="#8e44ad")
THRESH <- 3.63

adsl <- data.frame(
  USUBJID   = paste0("PT",sprintf("%03d",1:80)),
  TUMORTYPE = rep(names(TT_COLORS), c(25,13,12,16,14)),
  MSI_SCORE = c(runif(76,.5,3.5), runif(4,3.7,4.4))
) %>% mutate(MSISTS=ifelse(MSI_SCORE>=THRESH,"MSI-H","MSS"),
             TUMORTYPE=factor(TUMORTYPE,levels=names(TT_COLORS))) %>%
  arrange(desc(MSI_SCORE)) %>% mutate(rank=row_number())

p14a <- ggplot(adsl, aes(rank, MSI_SCORE, fill=MSISTS)) +
  geom_col(width=.92, alpha=.88) +
  geom_hline(yintercept=THRESH, color="#c0392b", linetype="dashed", linewidth=1.3) +
  annotate("text", x=45, y=THRESH+.08,
           label=sprintf("MSI-H threshold = %.2f  |  MSI-H: %d/80",
                         THRESH, sum(adsl$MSISTS=="MSI-H")),
           size=3, color="#c0392b", hjust=0) +
  scale_fill_manual(values=c("MSI-H"="#c0392b","MSS"="#5b9bd5"), name="MSI Status") +
  theme_classic(base_size=10) + labs(x="Patients (ranked)", y="MSI Score",
                                      title="MSI Score — Ranked Across Cohort")

p14b <- ggplot(adsl, aes(TUMORTYPE, MSI_SCORE, fill=TUMORTYPE)) +
  geom_violin(alpha=.38, color=NA) +
  geom_boxplot(width=.12, outlier.shape=NA, color="#222222") +
  geom_jitter(aes(color=MSISTS), width=.10, size=2, alpha=.8) +
  geom_hline(yintercept=THRESH, color="#c0392b", linetype="dashed") +
  scale_fill_manual(values=TT_COLORS, guide="none") +
  scale_color_manual(values=c("MSI-H"="#c0392b","MSS"="#5b9bd5"), name="MSI Status") +
  theme_classic(base_size=10) + theme(axis.title.x=element_blank()) +
  labs(title="By Tumor Type")

msi_tt <- adsl %>% group_by(TUMORTYPE) %>%
  summarise(MSS=sum(MSISTS=="MSS")/n()*100,
            `MSI-H`=sum(MSISTS=="MSI-H")/n()*100, .groups="drop") %>%
  pivot_longer(c(MSS,`MSI-H`), names_to="status", values_to="pct") %>%
  mutate(TUMORTYPE=factor(TUMORTYPE,levels=rev(names(TT_COLORS))),
         status=factor(status,levels=c("MSS","MSI-H")))

p14c <- ggplot(msi_tt, aes(pct, TUMORTYPE, fill=status)) +
  geom_col(width=.65, alpha=.88) +
  geom_text(data=msi_tt%>%filter(pct>0),
            aes(label=paste0(round(pct),"%")),
            position=position_stack(vjust=.5), color="white", fontface="bold", size=3.2) +
  scale_fill_manual(values=c(MSS="#5b9bd5","MSI-H"="#c0392b"), name="MSI Status") +
  scale_x_continuous(name="% of Patients", limits=c(0,108)) +
  theme_classic(base_size=10) + theme(axis.title.y=element_blank()) +
  labs(title="MSI Status by Tumor Type")

p <- p14a / (p14b | p14c) + plot_layout(heights=c(1.2,1)) +
  plot_annotation(title="MSI Plot — Microsatellite Instability Across Cohort",
                  subtitle="n = 80 patients")

ggsave("plots/14_msi.png", p, width=16, height=10, dpi=180, bg="white")
cat("14 done\n")
