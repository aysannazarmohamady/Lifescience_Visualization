library(ggplot2); library(dplyr); library(tidyr); library(patchwork)

set.seed(11)
SIGS <- c("SBS3","SBS13","SBS40","SBS4","SBS1","SBS44",
          "SBS6","SBS5","SBS18","SBS22","SBS2","SBS17")
SIG_COLORS <- c(SBS3="#3498db",SBS13="#aed6f1",SBS40="#f39c12",SBS4="#e74c3c",
                SBS1="#f1948a",SBS44="#8e44ad",SBS6="#d7bde2",SBS5="#1abc9c",
                SBS18="#d5dbdb",SBS22="#f9e79f",SBS2="#85c1e9",SBS17="#922b21")
TT_COLORS <- c(NSCLC="#2c5f8a",BRCA="#c0392b",HCC="#27ae60",CRC="#e67e22",PDAC="#8e44ad")
TT <- rep(names(TT_COLORS), c(25,13,12,16,14))

make_weights <- function(n, k) {
  w <- abs(rnorm(k)); w/sum(w)
}

adsig <- do.call(rbind, lapply(1:80, function(i) {
  w <- make_weights(1, length(SIGS))
  data.frame(USUBJID=paste0("PT",sprintf("%03d",i)),
             TUMORTYPE=TT[i], SIG_NAME=SIGS, SIG_WEIGHT=w)
})) %>% mutate(SIG_NAME=factor(SIG_NAME,levels=SIGS),
               TUMORTYPE=factor(TUMORTYPE,levels=names(TT_COLORS)))

dom <- adsig %>% group_by(USUBJID) %>% slice_max(SIG_WEIGHT,n=1) %>%
  arrange(TUMORTYPE, desc(SIG_WEIGHT))
pt_ord <- dom$USUBJID
adsig  <- adsig %>% mutate(USUBJID=factor(USUBJID,levels=pt_ord))

p11a <- ggplot(adsig, aes(USUBJID, SIG_WEIGHT, fill=SIG_NAME)) +
  geom_col(width=1, position="stack") +
  scale_fill_manual(values=SIG_COLORS, name="SBS Signature") +
  scale_y_continuous(name="Signature Fraction", limits=c(0,1.02)) +
  theme_classic(base_size=8) +
  theme(axis.text.x=element_blank(), axis.ticks.x=element_blank(),
        axis.title.x=element_blank(), legend.position="bottom") +
  guides(fill=guide_legend(ncol=4)) +
  labs(title="Per-Patient SBS Decomposition")

tt_mean <- adsig %>% group_by(SIG_NAME,TUMORTYPE) %>%
  summarise(mean_w=mean(SIG_WEIGHT), .groups="drop")

p11b <- ggplot(tt_mean, aes(mean_w, factor(SIG_NAME,levels=rev(SIGS)), fill=TUMORTYPE)) +
  geom_col(position="stack", width=.7, alpha=.85) +
  scale_fill_manual(values=TT_COLORS, name="Tumor Type") +
  scale_x_continuous(name="Mean Weight") +
  theme_classic(base_size=8) +
  theme(axis.title.y=element_blank()) +
  labs(title="By Tumor Type")

p <- p11a + p11b + plot_layout(widths=c(13,4)) +
  plot_annotation(title="Mutational Signature Plot — COSMIC SBS Contributions",
                  subtitle="n = 80 patients · 12 signatures")

ggsave("plots/11_signature.png", p, width=18, height=10, dpi=180, bg="white")
cat("11 done\n")
