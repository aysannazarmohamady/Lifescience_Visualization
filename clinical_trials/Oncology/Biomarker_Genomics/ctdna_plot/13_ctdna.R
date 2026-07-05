library(ggplot2); library(dplyr)

set.seed(13)
BOR_COLORS <- c(CR="#1a3f8f",PR="#4f8fd4",SD="#e8a020",PD="#c0392b")
DAYS <- c(0,21,42,63,84,126,168)

ctdna <- do.call(rbind, lapply(paste0("PT",sprintf("%03d",1:62)), function(id) {
  bor <- sample(c("CR","PR","SD","PD"),1,prob=c(.1,.35,.35,.2))
  base_trend <- c(CR=-80,PR=-50,SD=5,PD=60)[bor]
  do.call(rbind, lapply(DAYS, function(d) {
    data.frame(USUBJID=id, BESTRSPC=bor, ADTN=d,
               PCHG=base_trend*(d/168) + rnorm(1,0,20))
  }))
})) %>% mutate(BESTRSPC=factor(BESTRSPC, levels=c("CR","PR","SD","PD")))

ctdna_sum <- ctdna %>% group_by(BESTRSPC, ADTN) %>%
  summarise(mean=mean(PCHG), se=sd(PCHG)/sqrt(n()), n=n(), .groups="drop") %>%
  mutate(lo=mean-qt(.975,n-1)*se, hi=mean+qt(.975,n-1)*se)

bor_n <- ctdna %>% group_by(BESTRSPC) %>% summarise(n=n_distinct(USUBJID))
bor_lbl <- setNames(paste0(bor_n$BESTRSPC,"\n(n=",bor_n$n,")"), bor_n$BESTRSPC)

p <- ggplot() +
  geom_hline(yintercept=0, color="#888888", linetype="dashed") +
  geom_line(data=ctdna, aes(ADTN,PCHG,group=USUBJID,color=BESTRSPC), alpha=.18) +
  geom_ribbon(data=ctdna_sum, aes(ADTN,ymin=lo,ymax=hi,fill=BESTRSPC), alpha=.45) +
  geom_line(data=ctdna_sum, aes(ADTN,mean,color=BESTRSPC), linewidth=2.2) +
  geom_point(data=ctdna_sum, aes(ADTN,mean,color=BESTRSPC), size=2.8) +
  facet_wrap(~BESTRSPC, nrow=1, labeller=labeller(BESTRSPC=bor_lbl)) +
  scale_color_manual(values=BOR_COLORS, guide="none") +
  scale_fill_manual(values=BOR_COLORS, guide="none") +
  scale_x_continuous(breaks=DAYS, name="Day") +
  scale_y_continuous(name="ctDNA VAF % Change from Baseline", limits=c(-155,185)) +
  theme_classic(base_size=10) + theme(strip.text=element_text(face="bold",size=12)) +
  labs(title="ctDNA Dynamics — Longitudinal Change by Best Overall Response",
       subtitle="Treatment Arm · n = 62 patients · Mean ± 95% CI")

ggsave("plots/13_ctdna.png", p, width=17, height=7, dpi=180, bg="white")
cat("13 done\n")
