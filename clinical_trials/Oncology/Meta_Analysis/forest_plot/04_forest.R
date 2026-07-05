library(ggplot2); library(dplyr); library(survival)

set.seed(4)
n <- 80
adsl <- data.frame(
  USUBJID   = paste0("PT",sprintf("%03d",1:n)),
  ARM       = c(rep("TREATMENT",62), rep("CONTROL",18)),
  TUMORTYPE = sample(c("NSCLC","CRC","BRCA","HCC","PDAC"), n, TRUE),
  AGEGR1    = sample(c("<65",">=65"), n, TRUE),
  SEX       = sample(c("M","F"), n, TRUE),
  TMBHIGH   = sample(c("Y","N"), n, TRUE, prob=c(.4,.6)),
  KRASMUT   = sample(c("Y","N"), n, TRUE, prob=c(.35,.65)),
  TP53MUT   = sample(c("Y","N"), n, TRUE, prob=c(.55,.45))
)
adtte <- adsl %>% mutate(
  AVAL  = rexp(n, rate=ifelse(ARM=="TREATMENT",1/18,1/12)),
  CNSR  = sample(c(0,1),n,TRUE,prob=c(.7,.3)),
  EVENT = 1-CNSR,
  TRT   = as.integer(ARM=="TREATMENT")
)

cox_hr <- function(df) {
  nt <- sum(df$ARM=="TREATMENT"); nc <- sum(df$ARM=="CONTROL")
  if(nt<2 | nc<2) return(list(hr=NA,lo=NA,hi=NA,nt=nt,nc=nc))
  fit <- tryCatch(coxph(Surv(AVAL,EVENT)~TRT,data=df),error=function(e)NULL)
  if(is.null(fit)) return(list(hr=NA,lo=NA,hi=NA,nt=nt,nc=nc))
  ci <- exp(confint(fit))
  list(hr=exp(coef(fit)), lo=ci[1], hi=ci[2], nt=nt, nc=nc)
}

rows <- list(
  list(l="Overall",    d=adtte,                                ov=TRUE),
  list(l="NSCLC",      d=adtte%>%filter(TUMORTYPE=="NSCLC"),   ov=FALSE),
  list(l="CRC",        d=adtte%>%filter(TUMORTYPE=="CRC"),     ov=FALSE),
  list(l="BRCA",       d=adtte%>%filter(TUMORTYPE=="BRCA"),    ov=FALSE),
  list(l="< 65 yrs",   d=adtte%>%filter(AGEGR1=="<65"),        ov=FALSE),
  list(l=">= 65 yrs",  d=adtte%>%filter(AGEGR1==">=65"),       ov=FALSE),
  list(l="Male",       d=adtte%>%filter(SEX=="M"),             ov=FALSE),
  list(l="Female",     d=adtte%>%filter(SEX=="F"),             ov=FALSE),
  list(l="TMB-High",   d=adtte%>%filter(TMBHIGH=="Y"),         ov=FALSE),
  list(l="TMB-Low",    d=adtte%>%filter(TMBHIGH=="N"),         ov=FALSE),
  list(l="KRAS Mut",   d=adtte%>%filter(KRASMUT=="Y"),         ov=FALSE),
  list(l="KRAS WT",    d=adtte%>%filter(KRASMUT=="N"),         ov=FALSE)
)

fdf <- do.call(rbind, lapply(seq_along(rows), function(i) {
  r <- cox_hr(rows[[i]]$d)
  data.frame(y=i, label=rows[[i]]$l, hr=r$hr, lo=r$lo, hi=r$hi,
             overall=rows[[i]]$ov,
             sig=!is.na(r$hr) & !is.na(r$hi) & r$hi<1)
}))

p <- ggplot(fdf, aes(y=reorder(label,y))) +
  geom_vline(xintercept=1, linetype="dashed", color="#888888") +
  geom_errorbarh(data=fdf%>%filter(!is.na(hr)),
                 aes(xmin=pmax(lo,.05,na.rm=TRUE), xmax=pmin(hi,10,na.rm=TRUE), height=0),
                 linewidth=0.9, color="#333333") +
  geom_point(data=fdf%>%filter(!overall,!sig,!is.na(hr)), aes(x=hr),
             shape=16, size=2.5, color="#555555") +
  geom_point(data=fdf%>%filter(!overall,sig,!is.na(hr)),  aes(x=hr),
             shape=16, size=2.8, color="#c0392b") +
  geom_point(data=fdf%>%filter(overall,!is.na(hr)),       aes(x=hr),
             shape=18, size=5,   color="#1a3f8f") +
  scale_x_log10(breaks=c(.1,.25,.5,1,2,4),
                name="HR (95% CI)\n← Favors Treatment   Favors Control →") +
  theme_classic(base_size=10) +
  theme(axis.title.y=element_blank()) +
  labs(title="Forest Plot — Overall Survival Subgroup Analysis",
       subtitle="Treatment (n=62) vs Control (n=18)",
       caption="Diamond = overall; red = CI excludes 1")

ggsave("plots/04_forest.png", p, width=13, height=9, dpi=180, bg="white")
cat("04 done\n")
