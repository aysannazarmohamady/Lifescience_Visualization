library(rlang)
library(dplyr)
library(tidyr)
library(tibble)
library(ggplot2)
library(patchwork)

adtr <- read.csv("ADTR.csv")
adsl <- read.csv("ADSL.csv")

trt <- adsl |> filter(ARM == "TREATMENT") |> select(USUBJID, BESTRSPC)

df <- adtr |>
  filter(PARAMCD == "SUMDIAM", ANL01FL == "Y",
         USUBJID %in% trt$USUBJID) |>
  left_join(trt, by = "USUBJID") |>
  mutate(
    VL = factor(
      dplyr::recode(as.character(AVISITN),
        "0"="BL","1"="C1","2"="C2","3"="C3","4"="C4","5"="C5","6"="C6",
        "7"="C7","8"="C8","9"="C9","10"="C10","11"="C11","12"="C12"),
      levels = c("BL","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12")),
    BESTRSPC  = factor(BESTRSPC,  levels = c("CR","PR","SD","PD")),
    TUMORTYPE = factor(TUMORTYPE, levels = c("NSCLC","BRCA","CRC","HCC","PDAC")))

BOR_COLOR <- c(CR="#2e8b57", PR="#4169e1", SD="#d4a017", PD="#b22222")
BOR_LONG  <- c(CR="Complete Response (CR)", PR="Partial Response (PR)",
               SD="Stable Disease (SD)",    PD="Progressive Disease (PD)")
TT_COLOR  <- c(NSCLC="#4169e1", BRCA="#b22222", CRC="#9467bd", HCC="#d4a017", PDAC="#2e8b57")

t_ci <- function(x) qt(0.975, df = length(x) - 1) * sd(x, na.rm=TRUE) / sqrt(sum(!is.na(x)))

base_thm <- function() {
  theme_classic(base_size = 11) +
    theme(axis.title       = element_text(size = 11.5),
          axis.text        = element_text(size = 10),
          legend.text      = element_text(size = 10),
          legend.title     = element_blank(),
          legend.key.width = unit(1.6, "lines"),
          plot.title       = element_text(size = 13.5, face = "bold", hjust = 0),
          panel.grid.major.y = element_line(linewidth = 0.3, colour = "grey85"),
          panel.grid.major.x = element_blank())
}

vis_labs <- c("BL","C1","C2","C3","C4","C5","C6","C7","C8","C9","C10","C11","C12")

nar_plot <- function(data, grp_var, grp_levels, grp_colors,
                     vis = 0:12, x_labs = vis_labs) {
  counts <- data |>
    filter(AVISITN %in% vis) |>
    group_by(.data[[grp_var]], AVISITN) |>
    summarise(n = n_distinct(USUBJID), .groups = "drop") |>
    tidyr::complete(!!sym(grp_var) := grp_levels, AVISITN = vis,
                    fill = list(n = 0)) |>
    mutate(x_pos = match(AVISITN, vis) - 1)

  n_grp <- length(grp_levels)

  ggplot(counts,
         aes(x = x_pos, y = .data[[grp_var]], label = n,
             colour = .data[[grp_var]])) +
    geom_text(size = 3.2, fontface = "bold") +
    annotate("text", x = -0.9, y = n_grp + 0.6,
             label = "Number at risk", hjust = 0,
             size = 3.3, fontface = "bold.italic", colour = "grey30") +
    scale_colour_manual(values = grp_colors, breaks = grp_levels) +
    scale_y_discrete(limits = rev(grp_levels)) +
    scale_x_continuous(breaks = seq_along(vis) - 1,
                       labels = x_labs[seq_along(vis)],
                       expand = expansion(add = c(1.2, 0.5))) +
    coord_cartesian(clip = "off") +
    theme_void() +
    theme(axis.text.y     = element_text(size = 9.5, face = "bold",
                                         colour = grp_colors[rev(grp_levels)],
                                         hjust = 1, margin = margin(r = 4)),
          axis.text.x     = element_text(size = 8.5, colour = "grey40",
                                         face = "bold"),
          plot.margin     = margin(2, 10, 2, 10),
          legend.position = "none")
}

bor_n <- df |> distinct(USUBJID, BESTRSPC) |> count(BESTRSPC) |> deframe()
bor_leg <- setNames(paste0(BOR_LONG, " (n=", bor_n[names(BOR_LONG)], ")"), names(BOR_LONG))

agg_A <- df |> filter(AVISITN %in% 0:12) |>
  group_by(BESTRSPC, AVISITN) |>
  summarise(m=mean(AVAL,na.rm=T), ci=t_ci(AVAL), .groups="drop") |>
  mutate(lo=m-ci, hi=m+ci)

bl_mean <- df |> filter(AVISITN==0) |> pull(AVAL) |> mean(na.rm=T)

pA <- ggplot(agg_A, aes(AVISITN, colour=BESTRSPC, fill=BESTRSPC)) +
  geom_ribbon(aes(ymin=lo, ymax=hi), alpha=0.18, colour=NA) +
  geom_line(aes(y=m), linewidth=1.1) + geom_point(aes(y=m), size=2.5) +
  geom_hline(yintercept=bl_mean, linetype="dashed", colour="grey50", linewidth=0.6) +
  annotate("text", x=12.7, y=bl_mean,
           label=sprintf("BL mean = %.1f mm", bl_mean),
           size=3.2, colour="grey45", hjust=0) +
  scale_x_continuous(breaks=0:12, labels=vis_labs,
                     expand=expansion(add=c(0.3,1.8))) +
  scale_y_continuous(limits=c(0,100)) +
  scale_colour_manual(values=BOR_COLOR, labels=bor_leg) +
  scale_fill_manual(values=BOR_COLOR,   labels=bor_leg) +
  labs(title="A   Mean SLD by Best Overall Response  (95% CI ribbon)",
       y="Mean SLD (mm)", x=NULL) +
  base_thm() +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(ncol=2), fill=guide_legend(ncol=2))

narA <- nar_plot(df, "BESTRSPC", c("CR","PR","SD","PD"), BOR_COLOR)

agg_B <- df |> filter(AVISITN %in% 1:12) |>
  group_by(BESTRSPC, AVISITN) |>
  summarise(m=mean(PCHG,na.rm=T),
            se=sd(PCHG,na.rm=T)/sqrt(sum(!is.na(PCHG))),
            .groups="drop") |>
  mutate(lo=m-se, hi=m+se)

bor_short <- setNames(paste0(names(bor_n)," (n=",bor_n,")"), names(bor_n))

pB <- ggplot(agg_B, aes(AVISITN, colour=BESTRSPC, fill=BESTRSPC)) +
  geom_ribbon(aes(ymin=lo,ymax=hi), alpha=0.18, colour=NA) +
  geom_line(aes(y=m), linewidth=1.1) + geom_point(aes(y=m), size=2.5) +
  geom_hline(yintercept=c(0), colour="grey60", linewidth=0.5) +
  geom_hline(yintercept=20,  colour="#e74c3c", linetype="dashed", linewidth=0.8) +
  geom_hline(yintercept=-30, colour="#27ae60", linetype="dashed", linewidth=0.8) +
  annotate("text",x=12.7,y=20, label="+20%\nPD threshold",
           size=2.9,colour="#e74c3c",hjust=0) +
  annotate("text",x=12.7,y=-30,label="-30%\nPR threshold",
           size=2.9,colour="#27ae60",hjust=0) +
  scale_x_continuous(breaks=1:12, labels=paste0("C",1:12),
                     expand=expansion(add=c(0.3,2.0))) +
  scale_y_continuous(limits=c(-85,40)) +
  scale_colour_manual(values=BOR_COLOR, labels=bor_short) +
  scale_fill_manual(values=BOR_COLOR,   labels=bor_short) +
  labs(title="B   Mean % Change in SLD by BOR  (\u00b11 SE ribbon)",
       y="Mean % Change from Baseline", x=NULL) +
  base_thm() +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(ncol=2), fill=guide_legend(ncol=2))

narB <- nar_plot(df, "BESTRSPC", c("CR","PR","SD","PD"), BOR_COLOR,
                 vis=1:12, x_labs=vis_labs[-1])

bor_leg_dsh <- setNames(paste0(BOR_LONG," (n=",bor_n[names(BOR_LONG)],")"), names(BOR_LONG))

agg_C <- df |> filter(AVISITN %in% 0:12) |>
  group_by(BESTRSPC, AVISITN) |>
  summarise(m=median(AVAL,na.rm=T),
            lo=quantile(AVAL,.25,na.rm=T),
            hi=quantile(AVAL,.75,na.rm=T), .groups="drop")

bl_med <- df |> filter(AVISITN==0) |> pull(AVAL) |> median(na.rm=T)

pC <- ggplot(agg_C, aes(AVISITN, colour=BESTRSPC, fill=BESTRSPC)) +
  geom_ribbon(aes(ymin=lo,ymax=hi), alpha=0.18, colour=NA) +
  geom_line(aes(y=m), linewidth=1.1, linetype="dashed") +
  geom_point(aes(y=m), size=2.5) +
  geom_hline(yintercept=bl_med, linetype="dashed", colour="grey50", linewidth=0.6) +
  annotate("text", x=12.7, y=bl_med,
           label=sprintf("BL median = %.1f mm", bl_med),
           size=3.2, colour="grey45", hjust=0) +
  scale_x_continuous(breaks=0:12, labels=vis_labs,
                     expand=expansion(add=c(0.3,1.8))) +
  scale_y_continuous(limits=c(0,75)) +
  scale_colour_manual(values=BOR_COLOR, labels=bor_leg_dsh) +
  scale_fill_manual(values=BOR_COLOR,   labels=bor_leg_dsh) +
  labs(title="C   Median SLD by Best Overall Response  (IQR ribbon)",
       y="Median SLD (mm)", x=NULL) +
  base_thm() +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(ncol=2), fill=guide_legend(ncol=2))

narC <- nar_plot(df, "BESTRSPC", c("CR","PR","SD","PD"), BOR_COLOR)

tt_n <- df |> distinct(USUBJID, TUMORTYPE) |> count(TUMORTYPE) |> deframe()
tt_leg <- setNames(paste0(names(tt_n)," (n=",tt_n,")"), names(tt_n))

agg_D <- df |> filter(AVISITN %in% 0:12) |>
  group_by(TUMORTYPE, AVISITN) |>
  summarise(m=mean(AVAL,na.rm=T), ci=t_ci(AVAL), .groups="drop") |>
  mutate(lo=m-ci, hi=m+ci)

pD <- ggplot(agg_D, aes(AVISITN, colour=TUMORTYPE, fill=TUMORTYPE)) +
  geom_ribbon(aes(ymin=lo,ymax=hi), alpha=0.15, colour=NA) +
  geom_line(aes(y=m), linewidth=1.1) + geom_point(aes(y=m), size=2.5) +
  scale_x_continuous(breaks=0:12, labels=vis_labs,
                     expand=expansion(add=c(0.3,0.5))) +
  scale_y_continuous(limits=c(0,80)) +
  scale_colour_manual(values=TT_COLOR, labels=tt_leg) +
  scale_fill_manual(values=TT_COLOR,   labels=tt_leg) +
  labs(title="D   Mean SLD by Tumor Type  (95% CI ribbon)",
       y="Mean SLD (mm)", x=NULL) +
  base_thm() +
  theme(legend.position="bottom") +
  guides(colour=guide_legend(ncol=3), fill=guide_legend(ncol=3))

narD <- nar_plot(df, "TUMORTYPE", c("NSCLC","BRCA","CRC","HCC","PDAC"), TT_COLOR)

wrap <- function(p, nar, nar_h=0.22)
  p / nar + plot_layout(heights=c(1, nar_h))

final <- (wrap(pA, narA) | wrap(pB, narB)) /
         (wrap(pC, narC) | wrap(pD, narD, nar_h=0.30)) +
  plot_annotation(
    title    = "Tumor Burden Plot  \u00b7  ONCVIZ-001  \u00b7  Vizatinib 300 mg QD  \u00b7  Treatment Arm (N=62)",
    subtitle = "RECIST 1.1  \u00b7  Sum of Longest Diameters (SLD, mm)  \u00b7  Data cutoff: 05 Mar 2026",
    caption  = "SLD = Sum of Longest Diameters  \u00b7  BOR = Best Overall Response  \u00b7  CI = Confidence Interval  \u00b7  SE = Standard Error  \u00b7  IQR = Interquartile Range  \u00b7  RECIST 1.1  \u00b7  ONCVIZ-001 ADaM v1",
    theme = theme(
      plot.title    = element_text(size=15.5, face="bold", hjust=0.5),
      plot.subtitle = element_text(size=9.5,  hjust=0.5, colour="grey40"),
      plot.caption  = element_text(size=8,    hjust=0.5, colour="grey50")))

ggsave("/mnt/user-data/outputs/tumor_burden_plot_R.png",
       final, width=22, height=22, dpi=180, bg="white")
message("Saved.")
