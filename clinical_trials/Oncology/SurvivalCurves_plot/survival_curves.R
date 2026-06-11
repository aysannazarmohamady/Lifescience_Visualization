suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(grid)
  library(gridExtra)
  library(scales)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR   <- "./Data/V1"
OUTPUT_DIR <- "./Outputs"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Load data ─────────────────────────────────────────────────────────────────
adtte <- read.csv(file.path(DATA_DIR, "ADTTE.csv"), stringsAsFactors = FALSE)
adsl  <- read.csv(file.path(DATA_DIR, "ADSL.csv"),  stringsAsFactors = FALSE)

adtte$EVENT <- 1L - as.integer(adtte$CNSR)

# ── Global palette (matches Python outputs) ───────────────────────────────────
C_TRT    <- "#1A73E8"   # Treatment  – blue
C_CTRL   <- "#D44C36"   # Control    – red
C_DFS    <- "#20776A"   # DFS single arm – teal
C_GRID   <- "#E8E8E8"
C_TEXT   <- "#161616"
C_MUT    <- "#6F6F6F"
C_PURPLE <- "#8E44AD"

# ── Shared ggplot2 theme ──────────────────────────────────────────────────────
theme_oncoviz <- function(base_size = 11) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      # panel
      panel.grid.major   = element_line(color = C_GRID, linewidth = 0.7),
      panel.grid.minor.y = element_line(color = "#F5F5F5", linewidth = 0.35),
      panel.grid.major.x = element_line(color = C_GRID,   linewidth = 0.7),
      panel.grid.minor.x = element_blank(),
      panel.background   = element_rect(fill = "white", color = NA),
      plot.background    = element_rect(fill = "white", color = NA),
      # axes
      axis.line          = element_line(color = "#CCCCCC"),
      axis.ticks         = element_blank(),
      axis.text          = element_text(color = C_TEXT, size = rel(0.95)),
      axis.title         = element_text(color = C_TEXT, face = "bold",
                                        size = rel(1.05)),
      axis.title.y       = element_text(margin = margin(r = 10)),
      # title / subtitle
      plot.title         = element_text(face = "bold", size = rel(1.25),
                                        color = C_TEXT,
                                        margin = margin(b = 3)),
      plot.subtitle      = element_text(size = rel(0.82), color = C_MUT,
                                        margin = margin(b = 4)),
      plot.caption       = element_text(size = rel(0.65), color = "#AAAAAA",
                                        hjust = 0),
      plot.margin        = margin(14, 14, 6, 14),
      # legend
      legend.background  = element_rect(fill = "white", color = "#DDDDDD",
                                        linewidth = 0.8),
      legend.key         = element_rect(fill = "white"),
      legend.text        = element_text(size = rel(0.92)),
      legend.title       = element_text(face = "bold", size = rel(0.92)),
      legend.margin      = margin(6, 8, 6, 8)
    )
}

# ── Helper: compute HR + 95 % CI from coxph ──────────────────────────────────
get_hr <- function(df_param, time_col = "AVALM", event_col = "EVENT",
                   arm_col = "ARM", ref = "CONTROL") {
  df <- df_param
  df[[arm_col]] <- relevel(factor(df[[arm_col]]), ref = ref)
  fit  <- coxph(
    as.formula(paste0("Surv(", time_col, ",", event_col, ") ~ ", arm_col)),
    data = df
  )
  smry <- summary(fit)
  hr   <- smry$conf.int[1, "exp(coef)"]
  lo   <- smry$conf.int[1, "lower .95"]
  hi   <- smry$conf.int[1, "upper .95"]
  pv   <- smry$sctest["pvalue"]
  list(hr = hr, lo = lo, hi = hi, pval = pv)
}

# ── Helper: landmark survival rate ───────────────────────────────────────────
lm_rate <- function(km_fit, t) {
  sv <- summary(km_fit, times = t, extend = TRUE)
  round(sv$surv * 100, 1)
}

# ── Helper: format p-value ────────────────────────────────────────────────────
fmt_p <- function(p) {
  if (p < 0.001) "P < 0.001" else sprintf("P = %.3f", p)
}

# ── Helper: stats annotation grob (bottom-right inset) ───────────────────────
stats_grob <- function(lines, x = 0.97, y = 0.03,
                       hjust = 1, vjust = 0,
                       fontsize = 9, col_val = C_TRT) {
  txt <- paste(lines, collapse = "\n")
  annotation_custom(
    grob = textGrob(
      txt,
      x = unit(x, "npc"), y = unit(y, "npc"),
      hjust = hjust, vjust = vjust,
      gp = gpar(
        fontsize  = fontsize,
        fontfamily = "mono",
        col        = C_TEXT,
        lineheight = 1.45
      )
    )
  )
}

# ══════════════════════════════════════════════════════════════════════════════
#  1.  KAPLAN-MEIER CURVE  (generic KM, Treatment vs Control, OS endpoint)
# ══════════════════════════════════════════════════════════════════════════════
cat("Building KM Curve...\n")

km_df <- adtte |>
  filter(PARAMCD == "OS") |>
  mutate(ARM = trimws(ARM))

km_fit <- survfit(
  Surv(AVALM, EVENT) ~ ARM,
  data    = km_df,
  conf.type = "log-log"
)

hr_km  <- get_hr(km_df)
pstr   <- fmt_p(hr_km$pval)

# survminer base plot
km_plot <- ggsurvplot(
  km_fit,
  data           = km_df,
  palette        = c(C_CTRL, C_TRT),           # Control first (alphabetical)
  size           = 1.1,
  conf.int       = TRUE,
  conf.int.alpha = 0.10,
  censor.shape   = "|",
  censor.size    = 5,
  risk.table     = TRUE,
  risk.table.height = 0.22,
  risk.table.col    = "strata",
  risk.table.fontsize = 3.6,
  tables.theme   = theme_cleantable(),
  break.x.by    = 6,
  xlim          = c(0, 44),
  ylim          = c(0, 1.0),
  xlab          = "Time (Months)",
  ylab          = "Kaplan–Meier Survival Probability",
  legend.labs   = c(
    sprintf("Control    (n=18, events=13)"),
    sprintf("Treatment  (n=62, events=31)")
  ),
  legend.title  = "",
  ggtheme       = theme_oncoviz()
)

# Add title / stats annotation on main plot
n_trt  <- sum(km_df$ARM == "TREATMENT")
n_ctrl <- sum(km_df$ARM == "CONTROL")
med_t  <- summary(km_fit)$table["ARM=TREATMENT", "median"]
med_c  <- summary(km_fit)$table["ARM=CONTROL",   "median"]

km_plot$plot <- km_plot$plot +
  labs(
    title    = "Kaplan–Meier Estimate of Overall Survival",
    subtitle = paste0(
      "Phase I/II  \u00b7  Vizatinib vs Control  \u00b7  ITT Population  ",
      "\u00b7  Data cut-off: 05 Mar 2026"
    ),
    caption  = paste0(
      "ITT = Intent-to-Treat  |  HR = Hazard Ratio  |  CI = Confidence Interval  ",
      "|  ONCVIZ-001 Synthetic Dataset"
    )
  ) +
  annotate(
    "label",
    x = 43, y = 0.04,
    label = paste0(
      "Log-rank:  ", pstr, "\n",
      sprintf("Median OS\n  Treatment : %.1f mo\n  Control   : %.1f mo",
              med_t, med_c)
    ),
    hjust = 1, vjust = 0,
    size = 3.3, family = "mono",
    fill = "white", color = C_TEXT,
    label.size = 0.4, label.padding = unit(0.4, "lines")
  ) +
  scale_color_manual(
    values = c("ARM=CONTROL" = C_CTRL, "ARM=TREATMENT" = C_TRT),
    labels = c(
      "ARM=CONTROL"   = sprintf("Control    (n=%d, events=13)", n_ctrl),
      "ARM=TREATMENT" = sprintf("Treatment  (n=%d, events=31)", n_trt)
    )
  ) +
  scale_fill_manual(
    values = c("ARM=CONTROL" = C_CTRL, "ARM=TREATMENT" = C_TRT)
  )

# Style risk table
km_plot$table <- km_plot$table +
  labs(title = "Number at risk", x = "Time (Months)") +
  theme(
    plot.title      = element_text(size = 9, color = C_MUT,
                                   face = "italic", hjust = 0),
    axis.text.y     = element_text(size = 9.5, face = "bold"),
    axis.title.x    = element_text(face = "bold", size = 11, color = C_TEXT),
    panel.grid      = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

png(file.path(OUTPUT_DIR, "km_overall_survival.png"),
    width = 13, height = 10.5, units = "in", res = 200)
print(km_plot)
dev.off()
cat("  Saved: km_overall_survival.png\n")

# ══════════════════════════════════════════════════════════════════════════════
#  2.  OVERALL SURVIVAL (OS) CURVE  — full clinical reporting standard
# ══════════════════════════════════════════════════════════════════════════════
cat("Building OS Curve...\n")

os_df <- adtte |>
  filter(PARAMCD == "OS") |>
  mutate(ARM = trimws(ARM))

os_fit <- survfit(
  Surv(AVALM, EVENT) ~ ARM,
  data = os_df,
  conf.type = "log-log"
)

hr_os  <- get_hr(os_df)
pstr_os <- fmt_p(hr_os$pval)

med_os_t <- summary(os_fit)$table["ARM=TREATMENT", "median"]
med_os_c <- summary(os_fit)$table["ARM=CONTROL",   "median"]

lm12_t  <- lm_rate(os_fit, 12)[2]   # TREATMENT is 2nd level
lm12_c  <- lm_rate(os_fit, 12)[1]
lm24_t  <- lm_rate(os_fit, 24)[2]
lm24_c  <- lm_rate(os_fit, 24)[1]

os_plot <- ggsurvplot(
  os_fit,
  data              = os_df,
  palette           = c(C_CTRL, C_TRT),
  size              = 1.1,
  conf.int          = TRUE,
  conf.int.alpha    = 0.10,
  censor.shape      = "|",
  censor.size       = 5,
  risk.table        = TRUE,
  risk.table.height = 0.22,
  risk.table.col    = "strata",
  risk.table.fontsize = 3.6,
  tables.theme      = theme_cleantable(),
  break.x.by        = 6,
  xlim              = c(0, 44),
  ylim              = c(0, 1.0),
  xlab              = "Time from Randomization (Months)",
  ylab              = "Probability of Overall Survival",
  legend.labs       = c(
    "Control    (n=18, events=13)",
    "Treatment  (n=62, events=31)"
  ),
  legend.title      = "",
  ggtheme           = theme_oncoviz()
)

# Landmark vertical guides
os_plot$plot <- os_plot$plot +
  geom_vline(xintercept = c(12, 24),
             color = "#D0D0D0", linewidth = 0.9, linetype = "dotted") +
  labs(
    title    = "Overall Survival",
    subtitle = paste0(
      "Phase I/II  \u00b7  Vizatinib 300 mg QD vs Control  \u00b7  ",
      "ITT Population  \u00b7  Data cut-off: 05 Mar 2026"
    ),
    caption  = paste0(
      "ITT = Intent-to-Treat  |  HR = Hazard Ratio  |  CI = Confidence Interval  ",
      "|  NR = Not Reached  |  ONCVIZ-001 Synthetic Dataset"
    )
  ) +
  # Stats inset (bottom-right)
  annotate(
    "label",
    x = 43, y = 0.04,
    label = paste0(
      sprintf("HR (95%% CI) : %.2f (%.2f\u2013%.2f)\n",
              hr_os$hr, hr_os$lo, hr_os$hi),
      "Log-rank    : ", pstr_os, "\n",
      "\nMedian OS\n",
      sprintf("  Treatment : %.1f mo\n", med_os_t),
      sprintf("  Control   : %.1f mo\n", med_os_c),
      "\n12-mo OS Rate\n",
      sprintf("  Treatment : %.1f%%\n", lm12_t),
      sprintf("  Control   : %.1f%%\n", lm12_c),
      "\n24-mo OS Rate\n",
      sprintf("  Treatment : %.1f%%\n", lm24_t),
      sprintf("  Control   : %.1f%%",   lm24_c)
    ),
    hjust = 1, vjust = 0,
    size = 3.1, family = "mono",
    fill = "#F9F9F9", color = C_TEXT,
    label.size = 0.4, label.padding = unit(0.45, "lines")
  ) +
  scale_color_manual(
    values = c("ARM=CONTROL" = C_CTRL, "ARM=TREATMENT" = C_TRT)
  ) +
  scale_fill_manual(
    values = c("ARM=CONTROL" = C_CTRL, "ARM=TREATMENT" = C_TRT)
  )

os_plot$table <- os_plot$table +
  labs(title = "Number at risk", x = "Time from Randomization (Months)") +
  theme(
    plot.title      = element_text(size = 9, color = C_MUT,
                                   face = "italic", hjust = 0),
    axis.text.y     = element_text(size = 9.5, face = "bold"),
    axis.title.x    = element_text(face = "bold", size = 11, color = C_TEXT),
    panel.grid      = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

png(file.path(OUTPUT_DIR, "os_curve.png"),
    width = 13, height = 10.5, units = "in", res = 200)
print(os_plot)
dev.off()
cat("  Saved: os_curve.png\n")

# ══════════════════════════════════════════════════════════════════════════════
#  3.  PROGRESSION-FREE SURVIVAL (PFS) CURVE
# ══════════════════════════════════════════════════════════════════════════════
cat("Building PFS Curve...\n")

pfs_df <- adtte |>
  filter(PARAMCD == "PFS") |>
  mutate(ARM = trimws(ARM))

pfs_fit <- survfit(
  Surv(AVALM, EVENT) ~ ARM,
  data = pfs_df,
  conf.type = "log-log"
)

hr_pfs   <- get_hr(pfs_df)
pstr_pfs <- fmt_p(hr_pfs$pval)

med_pfs_t <- summary(pfs_fit)$table["ARM=TREATMENT", "median"]
med_pfs_c <- summary(pfs_fit)$table["ARM=CONTROL",   "median"]

lm3_t  <- lm_rate(pfs_fit, 3)[2];  lm3_c  <- lm_rate(pfs_fit, 3)[1]
lm6_t  <- lm_rate(pfs_fit, 6)[2];  lm6_c  <- lm_rate(pfs_fit, 6)[1]
lm12_t <- lm_rate(pfs_fit, 12)[2]; lm12_c <- lm_rate(pfs_fit, 12)[1]

e_pfs_t <- sum(pfs_df$ARM == "TREATMENT" & pfs_df$EVENT == 1)
e_pfs_c <- sum(pfs_df$ARM == "CONTROL"   & pfs_df$EVENT == 1)

pfs_plot <- ggsurvplot(
  pfs_fit,
  data              = pfs_df,
  palette           = c(C_CTRL, C_TRT),
  size              = 1.1,
  conf.int          = TRUE,
  conf.int.alpha    = 0.10,
  censor.shape      = "|",
  censor.size       = 5,
  risk.table        = TRUE,
  risk.table.height = 0.22,
  risk.table.col    = "strata",
  risk.table.fontsize = 3.4,
  tables.theme      = theme_cleantable(),
  break.x.by        = 3,
  xlim              = c(0, 43),
  ylim              = c(0, 1.0),
  xlab              = "Time from Randomization (Months)",
  ylab              = "Probability of Progression-Free Survival",
  legend.labs       = c(
    sprintf("Control    (n=18, events=%d)", e_pfs_c),
    sprintf("Treatment  (n=62, events=%d)", e_pfs_t)
  ),
  legend.title      = "",
  ggtheme           = theme_oncoviz()
)

pfs_plot$plot <- pfs_plot$plot +
  geom_vline(xintercept = c(3, 6, 12),
             color = "#D0D0D0", linewidth = 0.9, linetype = "dotted") +
  labs(
    title    = "Progression-Free Survival",
    subtitle = paste0(
      "Phase I/II  \u00b7  Vizatinib 300 mg QD vs Control  \u00b7  ",
      "ITT Population  \u00b7  Data cut-off: 05 Mar 2026"
    ),
    caption  = paste0(
      "ITT = Intent-to-Treat  |  PFS = Progression-Free Survival  |  ",
      "HR = Hazard Ratio  |  CI = Confidence Interval  |  ONCVIZ-001 Synthetic Dataset"
    )
  ) +
  annotate(
    "label",
    x = 42, y = 0.04,
    label = paste0(
      sprintf("HR (95%% CI) : %.2f (%.2f\u2013%.2f)\n",
              hr_pfs$hr, hr_pfs$lo, hr_pfs$hi),
      "Log-rank    : ", pstr_pfs, "\n",
      "\nMedian PFS\n",
      sprintf("  Treatment : %.1f mo\n", med_pfs_t),
      sprintf("  Control   : %.1f mo\n", med_pfs_c),
      "\n3-mo PFS Rate\n",
      sprintf("  Treatment : %.1f%%\n",  lm3_t),
      sprintf("  Control   : %.1f%%\n",  lm3_c),
      "\n6-mo PFS Rate\n",
      sprintf("  Treatment : %.1f%%\n",  lm6_t),
      sprintf("  Control   : %.1f%%\n",  lm6_c),
      "\n12-mo PFS Rate\n",
      sprintf("  Treatment : %.1f%%\n",  lm12_t),
      sprintf("  Control   : %.1f%%",    lm12_c)
    ),
    hjust = 1, vjust = 0,
    size = 3.1, family = "mono",
    fill = "#F9F9F9", color = C_TEXT,
    label.size = 0.4, label.padding = unit(0.45, "lines")
  ) +
  scale_color_manual(
    values = c("ARM=CONTROL" = C_CTRL, "ARM=TREATMENT" = C_TRT)
  ) +
  scale_fill_manual(
    values = c("ARM=CONTROL" = C_CTRL, "ARM=TREATMENT" = C_TRT)
  )

pfs_plot$table <- pfs_plot$table +
  labs(title = "Number at risk", x = "Time from Randomization (Months)") +
  theme(
    plot.title      = element_text(size = 9, color = C_MUT,
                                   face = "italic", hjust = 0),
    axis.text.y     = element_text(size = 9.5, face = "bold"),
    axis.title.x    = element_text(face = "bold", size = 11, color = C_TEXT),
    panel.grid      = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

png(file.path(OUTPUT_DIR, "pfs_curve.png"),
    width = 13, height = 10.5, units = "in", res = 200)
print(pfs_plot)
dev.off()
cat("  Saved: pfs_curve.png\n")

# ══════════════════════════════════════════════════════════════════════════════
#  4.  EVENT-FREE SURVIVAL (EFS) CURVE
# ══════════════════════════════════════════════════════════════════════════════
cat("Building EFS Curve...\n")

efs_df <- adtte |>
  filter(PARAMCD == "EFS") |>
  mutate(ARM = trimws(ARM))

efs_fit <- survfit(
  Surv(AVALM, EVENT) ~ ARM,
  data = efs_df,
  conf.type = "log-log"
)

hr_efs   <- get_hr(efs_df)
pstr_efs <- fmt_p(hr_efs$pval)

med_efs_t <- summary(efs_fit)$table["ARM=TREATMENT", "median"]
med_efs_c <- summary(efs_fit)$table["ARM=CONTROL",   "median"]

lm3_efs_t  <- lm_rate(efs_fit, 3)[2];  lm3_efs_c  <- lm_rate(efs_fit, 3)[1]
lm6_efs_t  <- lm_rate(efs_fit, 6)[2];  lm6_efs_c  <- lm_rate(efs_fit, 6)[1]
lm12_efs_t <- lm_rate(efs_fit, 12)[2]; lm12_efs_c <- lm_rate(efs_fit, 12)[1]

e_efs_t    <- sum(efs_df$ARM == "TREATMENT" & efs_df$EVENT == 1)
e_efs_c    <- sum(efs_df$ARM == "CONTROL"   & efs_df$EVENT == 1)

# Event breakdown (EFS-specific)
n_prog  <- sum(efs_df$COMPTYPE == "PROGRESSION", na.rm = TRUE)
n_death <- sum(efs_df$COMPTYPE == "DEATH_WITHOUT_PROGRESSION", na.rm = TRUE)

efs_plot <- ggsurvplot(
  efs_fit,
  data              = efs_df,
  palette           = c(C_CTRL, C_TRT),
  size              = 1.1,
  conf.int          = TRUE,
  conf.int.alpha    = 0.10,
  censor.shape      = "|",
  censor.size       = 5,
  risk.table        = TRUE,
  risk.table.height = 0.22,
  risk.table.col    = "strata",
  risk.table.fontsize = 3.4,
  tables.theme      = theme_cleantable(),
  break.x.by        = 3,
  xlim              = c(0, 43),
  ylim              = c(0, 1.0),
  xlab              = "Time from Randomization (Months)",
  ylab              = "Probability of Event-Free Survival",
  legend.labs       = c(
    sprintf("Control    (n=18, events=%d)", e_efs_c),
    sprintf("Treatment  (n=62, events=%d)", e_efs_t)
  ),
  legend.title      = "",
  ggtheme           = theme_oncoviz()
)

# Mark "Death without progression" patients with triangle
dwp_pts <- efs_df |> filter(COMPTYPE == "DEATH_WITHOUT_PROGRESSION")

# Get survival probability at each DWP patient's time
sf_trt <- data.frame(
  time = efs_fit$time[1:efs_fit$strata["ARM=TREATMENT"]],
  surv = efs_fit$surv[1:efs_fit$strata["ARM=TREATMENT"]]
)

efs_plot$plot <- efs_plot$plot +
  geom_vline(xintercept = c(3, 6, 12),
             color = "#D0D0D0", linewidth = 0.9, linetype = "dotted") +
  # DWP marker (triangle) – approximate y at event time
  {
    if (nrow(dwp_pts) > 0) {
      list(
        geom_point(
          data  = dwp_pts,
          aes(x = AVALM, y = 0.02),
          shape = 17, size = 3.5,
          color = "#6929C4", inherit.aes = FALSE
        ),
        annotate(
          "text",
          x = dwp_pts$AVALM[1] + 0.5, y = 0.06,
          label = "Death w/o\nprogression",
          size = 2.8, color = "#6929C4",
          fontface = "italic", hjust = 0
        )
      )
    }
  } +
  labs(
    title    = "Event-Free Survival",
    subtitle = paste0(
      "Phase I/II  \u00b7  Vizatinib 300 mg QD vs Control  \u00b7  ITT Population  ",
      "\u00b7  Data cut-off: 05 Mar 2026\n",
      "EFS events: disease progression, death from any cause, ",
      "or treatment discontinuation due to toxicity"
    ),
    caption  = paste0(
      "ITT = Intent-to-Treat  |  EFS = Event-Free Survival  |  HR = Hazard Ratio  ",
      "|  CI = Confidence Interval  |  ONCVIZ-001 Synthetic Dataset"
    )
  ) +
  annotate(
    "label",
    x = 42, y = 0.04,
    label = paste0(
      sprintf("HR (95%% CI) : %.2f (%.2f\u2013%.2f)\n",
              hr_efs$hr, hr_efs$lo, hr_efs$hi),
      "Log-rank    : ", pstr_efs, "\n",
      "\nMedian EFS\n",
      sprintf("  Treatment : %.1f mo\n", med_efs_t),
      sprintf("  Control   : %.1f mo\n", med_efs_c),
      "\n3-mo EFS Rate\n",
      sprintf("  Treatment : %.1f%%\n",  lm3_efs_t),
      sprintf("  Control   : %.1f%%\n",  lm3_efs_c),
      "\n6-mo / 12-mo EFS Rate\n",
      sprintf("  Treatment : %.1f%% / %.1f%%\n", lm6_efs_t,  lm12_efs_t),
      sprintf("  Control   : %.1f%% / %.1f%%\n", lm6_efs_c,  lm12_efs_c),
      "\nEvent breakdown\n",
      sprintf("  Progression            : %d\n", n_prog),
      sprintf("  Death w/o progression  : %d",   n_death)
    ),
    hjust = 1, vjust = 0,
    size = 3.0, family = "mono",
    fill = "#F9F9F9", color = C_TEXT,
    label.size = 0.4, label.padding = unit(0.45, "lines")
  ) +
  scale_color_manual(
    values = c("ARM=CONTROL" = C_CTRL, "ARM=TREATMENT" = C_TRT)
  ) +
  scale_fill_manual(
    values = c("ARM=CONTROL" = C_CTRL, "ARM=TREATMENT" = C_TRT)
  )

efs_plot$table <- efs_plot$table +
  labs(title = "Number at risk", x = "Time from Randomization (Months)") +
  theme(
    plot.title      = element_text(size = 9, color = C_MUT,
                                   face = "italic", hjust = 0),
    axis.text.y     = element_text(size = 9.5, face = "bold"),
    axis.title.x    = element_text(face = "bold", size = 11, color = C_TEXT),
    panel.grid      = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

png(file.path(OUTPUT_DIR, "efs_curve.png"),
    width = 13, height = 10.5, units = "in", res = 200)
print(efs_plot)
dev.off()
cat("  Saved: efs_curve.png\n")

# ══════════════════════════════════════════════════════════════════════════════
#  5.  DISEASE-FREE SURVIVAL (DFS) CURVE  — single arm, CR patients only
# ══════════════════════════════════════════════════════════════════════════════
cat("Building DFS Curve...\n")

dfs_raw <- adtte |> filter(PARAMCD == "DFS")
dfs_df  <- dfs_raw |>
  left_join(
    adsl |> select(USUBJID, PHASE, DOSELEVEL, BESTRSPC, PRIORSURG),
    by = "USUBJID"
  ) |>
  mutate(
    ARM      = trimws(ARM),
    TUMORTYPE = trimws(TUMORTYPE)
  )

# Single arm – Treatment / CR patients only
n_dfs      <- nrow(dfs_df)
n_events   <- sum(dfs_df$EVENT)
n_censored <- sum(dfs_df$CNSR)
n_surg     <- sum(trimws(dfs_df$PRIORSURG) == "Y", na.rm = TRUE)

dfs_fit <- survfit(
  Surv(AVALM, EVENT) ~ 1,
  data = dfs_df,
  conf.type = "log-log"
)

med_dfs    <- summary(dfs_fit)$table["median"]
lm6_dfs    <- lm_rate(dfs_fit, 6)[1]
lm12_dfs   <- lm_rate(dfs_fit, 12)[1]
lm24_dfs   <- lm_rate(dfs_fit, 24)[1]

# Tumour type colours matching Python output
tumor_cols <- c(
  HCC   = "#E74C3C",
  BRCA  = "#9B59B6",
  NSCLC = "#2980B9",
  CRC   = "#27AE60",
  PDAC  = "#E67E22"
)

# Build base DFS plot manually (ggsurvplot single-strata)
dfs_base <- ggsurvplot(
  dfs_fit,
  data              = dfs_df,
  palette           = C_DFS,
  size              = 1.2,
  conf.int          = TRUE,
  conf.int.alpha    = 0.12,
  censor.shape      = "|",
  censor.size       = 5,
  risk.table        = TRUE,
  risk.table.height = 0.15,
  risk.table.col    = "black",
  risk.table.fontsize = 3.8,
  tables.theme      = theme_cleantable(),
  break.x.by        = 3,
  xlim              = c(0, max(dfs_df$AVALM) + 3),
  ylim              = c(0, 1.0),
  xlab              = "Time from Complete Response Confirmation (Months)",
  ylab              = "Probability of Disease-Free Survival",
  legend            = "none",
  ggtheme           = theme_oncoviz()
)

# Add per-patient markers and tumour labels
dfs_base$plot <- dfs_base$plot +
  geom_vline(xintercept = c(6, 12, 24),
             color = "#D0D0D0", linewidth = 0.9, linetype = "dotted") +
  # Recurrence markers
  geom_point(
    data = dfs_df |> filter(EVENT == 1),
    aes(x = AVALM, y = 0.08),
    shape = 4, size = 4, stroke = 2,
    color = "#C0392B", inherit.aes = FALSE
  ) +
  # Censored markers
  geom_point(
    data = dfs_df |> filter(EVENT == 0),
    aes(x = AVALM, y = 0.08),
    shape = "|", size = 6,
    color = C_MUT, inherit.aes = FALSE
  ) +
  # Tumour type labels per patient
  geom_text(
    data = dfs_df,
    aes(x = AVALM + 0.3, y = 0.18, label = TUMORTYPE,
        color = TUMORTYPE),
    size = 2.7, fontface = "bold", hjust = 0,
    inherit.aes = FALSE
  ) +
  scale_color_manual(values = tumor_cols, guide = "none") +
  # Legend: KM line + CI + markers
  annotate("segment",
           x = 0.5, xend = 3, y = 0.97, yend = 0.97,
           color = C_DFS, linewidth = 1.2) +
  annotate("text",
           x = 3.3, y = 0.97,
           label = sprintf("Vizatinib 300 mg  (n=%d, CR patients)", n_dfs),
           size = 3.5, color = C_TEXT, hjust = 0) +
  annotate("point",
           x = 0.8, y = 0.91,
           shape = 4, size = 3.5, stroke = 1.8, color = "#C0392B") +
  annotate("text",
           x = 3.3, y = 0.91,
           label = sprintf("Recurrence / Progression  (n=%d)", n_events),
           size = 3.5, color = C_TEXT, hjust = 0) +
  annotate("text",
           x = 0.8, y = 0.85, label = "|",
           size = 4.5, color = C_MUT) +
  annotate("text",
           x = 3.3, y = 0.85,
           label = sprintf("Censored  (n=%d)", n_censored),
           size = 3.5, color = C_TEXT, hjust = 0) +
  labs(
    title    = "Disease-Free Survival",
    subtitle = paste0(
      "Phase I/II  \u00b7  Vizatinib 300 mg QD  \u00b7  ",
      "Complete Responders  \u00b7  Data cut-off: 05 Mar 2026\n",
      "DFS defined as time from CR confirmation to disease recurrence, ",
      "progression, or death from any cause"
    ),
    caption  = paste0(
      "DFS = Disease-Free Survival  |  CR = Complete Response  |  ",
      "Single-arm analysis  |  ONCVIZ-001 Synthetic Dataset"
    )
  ) +
  annotate(
    "label",
    x = max(dfs_df$AVALM) + 2.5, y = 0.04,
    label = paste0(
      "Vizatinib 300 mg  (CR patients)\n",
      sprintf("Median DFS      : %.1f mo\n",   med_dfs),
      sprintf("6-mo DFS rate   : %.1f%%\n",    lm6_dfs),
      sprintf("12-mo DFS rate  : %.1f%%\n",    lm12_dfs),
      sprintf("24-mo DFS rate  : %.1f%%\n\n",  lm24_dfs),
      "Patient breakdown\n",
      sprintf("  HCC : %d   BRCA : %d   NSCLC : %d\n",
              sum(dfs_df$TUMORTYPE == "HCC",   na.rm = TRUE),
              sum(dfs_df$TUMORTYPE == "BRCA",  na.rm = TRUE),
              sum(dfs_df$TUMORTYPE == "NSCLC", na.rm = TRUE)),
      sprintf("  Prior surgery : %d/%d\n\n", n_surg, n_dfs),
      "Note: Single-arm; no control arm\n",
      "HR analysis not applicable"
    ),
    hjust = 1, vjust = 0,
    size = 3.0, family = "mono",
    fill = "#F9F9F9", color = C_TEXT,
    label.size = 0.4, label.padding = unit(0.45, "lines")
  )

dfs_base$table <- dfs_base$table +
  labs(
    title = "Number at risk",
    x     = "Time from Complete Response Confirmation (Months)"
  ) +
  theme(
    plot.title      = element_text(size = 9, color = C_MUT,
                                   face = "italic", hjust = 0),
    axis.text.y     = element_text(size = 9.5, face = "bold",
                                   color = C_DFS),
    axis.title.x    = element_text(face = "bold", size = 11, color = C_TEXT),
    panel.grid      = element_blank(),
    plot.background = element_rect(fill = "white", color = NA)
  )

png(file.path(OUTPUT_DIR, "dfs_curve.png"),
    width = 13, height = 10.5, units = "in", res = 200)
print(dfs_base)
dev.off()
cat("  Saved: dfs_curve.png\n")

# ══════════════════════════════════════════════════════════════════════════════
cat("\n\u2714 All survival curves saved to:", OUTPUT_DIR, "\n")
cat("  km_overall_survival.png\n")
cat("  os_curve.png\n")
cat("  pfs_curve.png\n")
cat("  efs_curve.png\n")
cat("  dfs_curve.png\n")
# ══════════════════════════════════════════════════════════════════════════════
