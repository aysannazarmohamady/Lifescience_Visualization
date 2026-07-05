suppressPackageStartupMessages({
  library(survival)
  library(survminer)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(gridExtra)
  library(grid)
  library(scales)
  library(cmprsk)   # for CIF (competing risks)
  library(boot)     # for bootstrap RMST CI
})

#Paths
DATA_DIR   <- "./Data/V1"          # adjust if needed
OUTPUT_DIR <- "./Outputs"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

#Load data
adtte <- read.csv(file.path(DATA_DIR, "ADTTE.csv"), stringsAsFactors = FALSE)
adsl  <- read.csv(file.path(DATA_DIR, "ADSL.csv"),  stringsAsFactors = FALSE)

adtte$ARM   <- trimws(adtte$ARM)
adtte$EVENT <- 1L - as.integer(adtte$CNSR)

#Global palette
C_TRT    <- "#1A73E8"
C_CTRL   <- "#D44C36"
C_CR     <- "#2d6a4f"
C_PR     <- "#2196F3"
C_ALL    <- "#444444"
C_GRID   <- "#E8E8E8"
C_BG     <- "#F8F9FA"
C_TEXT   <- "#161616"
C_MUT    <- "#6F6F6F"
C_GREEN  <- "#27AE60"

#Shared ggplot2 theme
theme_oncoviz <- function(base_size = 11) {
  theme_classic(base_size = base_size) %+replace%
    theme(
      panel.grid.major   = element_line(color = C_GRID,  linewidth = 0.7),
      panel.grid.minor   = element_blank(),
      panel.background   = element_rect(fill = C_BG,    color = NA),
      plot.background    = element_rect(fill = "white", color = NA),
      axis.line          = element_line(color = "#CCCCCC"),
      axis.ticks         = element_blank(),
      axis.text          = element_text(color = C_TEXT, size = rel(0.95)),
      axis.title         = element_text(color = C_TEXT, face = "bold", size = rel(1.05)),
      axis.title.y       = element_text(margin = margin(r = 10)),
      plot.title         = element_text(face = "bold", size = rel(1.25),
                                        color = C_TEXT, margin = margin(b = 3)),
      plot.subtitle      = element_text(size = rel(0.82), color = C_MUT,
                                        margin = margin(b = 4)),
      plot.caption       = element_text(size = rel(0.65), color = "#AAAAAA", hjust = 0),
      plot.margin        = margin(14, 14, 6, 14),
      legend.background  = element_rect(fill = "white", color = "#DDDDDD", linewidth = 0.8),
      legend.key         = element_rect(fill = "white"),
      legend.text        = element_text(size = rel(0.92)),
      legend.title       = element_text(face = "bold", size = rel(0.92)),
      legend.margin      = margin(6, 8, 6, 8)
    )
}

#Helpers
get_hr <- function(df, time_col = "AVALM", event_col = "EVENT",
                   arm_col = "ARM", ref = "CONTROL") {
  df[[arm_col]] <- relevel(factor(df[[arm_col]]), ref = ref)
  fit  <- coxph(as.formula(
    paste0("Surv(", time_col, ",", event_col, ") ~ ", arm_col)), data = df)
  smry <- summary(fit)
  list(hr  = smry$conf.int[1, "exp(coef)"],
       lo  = smry$conf.int[1, "lower .95"],
       hi  = smry$conf.int[1, "upper .95"],
       pval= smry$sctest["pvalue"])
}

fmt_p <- function(p) {
  if (p < 0.001) "p < 0.001" else sprintf("p = %.4f", p)
}

lm_rate <- function(km_fit, t, strata_idx = NULL) {
  sv <- summary(km_fit, times = t, extend = TRUE)
  if (is.null(strata_idx)) round(sv$surv * 100, 1)
  else                      round(sv$surv[strata_idx] * 100, 1)
}

style_risk_table <- function(tbl) {
  tbl +
    labs(title = "Number at risk") +
    theme(
      plot.title      = element_text(size = 9, color = C_MUT, face = "italic", hjust = 0),
      axis.text.y     = element_text(size = 9.5, face = "bold"),
      axis.title.x    = element_text(face = "bold", size = 11, color = C_TEXT),
      panel.grid      = element_blank(),
      plot.background = element_rect(fill = "white", color = NA)
    )
}

#  IMAGE 1 — Time to Response (TTR)
#  Left:  KM curve (responders only, treatment vs control)
#  Right: Individual TTR horizontal bar chart
cat("Building Image 1: TTR...\n")

ttr_df <- adtte %>% filter(PARAMCD == "TTR")

n_trt_ttr <- sum(ttr_df$ARM == "TREATMENT")
n_ctl_ttr <- sum(ttr_df$ARM == "CONTROL")

km_ttr <- survfit(Surv(AVALM, EVENT) ~ ARM, data = ttr_df, conf.type = "log-log")

# Medians
med_trt_ttr <- summary(km_ttr)$table["ARM=TREATMENT", "median"]
med_ctl_ttr <- summary(km_ttr)$table["ARM=CONTROL",   "median"]

# Left: KM plot 
p_ttr_km <- ggsurvplot(
  km_ttr,
  data            = ttr_df,
  palette         = c(C_CTRL, C_TRT),
  size            = 1.1,
  conf.int        = TRUE,
  conf.int.alpha  = 0.10,
  censor.shape    = "|", censor.size = 5,
  risk.table      = FALSE,
  break.x.by     = 1,
  xlim            = c(0, 5),
  ylim            = c(0, 1.05),
  xlab            = "Time to Response (Months)",
  ylab            = "Probability of Not Yet Responding",
  legend.labs     = c("Control", "Treatment (Vizatinib)"),
  legend.title    = "",
  ggtheme         = theme_oncoviz()
)

# Add median lines + stats box
p_ttr_km$plot <- p_ttr_km$plot +
  geom_hline(yintercept = 0.5, color = "grey60", linetype = "dotted", linewidth = 0.8) +
  geom_vline(xintercept = med_trt_ttr, color = C_TRT,  linetype = "dashed",
             linewidth = 0.9, alpha = 0.6) +
  geom_vline(xintercept = med_ctl_ttr, color = C_CTRL, linetype = "dashed",
             linewidth = 0.9, alpha = 0.6) +
  annotate("label",
           x = med_trt_ttr + 0.6, y = 0.15,
           label = sprintf("0.5 mo  %.1f mo", med_trt_ttr),
           size = 3.2, color = C_CTRL, hjust = 0,
           fill = "white", label.size = 0.3) +
  annotate("label",
           x = 0.05, y = 0.47,
           label = sprintf(paste0(
             "Responder Summary\n\n",
             "Treatment  n = %d\n  Median TTR: %.1f mo\n\n",
             "Control    n = %d\n  Median TTR: %.1f mo"),
             n_trt_ttr, med_trt_ttr, n_ctl_ttr, med_ctl_ttr),
           hjust = 0, vjust = 1, size = 3.0,
           family = "mono", fill = "white",
           color = C_TEXT, label.size = 0.4,
           label.padding = unit(0.4, "lines")) +
  scale_color_manual(
    values = c("ARM=CONTROL" = C_CTRL, "ARM=TREATMENT" = C_TRT),
    labels = c("ARM=CONTROL" = "Control", "ARM=TREATMENT" = "Treatment (Vizatinib)")
  ) +
  scale_fill_manual(values = c("ARM=CONTROL" = C_CTRL, "ARM=TREATMENT" = C_TRT))

# Right: Individual TTR bars 
ttr_swim <- ttr_df %>%
  arrange(AVALM) %>%
  mutate(
    row_id = row_number(),
    col    = ifelse(ARM == "TREATMENT", C_TRT, C_CTRL)
  )

p_ttr_bar <- ggplot(ttr_swim, aes(x = AVALM, y = row_id, fill = ARM)) +
  geom_col(width = 0.75, alpha = 0.78, color = "white", linewidth = 0.3) +
  scale_fill_manual(values = c("TREATMENT" = C_TRT, "CONTROL" = C_CTRL),
                    labels = c("TREATMENT" = "Treatment", "CONTROL" = "Control"),
                    name   = "") +
  labs(title    = "Individual TTR\n(All Responders)",
       x        = "Weeks to First Response (Months)",
       y        = NULL) +
  theme_oncoviz() +
  theme(axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom") +
  xlim(0, NA)

# Combine 
png(file.path(OUTPUT_DIR, "Image1_TTR.png"),
    width = 18, height = 8, units = "in", res = 180)

grid.arrange(
  p_ttr_km$plot, p_ttr_bar,
  ncol  = 2,
  widths = c(1.1, 0.9),
  top   = textGrob(
    "Time to Response (TTR)\nKaplan-Meier Analysis \u2013 Responders Only",
    gp = gpar(fontsize = 15, fontface = "bold"))
)

dev.off()
cat("  Saved: Image1_TTR.png\n")

#  IMAGE 2 — Time to Progression (TTP)
#  KM curve with 95% CI + number-at-risk table
cat("Building Image 2: TTP...\n")

ttp_df  <- adtte %>% filter(PARAMCD == "TTP")
km_ttp  <- survfit(Surv(AVALM, EVENT) ~ ARM, data = ttp_df, conf.type = "log-log")
hr_ttp  <- get_hr(ttp_df)

n_trt_ttp <- sum(ttp_df$ARM == "TREATMENT")
n_ctl_ttp <- sum(ttp_df$ARM == "CONTROL")
med_trt_ttp <- summary(km_ttp)$table["ARM=TREATMENT", "median"]
med_ctl_ttp <- summary(km_ttp)$table["ARM=CONTROL",   "median"]

p_ttp <- ggsurvplot(
  km_ttp,
  data              = ttp_df,
  palette           = c(C_CTRL, C_TRT),
  size              = 1.2,
  conf.int          = TRUE,
  conf.int.alpha    = 0.10,
  censor.shape      = "|", censor.size = 5,
  risk.table        = TRUE,
  risk.table.height = 0.22,
  risk.table.col    = "strata",
  risk.table.fontsize = 3.6,
  tables.theme      = theme_cleantable(),
  break.x.by        = 5,
  xlim               = c(0, 42),
  ylim               = c(0, 1.05),
  xlab               = "Time (Months)",
  ylab               = "Probability of No Progression",
  legend.labs        = c("Control", "Treatment (Vizatinib)"),
  legend.title       = "",
  ggtheme            = theme_oncoviz()
)

p_ttp$plot <- p_ttp$plot +
  geom_hline(yintercept = 0.5, color = "grey70", linetype = "dotted", linewidth = 0.8) +
  geom_vline(xintercept = med_trt_ttp, color = C_TRT,  linetype = "dotted",
             linewidth = 0.9, alpha = 0.55) +
  geom_vline(xintercept = med_ctl_ttp, color = C_CTRL, linetype = "dotted",
             linewidth = 0.9, alpha = 0.55) +
  labs(title = "Time to Progression (TTP)",
       subtitle = "Kaplan-Meier Curve with 95% Confidence Intervals") +
  annotate("label",
           x = 41, y = 0.98,
           label = sprintf(paste0(
             "Efficacy Summary\n\n",
             "Treatment  (n=%d)\n  Median TTP: %.1f mo\n",
             "Control    (n=%d)\n  Median TTP: %.1f mo\n\n",
             "Log-rank %s"),
             n_trt_ttp, med_trt_ttp,
             n_ctl_ttp, med_ctl_ttp,
             fmt_p(hr_ttp$pval)),
           hjust = 1, vjust = 1, size = 3.1, family = "mono",
           fill = "white", color = C_TEXT,
           label.size = 0.4, label.padding = unit(0.45, "lines")) +
  scale_color_manual(
    values = c("ARM=CONTROL" = C_CTRL, "ARM=TREATMENT" = C_TRT),
    labels = c("ARM=CONTROL" = "Control", "ARM=TREATMENT" = "Treatment (Vizatinib)")
  ) +
  scale_fill_manual(values = c("ARM=CONTROL" = C_CTRL, "ARM=TREATMENT" = C_TRT))

p_ttp$table <- style_risk_table(p_ttp$table) +
  labs(x = "Time (Months)")

png(file.path(OUTPUT_DIR, "Image2_TTP.png"),
    width = 13, height = 11, units = "in", res = 180)
print(p_ttp)
dev.off()
cat("  Saved: Image2_TTP.png\n")

#  IMAGE 3 — Duration of Response (DOR)
#  Left:  KM by CR / PR / Overall (treatment responders)
#  Right: Individual DOR swimmer plot
cat("Building Image 3: DOR...\n")

dor_df <- adtte %>%
  filter(PARAMCD == "DOR", ARM == "TREATMENT") %>%
  left_join(adsl %>% select(USUBJID, BESTRSPC), by = "USUBJID")

cr_df  <- dor_df %>% filter(BESTRSPC == "CR")
pr_df  <- dor_df %>% filter(BESTRSPC == "PR")
all_df <- dor_df %>% filter(BESTRSPC %in% c("CR","PR"))

km_cr  <- survfit(Surv(AVALM, EVENT) ~ 1, data = cr_df,  conf.type = "log-log")
km_pr  <- survfit(Surv(AVALM, EVENT) ~ 1, data = pr_df,  conf.type = "log-log")
km_dor <- survfit(Surv(AVALM, EVENT) ~ 1, data = all_df, conf.type = "log-log")

med_cr  <- summary(km_cr)$table["median"]
med_pr  <- summary(km_pr)$table["median"]
med_all <- summary(km_dor)$table["median"]

n_cr <- nrow(cr_df); n_pr <- nrow(pr_df); n_all <- nrow(all_df)

# Build combined KM data frame for manual ggplot
build_km_df <- function(km, label, col) {
  sv <- data.frame(
    time  = km$time,
    surv  = km$surv,
    lower = km$lower,
    upper = km$upper,
    label = label,
    col   = col
  )
  rbind(data.frame(time=0, surv=1, lower=1, upper=1, label=label, col=col), sv)
}

km_data <- bind_rows(
  build_km_df(km_cr,  sprintf("CR (n=%d)",      n_cr),  C_CR),
  build_km_df(km_pr,  sprintf("PR (n=%d)",      n_pr),  C_PR),
  build_km_df(km_dor, sprintf("Overall (n=%d)", n_all), C_ALL)
)
km_data$label <- factor(km_data$label,
                         levels = c(sprintf("CR (n=%d)", n_cr),
                                    sprintf("PR (n=%d)", n_pr),
                                    sprintf("Overall (n=%d)", n_all)))

col_map <- setNames(c(C_CR, C_PR, C_ALL), levels(km_data$label))
lty_map <- setNames(c("solid","solid","dashed"), levels(km_data$label))

p_dor_km <- ggplot(km_data, aes(x = time, y = surv,
                                  color = label, fill = label,
                                  linetype = label)) +
  geom_step(linewidth = 2.0) +
  geom_stepribbon(aes(ymin = lower, ymax = upper), alpha = 0.10,
                  color = NA, show.legend = FALSE) +
  geom_hline(yintercept = 0.5, color = "grey70", linetype = "dotted", linewidth = 0.7) +
  geom_vline(xintercept = med_cr,  color = C_CR,  linetype = "dotted",
             linewidth = 0.9, alpha = 0.5) +
  geom_vline(xintercept = med_pr,  color = C_PR,  linetype = "dotted",
             linewidth = 0.9, alpha = 0.5) +
  geom_vline(xintercept = med_all, color = C_ALL, linetype = "dotted",
             linewidth = 0.8, alpha = 0.5) +
  annotate("label",
           x = 0.3, y = 0.50,
           label = sprintf(paste0(
             "Response Duration Stats\n\n",
             "CR (n=%d)\n  Median: %.1f mo\n",
             "PR (n=%d)\n  Median: %.1f mo\n\n",
             "Overall median: %.1f mo"),
             n_cr, med_cr, n_pr, med_pr, med_all),
           hjust = 0, vjust = 1, size = 3.0, family = "mono",
           fill = "#f0f7f0", color = C_TEXT,
           label.size = 0.4, label.padding = unit(0.45, "lines")) +
  scale_color_manual(values = col_map, name = "") +
  scale_fill_manual(values  = col_map, name = "") +
  scale_linetype_manual(values = lty_map, name = "") +
  scale_y_continuous(labels = percent_format(accuracy = 1), limits = c(0, 1.05)) +
  scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
  labs(title    = "Duration of Response (DOR)",
       subtitle = "Kaplan-Meier by Best Overall Response",
       x        = "Duration of Response (Months)",
       y        = "Probability of Maintaining Response") +
  theme_oncoviz() +
  theme(legend.position = c(0.97, 0.97),
        legend.justification = c(1, 1))

# Right: Individual swimmer
swim_df <- all_df %>%
  arrange(AVALM) %>%
  mutate(
    row_id  = row_number(),
    label   = BESTRSPC,
    is_event= EVENT == 1
  )

p_dor_swim <- ggplot(swim_df, aes(x = AVALM, y = row_id, fill = label)) +
  geom_col(width = 0.75, alpha = 0.78, color = "white", linewidth = 0.2) +
  # Event marker (circle = event, triangle = censored)
  geom_point(data = swim_df %>% filter(EVENT == 1),
             aes(x = AVALM, y = row_id, shape = "Event"),
             color = "grey20", size = 4, fill = "grey20",
             inherit.aes = FALSE) +
  geom_point(data = swim_df %>% filter(EVENT == 0),
             aes(x = AVALM, y = row_id, shape = "Censored"),
             color = "grey20", size = 4,
             inherit.aes = FALSE) +
  scale_fill_manual(values = c("CR" = C_CR, "PR" = C_PR),
                    name   = "") +
  scale_shape_manual(values = c("Event" = 16, "Censored" = 17),
                     name   = "") +
  labs(title = "Individual DOR\n(Responders)",
       x     = "Duration of Response (Months)",
       y     = NULL) +
  theme_oncoviz() +
  theme(axis.text.y  = element_blank(),
        axis.ticks.y = element_blank(),
        legend.position = "bottom") +
  xlim(0, NA)

# Combine 
png(file.path(OUTPUT_DIR, "Image3_DOR.png"),
    width = 18, height = 9, units = "in", res = 180)
grid.arrange(
  p_dor_km, p_dor_swim,
  ncol   = 2,
  widths = c(1.1, 0.9)
)
dev.off()
cat("  Saved: Image3_DOR.png\n")

#  IMAGE 4 — Landmark Analysis: Overall Survival
#  Top:          KM with landmark lines at 6 / 12 / 24 months
#  Bottom-left:  Survival rate bar chart at each landmark
#  Bottom-right: Absolute difference with 95% bootstrap CI

cat("Building Image 4: Landmark OS...\n")

os_df   <- adtte %>% filter(PARAMCD == "OS")
os_trt  <- os_df %>% filter(ARM == "TREATMENT")
os_ctrl <- os_df %>% filter(ARM == "CONTROL")
km_os   <- survfit(Surv(AVALM, EVENT) ~ ARM, data = os_df, conf.type = "log-log")

# Survival rates at landmarks
lm_times <- c(6, 12, 24)
lm_t <- sapply(lm_times, function(t) lm_rate(km_os, t, strata_idx = 2))  # TRT=2nd
lm_c <- sapply(lm_times, function(t) lm_rate(km_os, t, strata_idx = 1))  # CTL=1st
diffs <- lm_t - lm_c

# Bootstrap CI for differences
boot_surv_diff <- function(data_t, data_c, t_lm, R = 300) {
  f <- function(d, i) {
    st <- d[i[i <= nrow(data_t)],  ]
    sc <- d[i[i >  nrow(data_t)] - nrow(data_t), ]
    k1 <- survfit(Surv(AVALM, EVENT) ~ 1, data = st)
    k2 <- survfit(Surv(AVALM, EVENT) ~ 1, data = sc)
    sv1 <- summary(k1, times = t_lm, extend = TRUE)$surv * 100
    sv2 <- summary(k2, times = t_lm, extend = TRUE)$surv * 100
    if (length(sv1) == 0) sv1 <- 0
    if (length(sv2) == 0) sv2 <- 0
    sv1 - sv2
  }
  combined <- bind_rows(data_t, data_c)
  idx_pool <- c(seq_len(nrow(data_t)), seq_len(nrow(data_c)) + nrow(data_t))
  boot_res <- boot(combined, statistic = f, R = R,
                   sim = "ordinary", stype = "i")
  boot.ci(boot_res, type = "perc", conf = 0.95)$percent[4:5]
}

set.seed(42)
ci_list <- lapply(lm_times, function(t)
  boot_surv_diff(os_trt, os_ctrl, t, R = 300))

ci_lo <- sapply(ci_list, `[`, 1)
ci_hi <- sapply(ci_list, `[`, 2)

# Top: KM with landmarks 
p_lm_km <- ggsurvplot(
  km_os,
  data              = os_df,
  palette           = c(C_CTRL, C_TRT),
  size              = 1.1,
  conf.int          = TRUE,
  conf.int.alpha    = 0.10,
  censor.shape      = "+", censor.size = 4,
  risk.table        = FALSE,
  break.x.by        = 5,
  xlim               = c(0, 46),
  ylim               = c(0, 1.05),
  xlab               = "Time (Months)",
  ylab               = "Overall Survival Probability",
  legend.labs        = c("Control", "Treatment (Vizatinib)"),
  legend.title       = "",
  ggtheme            = theme_oncoviz()
)

lm_cols   <- c("#E6A817", "#8E44AD", "#1ABC9C")
lm_labels <- c("6m", "12m", "24m")

p_lm_km$plot <- p_lm_km$plot +
  geom_vline(xintercept = lm_times,
             color    = lm_cols,
             linetype = "dashed", linewidth = 1.4, alpha = 0.85) +
  annotate("text",
           x = lm_times + 0.5, y = 1.03,
           label = lm_labels, color = lm_cols,
           fontface = "bold", size = 3.8) +
  scale_color_manual(
    values = c("ARM=CONTROL" = C_CTRL, "ARM=TREATMENT" = C_TRT),
    labels = c("ARM=CONTROL" = "Control", "ARM=TREATMENT" = "Treatment (Vizatinib)")
  ) +
  scale_fill_manual(values = c("ARM=CONTROL" = C_CTRL, "ARM=TREATMENT" = C_TRT))

# Bottom-left: Bar chart 
bar_df <- data.frame(
  Landmark = rep(paste0(lm_times, "-Month\nLandmark"), each = 2),
  Group    = rep(c("Treatment", "Control"), 3),
  Rate     = c(rbind(lm_t, lm_c))
)
bar_df$Landmark <- factor(bar_df$Landmark,
                           levels = paste0(lm_times, "-Month\nLandmark"))
bar_df$Group    <- factor(bar_df$Group, levels = c("Treatment", "Control"))

p_bar <- ggplot(bar_df, aes(x = Landmark, y = Rate, fill = Group)) +
  geom_col(position = position_dodge(width = 0.7),
           width = 0.65, alpha = 0.82) +
  geom_text(aes(label = paste0(round(Rate), "%"),
                color = Group),
            position = position_dodge(width = 0.7),
            vjust = -0.5, fontface = "bold", size = 3.5) +
  scale_fill_manual(values  = c("Treatment" = C_TRT, "Control" = C_CTRL)) +
  scale_color_manual(values = c("Treatment" = C_TRT, "Control" = C_CTRL),
                     guide = "none") +
  scale_y_continuous(labels = function(x) paste0(x, "%"), limits = c(0, 110)) +
  labs(title  = "Survival Rate at Each Landmark",
       x      = NULL, y = "Survival Probability",
       fill   = "") +
  theme_oncoviz() +
  theme(legend.position = "top")

# Bottom-right: Difference + CI 
diff_df <- data.frame(
  Landmark = factor(paste0(lm_times, "-Month\nLandmark"),
                    levels = paste0(lm_times, "-Month\nLandmark")),
  diff  = diffs,
  ci_lo = ci_lo,
  ci_hi = ci_hi
)

p_diff <- ggplot(diff_df, aes(x = Landmark, y = diff)) +
  geom_col(fill = C_GREEN, alpha = 0.82, width = 0.5) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.18, color = "black", linewidth = 1.2) +
  geom_text(aes(label = sprintf("+%.0f%%", diff), y = ci_hi + 1.5),
            fontface = "bold", size = 3.5, color = C_GREEN) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.9) +
  scale_y_continuous(labels = function(x) sprintf("%+.0f%%", x)) +
  labs(title = "Absolute Difference with 95% CI",
       x     = NULL,
       y     = "Survival Difference\n(Treatment \u2212 Control)") +
  theme_oncoviz()

# Combine 
png(file.path(OUTPUT_DIR, "Image4_Landmark_OS.png"),
    width = 14, height = 13, units = "in", res = 180)

grid.arrange(
  p_lm_km$plot,
  arrangeGrob(p_bar, p_diff, ncol = 2),
  nrow         = 2,
  heights      = c(1.6, 1),
  top          = textGrob(
    "Landmark Analysis \u2013 Overall Survival\nSurvival Rates at Pre-Specified Landmark Time Points (6, 12, 24 months)",
    gp = gpar(fontsize = 14, fontface = "bold"))
)

dev.off()
cat("  Saved: Image4_Landmark_OS.png\n")


#  IMAGE 5 — Competing Risks CIF
#  Aalen-Johansen CIF for Progression vs Death without Progression
#  With 95% bootstrap CI — one panel per arm

cat("Building Image 5: Competing Risks CIF...\n")

os_cmp <- adtte %>%
  filter(PARAMCD == "OS") %>%
  mutate(
    # etype: 1=Progression (event of interest), 2=Death w/o Progression, 0=Censored
    etype = case_when(
      COMPTYPE == "PROGRESSION"               ~ 1L,
      COMPTYPE == "DEATH_WITHOUT_PROGRESSION" ~ 2L,
      TRUE                                    ~ 0L
    )
  )

trt_cmp  <- os_cmp %>% filter(ARM == "TREATMENT")
ctrl_cmp <- os_cmp %>% filter(ARM == "CONTROL")

# Fit CIF via cmprsk::cuminc
fit_cif_trt  <- cuminc(trt_cmp$AVALM,  trt_cmp$etype)
fit_cif_ctrl <- cuminc(ctrl_cmp$AVALM, ctrl_cmp$etype)

# Extract CIF curves into data frames
extract_cif <- function(fit, group_name, event_label, event_code) {
  comp_name <- paste(group_name, event_code)
  if (!comp_name %in% names(fit)) return(NULL)
  data.frame(
    time  = fit[[comp_name]]$time,
    est   = fit[[comp_name]]$est * 100,
    var   = fit[[comp_name]]$var,
    event = event_label
  )
}

# Bootstrap CI for CIF
boot_cif <- function(df, etype_col = "etype", time_col = "AVALM",
                     n_boot = 300, seed = 42, tgrid = seq(0, 42, length.out = 200)) {
  set.seed(seed)
  mat1 <- matrix(NA, n_boot, length(tgrid))
  mat2 <- matrix(NA, n_boot, length(tgrid))
  for (b in seq_len(n_boot)) {
    samp <- df[sample(nrow(df), replace = TRUE), ]
    fit  <- cuminc(samp[[time_col]], samp[[etype_col]])
    nms  <- names(fit)
    nm1  <- grep("1$", nms, value = TRUE)[1]
    nm2  <- grep("2$", nms, value = TRUE)[1]
    if (!is.na(nm1))
      mat1[b,] <- approx(fit[[nm1]]$time, fit[[nm1]]$est,
                         xout = tgrid, rule = 2, method = "constant")$y * 100
    if (!is.na(nm2))
      mat2[b,] <- approx(fit[[nm2]]$time, fit[[nm2]]$est,
                         xout = tgrid, rule = 2, method = "constant")$y * 100
  }
  list(tgrid = tgrid,
       lo1   = apply(mat1, 2, quantile, 0.025, na.rm = TRUE),
       hi1   = apply(mat1, 2, quantile, 0.975, na.rm = TRUE),
       lo2   = apply(mat2, 2, quantile, 0.025, na.rm = TRUE),
       hi2   = apply(mat2, 2, quantile, 0.975, na.rm = TRUE))
}

boot_trt  <- boot_cif(trt_cmp,  n_boot = 300)
boot_ctrl <- boot_cif(ctrl_cmp, n_boot = 300)

# Helper: landmark value from CIF
cif_at <- function(fit, group, ecode, t_lm) {
  comp_name <- paste(group, ecode)
  if (!comp_name %in% names(fit)) return(0)
  approx(fit[[comp_name]]$time, fit[[comp_name]]$est * 100,
         xout = t_lm, rule = 2, method = "constant")$y
}

# Build one CIF panel as ggplot
make_cif_panel <- function(fit, boot_res, group, n_pts, title, box_fill) {
  C_PROG  <- "#1565C0"
  C_DEATH <- "#C62828"

  # Progression curve
  nms <- names(fit)
  nm1 <- grep("1$", nms, value = TRUE)[1]
  nm2 <- grep("2$", nms, value = TRUE)[1]

  prog_df  <- data.frame(time = fit[[nm1]]$time, est = fit[[nm1]]$est * 100)
  death_df <- data.frame(time = fit[[nm2]]$time, est = fit[[nm2]]$est * 100)

  ci_df1 <- data.frame(time  = boot_res$tgrid,
                       lo    = boot_res$lo1,
                       hi    = boot_res$hi1)
  ci_df2 <- data.frame(time  = boot_res$tgrid,
                       lo    = boot_res$lo2,
                       hi    = boot_res$hi2)

  lm12p <- cif_at(fit, group, 1, 12); lm12d <- cif_at(fit, group, 2, 12)
  lm24p <- cif_at(fit, group, 1, 24); lm24d <- cif_at(fit, group, 2, 24)

  box_label <- sprintf(paste0(
    "Cumulative Incidence\n",
    "\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n",
    "At 12 months:\n",
    "  Progression:  %.1f%%\n",
    "  Death (comp): %.1f%%\n",
    "At 24 months:\n",
    "  Progression:  %.1f%%\n",
    "  Death (comp): %.1f%%\n",
    "\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\u2500\n",
    "  n = %d"),
    lm12p, lm12d, lm24p, lm24d, n_pts)

  ggplot() +
    # CI ribbons
    geom_ribbon(data = ci_df1, aes(x = time, ymin = lo, ymax = hi),
                fill = C_PROG, alpha = 0.11) +
    geom_ribbon(data = ci_df2, aes(x = time, ymin = lo, ymax = hi),
                fill = C_DEATH, alpha = 0.11) +
    # Step curves
    geom_step(data = prog_df,  aes(x = time, y = est,
                                    color = "Progression (event of interest)"),
              linewidth = 2.1, direction = "hv") +
    geom_step(data = death_df, aes(x = time, y = est,
                                    color = "Death w/o Progression (competing)"),
              linewidth = 2.1, direction = "hv") +
    annotate("label",
             x = Inf, y = 42,
             label     = box_label,
             hjust = 1, vjust = 1, size = 2.9, family = "mono",
             fill      = box_fill, color = C_TEXT,
             label.size = 0.4, label.padding = unit(0.4, "lines")) +
    scale_color_manual(
      values = c("Progression (event of interest)"    = C_PROG,
                 "Death w/o Progression (competing)" = C_DEATH),
      name   = "") +
    scale_y_continuous(labels = function(x) paste0(x, "%"),
                       limits = c(0, 105)) +
    scale_x_continuous(expand = expansion(mult = c(0, 0.02))) +
    labs(title = title,
         x     = "Time (Months)",
         y     = "Cumulative Incidence") +
    theme_oncoviz() +
    theme(legend.position = "top",
          legend.direction = "vertical")
}

p_cif_trt  <- make_cif_panel(fit_cif_trt,  boot_trt,  "1", nrow(trt_cmp),
                               "Competing Risks \u2013 Treatment (Vizatinib)\n(Cumulative Incidence)",
                               "#e3f0fb")
p_cif_ctrl <- make_cif_panel(fit_cif_ctrl, boot_ctrl, "1", nrow(ctrl_cmp),
                               "Competing Risks \u2013 Control\n(Cumulative Incidence)",
                               "#fce8e8")

png(file.path(OUTPUT_DIR, "Image5_CIF.png"),
    width = 16, height = 7, units = "in", res = 180)

grid.arrange(
  p_cif_trt, p_cif_ctrl,
  ncol = 2,
  top  = textGrob(
    "Competing Risks Analysis\nCumulative Incidence Functions (CIF) with 95% Bootstrap CI",
    gp = gpar(fontsize = 14, fontface = "bold"))
)

dev.off()
cat("  Saved: Image5_CIF.png\n")

#  IMAGE 6 — Restricted Mean Survival Time (RMST)
#  Top:          KM with shaded RMST area at τ = 24 months
#  Bottom-left:  RMST vs restriction time τ (12, 18, 24, 30)
#  Bottom-right: RMST difference with 95% bootstrap CI

cat("Building Image 6: RMST...\n")

os_df2  <- adtte %>% filter(PARAMCD == "OS")
os_trt2 <- os_df2 %>% filter(ARM == "TREATMENT")
os_ctl2 <- os_df2 %>% filter(ARM == "CONTROL")

km_ot <- survfit(Surv(AVALM, EVENT) ~ 1, data = os_trt2, conf.type = "log-log")
km_oc <- survfit(Surv(AVALM, EVENT) ~ 1, data = os_ctl2, conf.type = "log-log")

# RMST via numerical integration (trapezoid)
calc_rmst <- function(km, tau) {
  t <- c(0, km$time[km$time <= tau], tau)
  s <- c(1, km$surv[km$time <= tau],
         approx(km$time, km$surv, xout = tau, rule = 2)$y)
  # remove duplicates
  ord <- order(t); t <- t[ord]; s <- s[ord]
  dup <- duplicated(t); t <- t[!dup]; s <- s[!dup]
  sum(diff(t) * s[-length(s)])   # left-Riemann (step function)
}

taus     <- c(12, 18, 24, 30)
TAU      <- 24

rmst_t_vals <- sapply(taus, function(tau) calc_rmst(km_ot, tau))
rmst_c_vals <- sapply(taus, function(tau) calc_rmst(km_oc, tau))
rmst_diffs  <- rmst_t_vals - rmst_c_vals

# Bootstrap CI for RMST difference
boot_rmst_diff <- function(data_t, data_c, tau, R = 300, seed = 42) {
  set.seed(seed)
  diffs_b <- numeric(R)
  for (b in seq_len(R)) {
    st <- data_t[sample(nrow(data_t), replace = TRUE), ]
    sc <- data_c[sample(nrow(data_c), replace = TRUE), ]
    kt <- survfit(Surv(AVALM, EVENT) ~ 1, data = st)
    kc <- survfit(Surv(AVALM, EVENT) ~ 1, data = sc)
    diffs_b[b] <- calc_rmst(kt, tau) - calc_rmst(kc, tau)
  }
  quantile(diffs_b, c(0.025, 0.975))
}

set.seed(42)
ci_rmst <- lapply(taus, function(tau)
  boot_rmst_diff(os_trt2, os_ctl2, tau, R = 300))

ci_lo_rmst <- sapply(ci_rmst, `[`, 1)
ci_hi_rmst <- sapply(ci_rmst, `[`, 2)

r_t24 <- calc_rmst(km_ot, TAU)
r_c24 <- calc_rmst(km_oc, TAU)

# Top: KM with shaded RMST region 
km_os2  <- survfit(Surv(AVALM, EVENT) ~ ARM, data = os_df2, conf.type = "log-log")

p_rmst_km <- ggsurvplot(
  km_os2,
  data              = os_df2,
  palette           = c(C_CTRL, C_TRT),
  size              = 1.1,
  conf.int          = TRUE,
  conf.int.alpha    = 0.10,
  censor.shape      = "+", censor.size = 4,
  risk.table        = FALSE,
  break.x.by        = 5,
  xlim               = c(0, 46),
  ylim               = c(0, 1.05),
  xlab               = "Time (Months)",
  ylab               = "Overall Survival Probability",
  legend.labs        = c("Control", "Treatment (Vizatinib)"),
  legend.title       = "",
  ggtheme            = theme_oncoviz()
)

# Shade RMST area under treatment curve
shade_t <- data.frame(
  time = c(0, km_ot$time[km_ot$time <= TAU], TAU),
  surv = c(1, km_ot$surv[km_ot$time <= TAU],
           approx(km_ot$time, km_ot$surv, xout = TAU, rule = 2)$y)
)

p_rmst_km$plot <- p_rmst_km$plot +
  geom_ribbon(data  = shade_t,
              aes(x = time, ymin = 0, ymax = surv),
              fill  = C_TRT, alpha = 0.10, inherit.aes = FALSE) +
  geom_vline(xintercept = TAU, color = "#E67E22",
             linetype = "dashed", linewidth = 1.8, alpha = 0.9) +
  annotate("text", x = TAU + 0.5, y = 1.03,
           label = sprintf("\u03c4 = %d mo", TAU),
           color = "#E67E22", fontface = "bold", size = 3.8) +
  annotate("label",
           x = 1, y = 0.40,
           label = sprintf(paste0(
             "RMST at \u03c4=%dm\n\n",
             "Treatment: %.1f mo\n",
             "Control:   %.1f mo\n",
             "Difference: +%.1f mo"),
             TAU, r_t24, r_c24, r_t24 - r_c24),
           hjust = 0, vjust = 1, size = 3.0, family = "mono",
           fill = "white", color = C_TEXT,
           label.size = 0.4, label.padding = unit(0.4, "lines")) +
  scale_color_manual(
    values = c("ARM=CONTROL" = C_CTRL, "ARM=TREATMENT" = C_TRT),
    labels = c("ARM=CONTROL" = "Control", "ARM=TREATMENT" = "Treatment (Vizatinib)")
  ) +
  scale_fill_manual(values = c("ARM=CONTROL" = C_CTRL, "ARM=TREATMENT" = C_TRT))

# Bottom-left: RMST vs tau 
rmst_df <- data.frame(
  tau    = taus,
  trt    = rmst_t_vals,
  ctrl   = rmst_c_vals
)

p_rmst_line <- ggplot(rmst_df, aes(x = tau)) +
  geom_abline(slope = 1, intercept = 0, color = "grey60",
              linetype = "dashed", linewidth = 0.9) +
  geom_line(aes(y = trt,  color = "Treatment (Vizatinib)"),
            linewidth = 1.8) +
  geom_point(aes(y = trt, color = "Treatment (Vizatinib)"),
             size = 4) +
  geom_line(aes(y = ctrl, color = "Control"),
            linewidth = 1.8) +
  geom_point(aes(y = ctrl, color = "Control"),
             size = 4) +
  geom_text(aes(y = trt  + 0.5, label = sprintf("%.1f", trt)),
            size = 3.2, color = C_TRT, fontface = "bold", hjust = -0.2) +
  geom_text(aes(y = ctrl - 1.3, label = sprintf("%.1f", ctrl)),
            size = 3.2, color = C_CTRL, fontface = "bold", hjust = -0.2) +
  scale_color_manual(
    values = c("Treatment (Vizatinib)" = C_TRT, "Control" = C_CTRL),
    name   = "") +
  scale_x_continuous(limits = c(11, 31)) +
  labs(title    = "RMST vs. Restriction Time\nwith 95% Bootstrap CI",
       x        = "Restriction Time \u03c4 (Months)",
       y        = "RMST (Months)") +
  theme_oncoviz() +
  theme(legend.position = "top")

# Bottom-right: RMST difference 
rdiff_df <- data.frame(
  tau   = factor(taus),
  diff  = rmst_diffs,
  ci_lo = ci_lo_rmst,
  ci_hi = ci_hi_rmst
)

p_rmst_diff <- ggplot(rdiff_df, aes(x = tau, y = diff)) +
  geom_col(fill = C_GREEN, alpha = 0.82, width = 0.5) +
  geom_errorbar(aes(ymin = ci_lo, ymax = ci_hi),
                width = 0.18, color = "black", linewidth = 1.2) +
  geom_text(aes(label = sprintf("+%.1fmo", diff), y = ci_hi + 0.3),
            fontface = "bold", size = 3.5, color = C_GREEN) +
  geom_hline(yintercept = 0, color = "grey60", linewidth = 0.9) +
  labs(title = "RMST Difference with 95% Bootstrap CI",
       x     = "Restriction Time \u03c4 (Months)",
       y     = "RMST Difference (Months)\n(Treatment \u2212 Control)") +
  theme_oncoviz()

# Combine 
png(file.path(OUTPUT_DIR, "Image6_RMST.png"),
    width = 15, height = 13, units = "in", res = 180)

grid.arrange(
  p_rmst_km$plot,
  arrangeGrob(p_rmst_line, p_rmst_diff, ncol = 2),
  nrow    = 2,
  heights = c(1.6, 1),
  top     = textGrob(
    paste0("Restricted Mean Survival Time (RMST)\n",
           "Kaplan-Meier Curves with Shaded RMST Area (\u03c4 = 24 months)"),
    gp = gpar(fontsize = 14, fontface = "bold"))
)

dev.off()
cat("  Saved: Image6_RMST.png\n")

cat("\n\u2714 All 6 plots saved to:", OUTPUT_DIR, "\n")
cat("  Image1_TTR.png\n")
cat("  Image2_TTP.png\n")
cat("  Image3_DOR.png\n")
cat("  Image4_Landmark_OS.png\n")
cat("  Image5_CIF.png\n")
cat("  Image6_RMST.png\n")
