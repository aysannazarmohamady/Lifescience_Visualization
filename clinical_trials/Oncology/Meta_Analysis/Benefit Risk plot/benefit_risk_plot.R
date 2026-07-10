suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

DATA_DIR   <- "./Data/V1"
OUTPUT_DIR <- "./Out"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

adsl <- read.csv(file.path(DATA_DIR, "ADSL.csv"), stringsAsFactors = FALSE)
adae <- read.csv(file.path(DATA_DIR, "ADAE.csv"), stringsAsFactors = FALSE)

wilson_ci <- function(k, n) {
  if (n == 0) return(c(p = NA, lo = NA, hi = NA))
  p <- k / n; z <- 1.96
  denom <- 1 + z^2 / n
  center <- (p + z^2 / (2 * n)) / denom
  half <- z * sqrt(p * (1 - p) / n + z^2 / (4 * n^2)) / denom
  c(p = p * 100, lo = max(0, center - half) * 100, hi = min(1, center + half) * 100)
}

groups <- c(unique(adsl$TUMORTYPE), "Overall Treatment", "Overall Control")
rows <- lapply(groups, function(g) {
  sub <- if (g == "Overall Treatment") filter(adsl, ARM == "TREATMENT")
         else if (g == "Overall Control") filter(adsl, ARM == "CONTROL")
         else filter(adsl, TUMORTYPE == g, ARM == "TREATMENT")
  n <- nrow(sub); if (n == 0) return(NULL)
  orr_k <- sum(sub$BESTRSPC %in% c("CR", "PR"))
  orr <- wilson_ci(orr_k, n)
  pts <- sub$USUBJID
  g3_k <- n_distinct(adae$USUBJID[adae$USUBJID %in% pts & adae$AETOXGR >= 3])
  g3 <- wilson_ci(g3_k, n)
  data.frame(group = g, n = n, orr = orr["p"], orr_lo = orr["lo"], orr_hi = orr["hi"],
             g3 = g3["p"], g3_lo = g3["lo"], g3_hi = g3["hi"])
})
df <- bind_rows(rows)
df$group <- factor(df$group, levels = df$group)

orr_med <- median(df$orr); g3_med <- median(df$g3)

p_main <- ggplot(df, aes(x = g3, y = orr, color = group)) +
  geom_vline(xintercept = g3_med, color = "#cccccc") +
  geom_hline(yintercept = orr_med, color = "#cccccc") +
  geom_errorbarh(aes(xmin = g3_lo, xmax = g3_hi), height = 0, linewidth = 0.9) +
  geom_errorbar(aes(ymin = orr_lo, ymax = orr_hi), width = 0, linewidth = 0.9) +
  geom_point(aes(shape = grepl("Overall", group), size = grepl("Overall", group))) +
  scale_shape_manual(values = c(`TRUE` = 15, `FALSE` = 16), guide = "none") +
  scale_size_manual(values = c(`TRUE` = 5, `FALSE` = 4), guide = "none") +
  annotate("text", x = min(df$g3_lo) + 2, y = max(df$orr_hi) - 2, label = "More favorable\n(higher benefit, lower risk)",
           color = "#2166AC", hjust = 0, vjust = 1, size = 3) +
  annotate("text", x = median(df$g3), y = 3, label = "Less favorable\n(lower benefit, higher risk)",
           color = "#B2182B", hjust = 0.5, vjust = 0, size = 3) +
  labs(x = "Risk: Grade \u2265 3 Adverse Event Rate (%, 95% Wilson CI)",
       y = "Benefit: Objective Response Rate, ORR (%, 95% Wilson CI)", color = NULL) +
  theme_classic(base_size = 12) +
  theme(legend.position = "right", legend.spacing.y = unit(0.3, "cm")) +
  guides(color = guide_legend(byrow = TRUE))

p <- p_main +
  labs(title = "Benefit-Risk Plot \u2014 ORR vs. Grade \u2265 3 AE Rate by Tumor Type\n(Quadrant Framework)",
       caption = "Benefit = ORR (CR+PR) per RECIST 1.1. Risk = patients with \u22651 Grade \u22653 AE (CTCAE). Dashed reference lines = cross-group medians.\nFramework: benefit-risk quadrant scatter (informal MCDA); no formal weighting applied \u2014 confirm framework choice (e.g. BRAT, PrOACT-URL) for regulatory use.") +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13.5),
        plot.caption = element_text(size = 7.6, face = "italic", hjust = 0.5))

ggsave(file.path(OUTPUT_DIR, "benefit_risk_plot.png"), p, width = 11.5, height = 8, dpi = 300)
