suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(survival)
  library(patchwork)
})

DATA_DIR   <- "./Data/V1"
OUTPUT_DIR <- "./Out"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

adtte <- read.csv(file.path(DATA_DIR, "ADTTE.csv"), stringsAsFactors = FALSE)
os <- adtte %>% filter(PARAMCD == "OS") %>%
  mutate(EVENT = 1L - CNSR, TRT = as.integer(ARM == "TREATMENT"))

cox_hr <- function(df) {
  if (length(unique(df$TRT)) < 2 || sum(df$EVENT) < 2) return(c(hr = NA, lo = NA, hi = NA, se = NA))
  fit <- tryCatch(coxph(Surv(AVAL, EVENT) ~ TRT, data = df), error = function(e) NULL)
  if (is.null(fit)) return(c(hr = NA, lo = NA, hi = NA, se = NA))
  ci <- exp(confint(fit))
  c(hr = exp(coef(fit)), lo = ci[1], hi = ci[2], se = sqrt(vcov(fit)[1, 1]))
}

subgroups <- list(
  "Overall"       = os,
  "NSCLC"         = filter(os, TUMORTYPE == "NSCLC"),
  "CRC"           = filter(os, TUMORTYPE == "CRC"),
  "BRCA"          = filter(os, TUMORTYPE == "BRCA"),
  "HCC"           = filter(os, TUMORTYPE == "HCC"),
  "PDAC"          = filter(os, TUMORTYPE == "PDAC"),
  "Age < 65"      = filter(os, AGEGR1 == "<65"),
  "Age >= 65"     = filter(os, AGEGR1 == ">=65"),
  "Male"          = filter(os, SEX == "M"),
  "Female"        = filter(os, SEX == "F"),
  "ECOG 0"        = filter(os, ECOG == 0),
  "ECOG 1"        = filter(os, ECOG == 1),
  "TMB-High"      = filter(os, TMBHIGH == "Y"),
  "TMB-Low"       = filter(os, TMBHIGH == "N"),
  "PD-L1 High"    = filter(os, PDL1GRP == "HIGH"),
  "PD-L1 Medium"  = filter(os, PDL1GRP == "MED")
)

fdf <- bind_rows(lapply(names(subgroups), function(lab) {
  d <- subgroups[[lab]]
  r <- unname(cox_hr(d))
  data.frame(label = lab, hr = r[1], lo = r[2], hi = r[3], se = r[4],
             n = sum(d$ARM %in% c("TREATMENT", "CONTROL")), overall = (lab == "Overall"),
             stringsAsFactors = FALSE)
}))

fdf <- fdf %>%
  mutate(w = ifelse(!overall & !is.na(se) & se > 0, 1 / se^2, NA)) %>%
  mutate(weight_pct = ifelse(!overall, 100 * w / sum(w, na.rm = TRUE), NA))

fdf$label <- factor(paste0(fdf$label, " (n=", fdf$n, ")"), levels = rev(paste0(fdf$label, " (n=", fdf$n, ")")))
fdf <- fdf %>% mutate(sig = !overall & !is.na(hi) & hi < 1,
                       lo_c = pmax(lo, 0.05), hi_c = pmin(hi, 10))

p <- ggplot(fdf, aes(y = label)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#888888") +
  geom_errorbarh(data = filter(fdf, !overall, !is.na(hr)),
                  aes(xmin = lo_c, xmax = hi_c), height = 0, linewidth = 0.6,
                  color = ifelse(filter(fdf, !overall, !is.na(hr))$sig, "#B2182B", "#333333")) +
  geom_point(data = filter(fdf, !overall, !is.na(hr)),
             aes(x = hr, size = weight_pct, color = sig)) +
  scale_color_manual(values = c(`TRUE` = "#B2182B", `FALSE` = "#333333"), guide = "none") +
  scale_size_continuous(range = c(2, 9), name = "Weight %") +
  geom_point(data = filter(fdf, overall), aes(x = hr), shape = 18, size = 7, color = "#08306B") +
  geom_errorbarh(data = filter(fdf, overall), aes(xmin = lo_c, xmax = hi_c), height = 0, color = "#08306B", linewidth = 0.8) +
  scale_x_log10(breaks = c(0.1, 0.25, 0.5, 1, 2, 4), limits = c(0.08, 6)) +
  labs(x = "Hazard Ratio (95% CI), log scale", y = NULL,
       title = "Forest Plot \u2014 Overall Survival by Subgroup\nTreatment vs. Control (Cox Proportional Hazards)") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        legend.position = "right")

ggsave(file.path(OUTPUT_DIR, "forest_plot.png"), p, width = 11, height = 8.3, dpi = 300)

# Companion table (HR/CI/weight) written alongside the figure
write.csv(fdf %>% select(label, n, hr, lo, hi, weight_pct),
          file.path(OUTPUT_DIR, "forest_plot_table.csv"), row.names = FALSE)
