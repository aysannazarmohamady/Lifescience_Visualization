suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

DATA_DIR   <- "./Data/V1"
OUTPUT_DIR <- "./Out"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

adex <- read.csv(file.path(DATA_DIR, "ADEX.csv"), stringsAsFactors = FALSE)
adsl <- read.csv(file.path(DATA_DIR, "ADSL.csv"), stringsAsFactors = FALSE) %>% select(USUBJID, DOSELEVEL)

pt <- adex %>%
  filter(ARM == "TREATMENT") %>%
  group_by(USUBJID) %>%
  summarise(rdi = 100 * mean(DOSEINT, na.rm = TRUE)) %>%
  left_join(adsl, by = "USUBJID") %>%
  filter(DOSELEVEL > 0) %>%
  mutate(dose_f = factor(paste0(DOSELEVEL, " mg"), levels = paste0(sort(unique(DOSELEVEL)), " mg")))

med_lab <- pt %>% group_by(dose_f) %>% summarise(n = n(), med = median(rdi)) %>%
  mutate(lab = sprintf("%s\n(n=%d)", dose_f, n))
pt <- pt %>% left_join(med_lab %>% select(dose_f, lab), by = "dose_f")
lab_levels <- med_lab$lab[order(med_lab$dose_f)]
pt$lab <- factor(pt$lab, levels = lab_levels)

p <- ggplot(pt, aes(x = lab, y = rdi)) +
  geom_boxplot(fill = "#9ECAE1", color = "black", width = 0.5, outlier.color = "#B2182B") +
  geom_jitter(width = 0.06, alpha = 0.5, color = "#2166AC", size = 1.6) +
  geom_hline(yintercept = 100, color = "forestgreen", linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = 80, color = "#B2182B", linetype = "dotted", linewidth = 0.8) +
  annotate("text", x = 0.6, y = 101, label = "Planned dose intensity (100%)", hjust = 0, size = 3, color = "forestgreen") +
  annotate("text", x = 0.6, y = 78.5, label = "Clinically relevant reduction threshold (80%)", hjust = 0, size = 3, color = "#B2182B") +
  labs(x = NULL, y = "Relative Dose Intensity, RDI (%)\n= (Actual dose/time) \u00f7 (Planned dose/time) \u00d7 100",
       title = "Dose Intensity by Planned Dose Level",
       caption = sprintf("Per-patient RDI = mean of per-cycle dose intensity. Median RDI \u2014 %s",
                          paste(sprintf("%s: %.1f%%", med_lab$dose_f, med_lab$med), collapse = " | "))) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        plot.caption = element_text(size = 8, face = "italic", hjust = 0.5))

ggsave(file.path(OUTPUT_DIR, "dose_intensity_plot.png"), p, width = 8, height = 6.5, dpi = 300)
