suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

DATA_DIR   <- "./Data/V1"
OUTPUT_DIR <- "./Out"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

adtr <- read.csv(file.path(DATA_DIR, "ADTR.csv"), stringsAsFactors = FALSE) %>% filter(PARAMCD == "SUMDIAM")
adsl <- read.csv(file.path(DATA_DIR, "ADSL.csv"), stringsAsFactors = FALSE) %>% select(USUBJID, BESTRSPC)

tr <- adtr %>% mutate(weeks = ADTN / 7) %>% left_join(adsl, by = "USUBJID")

color_map <- c(CR = "#2166AC", PR = "#4393C3", SD = "#FDB863", PD = "#B2182B")

p <- ggplot(tr, aes(x = weeks, y = PCHG, group = USUBJID, color = BESTRSPC)) +
  geom_line(alpha = 0.6, linewidth = 0.5) +
  geom_point(size = 0.9, alpha = 0.6) +
  geom_hline(yintercept = 20, color = "#B2182B", linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = -30, color = "#2166AC", linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  scale_color_manual(values = color_map, name = "Best response") +
  labs(x = "Time Since Baseline (weeks)", y = "% Change in Sum of Longest Diameters (SLD) from Baseline",
       title = "Sum of Longest Diameters Over Time (Spaghetti Plot)\nper RECIST 1.1",
       caption = sprintf("n = %d evaluable patients with target-lesion assessments; markers = actual scan timepoints (unsmoothed).", n_distinct(tr$USUBJID))) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        plot.caption = element_text(size = 8, face = "italic", hjust = 0.5),
        legend.position = "top")

ggsave(file.path(OUTPUT_DIR, "sld_over_time.png"), p, width = 10, height = 7.3, dpi = 300)
