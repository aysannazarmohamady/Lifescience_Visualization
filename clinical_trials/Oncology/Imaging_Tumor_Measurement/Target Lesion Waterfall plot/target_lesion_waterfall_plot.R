suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

DATA_DIR   <- "./Data/V1"
OUTPUT_DIR <- "./Out"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

adtr <- read.csv(file.path(DATA_DIR, "ADTR.csv"), stringsAsFactors = FALSE) %>% filter(PARAMCD == "SUMDIAM")
adsl <- read.csv(file.path(DATA_DIR, "ADSL.csv"), stringsAsFactors = FALSE) %>% select(USUBJID, BESTRSPC)

best <- adtr %>% filter(AVISITN > 0) %>%
  group_by(USUBJID) %>% slice_min(PCHG, n = 1, with_ties = FALSE) %>% ungroup() %>%
  left_join(adsl, by = "USUBJID") %>%
  arrange(desc(PCHG)) %>%
  mutate(x = row_number())

color_map <- c(CR = "#2166AC", PR = "#4393C3", SD = "#FDB863", PD = "#B2182B")
n_counts <- best %>% count(BESTRSPC)

p <- ggplot(best, aes(x = x, y = PCHG, fill = BESTRSPC)) +
  geom_col(width = 0.85, color = "black", linewidth = 0.15) +
  geom_hline(yintercept = 20, color = "#B2182B", linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = -30, color = "#2166AC", linetype = "dashed", linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.4) +
  annotate("text", x = max(best$x) + 1, y = 20, label = "+20% (PD)", color = "#B2182B", hjust = 0, vjust = -0.3, size = 3) +
  annotate("text", x = max(best$x) + 1, y = -30, label = "\u221230% (PR)", color = "#2166AC", hjust = 0, vjust = 1.2, size = 3) +
  scale_fill_manual(values = color_map,
                     labels = sprintf("%s (n=%d)", n_counts$BESTRSPC, n_counts$n), name = NULL) +
  coord_cartesian(xlim = c(-1, max(best$x) + 9), clip = "off") +
  labs(x = sprintf("Individual Patients (n = %d, sorted by best %% change)", nrow(best)),
       y = "Best % Change in Target Lesion SLD from Baseline",
       title = "Target Lesion Change (Waterfall Plot)\nBest Overall Response per RECIST 1.1",
       caption = "Best % change = minimum post-baseline SLD change (nadir) per patient. RECIST 1.1; evaluable population (patients with \u22651 post-baseline target-lesion assessment).") +
  theme_classic(base_size = 12) +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        plot.caption = element_text(size = 8, face = "italic", hjust = 0.5),
        plot.margin = margin(5.5, 40, 5.5, 5.5))

ggsave(file.path(OUTPUT_DIR, "target_lesion_waterfall.png"), p, width = 12, height = 6.5, dpi = 300)
