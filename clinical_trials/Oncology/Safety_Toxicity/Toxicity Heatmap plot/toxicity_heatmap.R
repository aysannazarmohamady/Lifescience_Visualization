suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR   <- "./Data/V1"
OUTPUT_DIR <- "./Out"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Load data ─────────────────────────────────────────────────────────────────
adae <- read.csv(file.path(DATA_DIR, "ADAE.csv"), stringsAsFactors = FALSE)
adsl <- read.csv(file.path(DATA_DIR, "ADSL.csv"), stringsAsFactors = FALSE)

n_trt <- sum(adsl$ARM == "TREATMENT")

# ── Worst grade per patient per System Organ Class (Treatment arm) ────────────
worst <- adae %>%
  filter(ARM == "TREATMENT") %>%
  group_by(USUBJID, AESOC) %>%
  summarise(worst_grade = max(AETOXGR), .groups = "drop")

tab <- worst %>%
  count(AESOC, worst_grade, name = "n_pt") %>%
  complete(AESOC, worst_grade = 1:5, fill = list(n_pt = 0)) %>%
  mutate(pct = 100 * n_pt / n_trt)

soc_order <- tab %>% group_by(AESOC) %>% summarise(tot = sum(pct)) %>% arrange(tot) %>% pull(AESOC)
tab$AESOC <- factor(tab$AESOC, levels = soc_order)
tab$grade_lab <- factor(paste0("G", tab$worst_grade), levels = paste0("G", 1:5))

ctcae_v <- unique(na.omit(adae$CTCAE_V))[1]

caption_txt <- sprintf("Treatment arm, N = %d. Cell = %% of patients whose worst grade for that SOC equals the column grade (CTCAE v%s). Rows sorted by total burden.",
                        n_trt, ctcae_v)
caption_wrapped <- paste(strwrap(caption_txt, width = 90), collapse = "\n")

p <- ggplot(tab, aes(x = grade_lab, y = AESOC, fill = pct)) +
  geom_tile(color = NA) +
  geom_text(aes(label = ifelse(pct > 0, sprintf("%.0f", pct), "")),
            color = ifelse(tab$pct > 15, "white", "black"), size = 3.2) +
  scale_fill_gradient(low = "#FFF5F0", high = "#A50F15", name = "% patients\n(worst grade)") +
  scale_x_discrete(expand = c(0, 0)) +
  scale_y_discrete(expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       title = "Toxicity Heatmap \u2014 Worst Grade by\nSystem Organ Class (Treatment Arm)",
       caption = caption_wrapped) +
  theme_classic(base_size = 12) +
  theme(axis.line = element_blank(),
        axis.ticks = element_blank(),
        panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
        plot.title = element_text(hjust = 0.5, face = "bold", size = 13, margin = margin(b = 10)),
        plot.title.position = "plot",
        plot.caption = element_text(size = 8, face = "italic", hjust = 0.5, margin = margin(t = 10)),
        plot.caption.position = "plot",
        plot.margin = margin(t = 15, r = 20, b = 12, l = 12),
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(color = "black", linewidth = 0.5, fill = "white"),
        legend.box.margin = margin(t = 4, r = 4, b = 4, l = 4),
        legend.margin = margin(t = 8, r = 10, b = 8, l = 10))

ggsave(file.path(OUTPUT_DIR, "toxicity_heatmap.png"), p, width = 8, height = 6.3, dpi = 300)
