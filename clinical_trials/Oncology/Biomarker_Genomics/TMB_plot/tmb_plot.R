suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR   <- "./Data/V1"
OUTPUT_DIR <- "./clinical_trials/Oncology/Biomarker_Genomics/TMB_plot/output"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Load data ─────────────────────────────────────────────────────────────────
admut <- read.csv(file.path(DATA_DIR, "ADMUT.csv"), stringsAsFactors = FALSE)

tmb <- admut %>% distinct(USUBJID, TUMORTYPE, ARM, TMB, TMBHIGH, BESTRSPC)

# Order tumor types by median TMB (descending) for readability
tt_order <- tmb %>% group_by(TUMORTYPE) %>% summarise(med = median(TMB)) %>%
  arrange(desc(med)) %>% pull(TUMORTYPE)
tmb$TUMORTYPE <- factor(tmb$TUMORTYPE, levels = tt_order)

color_map <- c(CR = "#2166AC", PR = "#4393C3", SD = "#FDB863", PD = "#B2182B")
tmb_high_thresh <- 10  # mut/Mb, standard clinical TMB-High cutoff

n_high <- sum(tmb$TMBHIGH == "Y")
n_tot  <- nrow(tmb)

caption_txt <- sprintf("n = %d patients with mutation calls. TMB-High (\u2265%g mut/Mb): %d/%d patients (%.0f%%). One point per patient; jittered horizontally within tumor type.",
                        n_tot, tmb_high_thresh, n_high, n_tot, 100 * n_high / n_tot)
caption_wrapped <- paste(strwrap(caption_txt, width = 100), collapse = "\n")

p <- ggplot(tmb, aes(x = TUMORTYPE, y = TMB)) +
  geom_hline(aes(yintercept = tmb_high_thresh, linetype = "TMB-High threshold (10 mut/Mb)"), color = "grey40", linewidth = 0.9) +
  geom_jitter(aes(fill = BESTRSPC), shape = 21, color = "black", stroke = 0.4,
              width = 0.15, size = 3.2, alpha = 0.85) +
  scale_fill_manual(values = color_map, name = NULL) +
  scale_linetype_manual(values = c("TMB-High threshold (10 mut/Mb)" = "dashed"), name = NULL) +
  labs(x = NULL, y = "Tumor Mutational Burden (mut/Mb)",
       title = "Tumor Mutational Burden by Tumor Type and Best Response",
       caption = caption_wrapped) +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13, margin = margin(b = 10)),
        plot.title.position = "plot",
        plot.caption = element_text(size = 8, face = "italic", hjust = 0.5, margin = margin(t = 10)),
        plot.caption.position = "plot",
        plot.margin = margin(t = 15, r = 15, b = 12, l = 12),
        legend.position = "right",
        legend.key = element_blank(),
        legend.background = element_blank(),
        legend.box.background = element_rect(color = "black", linewidth = 0.5, fill = "white"),
        legend.box.margin = margin(t = 4, r = 4, b = 4, l = 4),
        legend.margin = margin(t = 8, r = 10, b = 8, l = 10))

ggsave(file.path(OUTPUT_DIR, "tmb_plot.png"), p, width = 9, height = 6, dpi = 300)
