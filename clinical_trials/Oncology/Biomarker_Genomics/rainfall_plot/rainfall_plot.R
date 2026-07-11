library(ggplot2); library(dplyr)

set.seed(10)

OUTPUT_DIR <- "./clinical_trials/Oncology/Biomarker_Genomics/rainfall_plot/output"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

GENES <- c("ARID1A","PIK3CA","EGFR","MET","BRAF","CDKN2A","RET","PTEN",
           "KRAS","RB1","TP53","NF1","SMAD4","KEAP1","STK11")
CHR_MAP   <- c(ARID1A="chr1",PIK3CA="chr3",EGFR="chr7",MET="chr7",BRAF="chr7",
               CDKN2A="chr9",RET="chr10",PTEN="chr10",KRAS="chr12",RB1="chr13",
               TP53="chr17",NF1="chr17",SMAD4="chr18",KEAP1="chr19",STK11="chr19")
CHR_ORDER <- c("chr1","chr3","chr7","chr9","chr10","chr12","chr13","chr17","chr18","chr19")

VARIANT_CLASSES <- c("Missense Mutation","Nonsense Mutation","Frame Shift (Del)",
                      "Frame Shift (Ins)","Splice Site","In-Frame Del","Other Mutation")
CLASS_COLORS <- c(
  "Missense Mutation"  = "#2c5f8a",
  "Nonsense Mutation"  = "#c0392b",
  "Frame Shift (Del)"  = "#e67e22",
  "Frame Shift (Ins)"  = "#f1c40f",
  "Splice Site"        = "#8e44ad",
  "In-Frame Del"       = "#27ae60",
  "Other Mutation"     = "#7f8c8d"
)

admut <- data.frame(
  PATIENT_ID    = sample(sprintf("PT-%03d", 1:66), 138, replace = TRUE),
  HUGO_SYMBOL   = sample(GENES, 138, replace = TRUE),
  VARIANT_CLASS = sample(VARIANT_CLASSES, 138, replace = TRUE,
                          prob = c(0.42, 0.14, 0.12, 0.08, 0.09, 0.08, 0.07))
) %>%
  mutate(
    CHR         = CHR_MAP[HUGO_SYMBOL],
    PROTEIN_POS = round(runif(n(), 10, 900))
  )

rainfall_df <- admut %>%
  mutate(CHR = factor(CHR, levels = CHR_ORDER)) %>%
  arrange(CHR, PROTEIN_POS) %>%
  group_by(CHR) %>%
  mutate(
    DIST     = abs(PROTEIN_POS - lag(PROTEIN_POS)),
    DIST     = pmax(DIST, 1),
    LOG_DIST = log10(DIST),
    xi       = row_number()
  ) %>%
  ungroup() %>%
  filter(!is.na(LOG_DIST))

rainfall_df <- rainfall_df %>%
  arrange(CHR, PROTEIN_POS) %>%
  mutate(x_global = row_number())

chr_mid <- rainfall_df %>%
  group_by(CHR) %>%
  summarise(mid = mean(x_global), .groups = "drop")

CLUSTER_THRESHOLD_AA <- 300
threshold_log <- log10(CLUSTER_THRESHOLD_AA)

p <- ggplot(rainfall_df, aes(x = x_global, y = LOG_DIST, color = VARIANT_CLASS)) +
  geom_hline(yintercept = threshold_log, color = "#c0392b",
             linetype = "dashed", linewidth = 0.6) +
  geom_point(size = 2.0, alpha = 0.85) +
  scale_color_manual(values = CLASS_COLORS, name = "Variant classification") +
  scale_x_continuous(breaks = chr_mid$mid, labels = chr_mid$CHR, name = "Chromosome") +
  scale_y_continuous(name = expression(log[10]*"(inter-mutation distance, AA)")) +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x      = element_text(face = "bold"),
    panel.background = element_rect(fill = "white", color = NA)
  ) +
  labs(
    title    = "Rainfall Plot — Inter-Mutation Distance Across the Sequenced Gene Panel",
    subtitle = paste0("n = ", nrow(rainfall_df), " mutations (of 138 total) · 15 genes · ",
                       "cluster threshold = ", CLUSTER_THRESHOLD_AA, " AA (red dashed line)"),
    caption  = "Distance is a protein-coordinate proxy (amino-acid position), not chromosomal base-pair distance."
  )

ggsave(file.path(OUTPUT_DIR, "rainfall_plot.png"), p, width = 15, height = 6, dpi = 180, bg = "white")
message("Saved: rainfall_plot.png")
