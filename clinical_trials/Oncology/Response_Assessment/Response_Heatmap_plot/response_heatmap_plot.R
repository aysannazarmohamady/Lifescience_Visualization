library(dplyr)
library(ggplot2)
library(tidyr)

DATA_DIR   <- "./Data/V1"
OUTPUT_DIR <- "./Outputs"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

adsl <- read.csv(file.path(DATA_DIR, "ADSL.csv"), stringsAsFactors = FALSE)

trt <- adsl |> filter(ARM == "TREATMENT") |>
  mutate(COHORT = case_when(
    DOSELEVEL == 100 ~ "Dose Level 1\n(100 mg)",
    DOSELEVEL == 200 ~ "Dose Level 2\n(200 mg)",
    DOSELEVEL == 300 ~ "RP2D\n(300 mg)",
    DOSELEVEL == 400 ~ "Dose Level 4\n(400 mg)"
  ))

COHORT_ORDER <- c("Dose Level 1\n(100 mg)","Dose Level 2\n(200 mg)",
                  "RP2D\n(300 mg)","Dose Level 4\n(400 mg)")
TUMOR_ORDER  <- c("NSCLC","BRCA","CRC","HCC","PDAC")
CUTOFF       <- "05 Mar 2026"

cell_df <- expand.grid(COHORT = COHORT_ORDER,
                       TUMORTYPE = TUMOR_ORDER,
                       stringsAsFactors = FALSE) |>
  left_join(
    trt |>
      group_by(COHORT, TUMORTYPE) |>
      summarise(
        n   = n(),
        n_cr = sum(BESTRSPC == "CR"),
        n_pr = sum(BESTRSPC == "PR"),
        n_sd = sum(BESTRSPC == "SD"),
        n_pd = sum(BESTRSPC == "PD"),
        orr  = round((sum(BESTRSPC %in% c("CR","PR")) / n()) * 100),
        .groups = "drop"
      ),
    by = c("COHORT","TUMORTYPE")
  ) |>
  mutate(
    n    = replace_na(n,    0L),
    n_cr = replace_na(n_cr, 0L),
    n_pr = replace_na(n_pr, 0L),
    n_sd = replace_na(n_sd, 0L),
    n_pd = replace_na(n_pd, 0L),
    orr  = ifelse(n == 0, NA_real_, orr),
    COHORT    = factor(COHORT,    levels = rev(COHORT_ORDER)),
    TUMORTYPE = factor(TUMORTYPE, levels = TUMOR_ORDER),
    counts_lbl = ifelse(n > 0,
      sprintf("CR:%d  PR:%d  SD:%d  PD:%d", n_cr, n_pr, n_sd, n_pd), NA),
    n_lbl = ifelse(n > 0, sprintf("n=%d", n), NA),
    orr_lbl   = ifelse(n > 0, sprintf("%d%%", orr), "n/a"),
    txt_color = ifelse(
      !is.na(orr) & (orr >= 65 | orr == 0),
      "white", "#222222"
    )
  )

pal <- c("#2166AC","#92C5DE","#F7F7F7","#F4A582","#B2182B")
breaks_pal <- c(0, 0.30, 0.50, 0.70, 1.00)

p <- ggplot(cell_df, aes(x = TUMORTYPE, y = COHORT)) +
  geom_tile(aes(fill = orr / 100),
            color = "white", linewidth = 1.2,
            width = 0.92, height = 0.92) +
  geom_tile(data = filter(cell_df, n == 0),
            aes(x = TUMORTYPE, y = COHORT),
            fill = "#EEEEEE", color = "white",
            linewidth = 1.2, width = 0.92, height = 0.92) +
  geom_text(aes(label = orr_lbl, color = txt_color),
            fontface = "bold", size = 5.8,
            nudge_y = 0.12) +
  geom_text(aes(label = counts_lbl, color = txt_color),
            size = 2.8, nudge_y = -0.10) +
  geom_text(aes(label = n_lbl, color = txt_color),
            size = 3.0, nudge_y = -0.28) +
  geom_text(data = filter(cell_df, n == 0),
            aes(label = "n/a"), color = "#AAAAAA",
            fontface = "bold.italic", size = 4.5) +
  scale_fill_gradientn(
    colors   = pal,
    values   = breaks_pal,
    limits   = c(0, 1),
    na.value = "#EEEEEE",
    labels   = c("0","25","50","75","100"),
    breaks   = c(0, 0.25, 0.50, 0.75, 1.00),
    name     = "ORR (%)"
  ) +
  scale_color_identity() +
  scale_x_discrete(position = "bottom", expand = c(0.02, 0.02)) +
  scale_y_discrete(expand = c(0.06, 0.06)) +
  guides(fill = guide_colorbar(
    barheight = unit(8, "cm"),
    barwidth  = unit(0.55, "cm"),
    ticks.linewidth = 0.5,
    frame.colour = "grey60",
    title.position = "top",
    title.hjust    = 0.5
  )) +
  labs(
    title    = "Objective Response Rate (ORR) by Dose Cohort and Tumor Type",
    subtitle = "Color intensity reflects ORR  ·  CR/PR/SD/PD counts per cell",
    x        = "Tumor Type",
    y        = "Dose Cohort",
    caption  = paste0("RECIST 1.1  ·  Treatment arm only  ·  Data cutoff: ", CUTOFF)
  ) +
  theme_minimal(base_size = 11) +
  theme(
    plot.title       = element_text(face = "bold", size = 14,
                                    hjust = 0.5, margin = margin(b = 4)),
    plot.subtitle    = element_text(size = 9, hjust = 0.5, color = "#555555",
                                    face = "italic", margin = margin(b = 10)),
    plot.caption     = element_text(size = 8.5, hjust = 0.5, color = "#666666",
                                    face = "italic", margin = margin(t = 8)),
    axis.title.x     = element_text(face = "bold", size = 12,
                                    margin = margin(t = 8)),
    axis.title.y     = element_text(face = "bold", size = 12,
                                    margin = margin(r = 8)),
    axis.text.x      = element_text(face = "bold", size = 11),
    axis.text.y      = element_text(face = "bold", size = 10,
                                    lineheight = 1.2),
    panel.grid       = element_blank(),
    legend.title     = element_text(face = "bold", size = 9),
    legend.text      = element_text(size = 9),
    legend.position  = "right",
    plot.margin      = margin(14, 14, 10, 14),
    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

ggsave(file.path(OUTPUT_DIR, "response_heatmap.png"),
       plot   = p,
       width  = 14,
       height = 9,
       dpi    = 150,
       bg     = "white")

message("Saved: response_heatmap.png")
