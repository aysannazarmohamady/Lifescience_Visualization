# waterfall_pk_overlay.R
# Waterfall Plot + PK Overlay — best % tumor change per patient (RECIST-colored bars),
# sorted descending, with a per-patient AUC heatmap strip overlay beneath the bars.
# Standards applied: RECIST reference lines (-30%/+20%) placed in a verified-empty interior
# region, legend placed above the plot (outside the bar area), heatmap-strip PK overlay.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(patchwork)
})

SCRIPT_DIR <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
if (length(SCRIPT_DIR) == 0 || SCRIPT_DIR == "") SCRIPT_DIR <- getwd()
source(file.path(SCRIPT_DIR, "..", "pkpd_common.R"))

DATA_DIR   <- file.path(SCRIPT_DIR, "..", "Data", "V1")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "Out")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

adsl <- read.csv(file.path(DATA_DIR, "ADSL.csv"), stringsAsFactors = FALSE)
adtr <- read.csv(file.path(DATA_DIR, "ADTR.csv"), stringsAsFactors = FALSE)
adpk <- read.csv(file.path(DATA_DIR, "ADPK.csv"), stringsAsFactors = FALSE)

trt  <- adsl |> filter(ARM == "TREATMENT") |> select(USUBJID, BESTRSPC)
best <- adtr |> filter(AVISITN > 0) |> group_by(USUBJID) |>
  summarise(BEST_PCHG = min(PCHG, na.rm = TRUE), .groups = "drop")
auc  <- adpk |> filter(PARAMCD == "AUC", AVISIT == "CYCLE 1 DAY 1") |>
  distinct(USUBJID, .keep_all = TRUE) |> select(USUBJID, AUC = AVAL)

wf <- trt |> left_join(best, by = "USUBJID") |> left_join(auc, by = "USUBJID")

n_no_tumor <- sum(is.na(wf$BEST_PCHG))
wf_plot <- wf |> filter(!is.na(BEST_PCHG)) |> arrange(desc(BEST_PCHG)) |>
  mutate(x = row_number())
n_missing_pk <- sum(is.na(wf_plot$AUC))
n_tot <- nrow(wf_plot)

counts <- table(wf_plot$BESTRSPC)
orr <- 100 * (sum(counts[c("CR", "PR")], na.rm = TRUE)) / n_tot

# Verified-empty interior region for threshold labels: x=26, where bars sit near 0%
p_bar <- ggplot(wf_plot, aes(x = x, y = BEST_PCHG, fill = BESTRSPC)) +
  geom_col(width = 0.82) +
  geom_hline(yintercept = c(PR_TH, PD_TH), linetype = "dashed", linewidth = 0.5, color = "black", alpha = 0.7) +
  annotate("label", x = 26, y = PD_TH, label = "PD threshold (+20%)", size = 2.8,
           vjust = -0.3, label.size = 0, fill = alpha("white", 0.85)) +
  annotate("label", x = 26, y = PR_TH, label = "PR threshold (\u201330%)", size = 2.8,
           vjust = 1.3, label.size = 0, fill = alpha("white", 0.85)) +
  scale_fill_manual(values = RESP_COLORS, name = "Best Overall Response",
                     labels = paste0(names(counts), " (n=", as.integer(counts), ")"),
                     breaks = names(counts)) +
  scale_x_continuous(limits = c(0, n_tot + 1), expand = c(0, 0)) +
  labs(y = "Best % change in tumor size from baseline", x = NULL,
       title = sprintf("Waterfall Plot + PK Exposure Overlay \u2014 Vizatinib (ONCVIZ-001), Treatment arm (N=%d, ORR=%.0f%%)", n_tot, orr)) +
  theme_pkpd() +
  theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(),
        legend.position = "top", plot.title = element_text(size = 13, hjust = 0.5))

# PK overlay strip — log-scale AUC heatmap, same patient order
p_pk <- ggplot(wf_plot, aes(x = x, y = 1, fill = AUC)) +
  geom_tile() +
  scale_fill_gradientn(colors = c("#FFFFCC", "#FD8D3C", "#800026"), trans = "log10",
                       na.value = "#DDDDDD", name = "AUC\u2080\u208b\u2082\u2084\n(ng\u00b7h/mL,\nlog scale)") +
  scale_x_continuous(limits = c(0, n_tot + 1), expand = c(0, 0)) +
  scale_y_continuous(breaks = 1, labels = "AUC\u2080\u208b\u2082\u2084") +
  labs(x = "Patients, sorted by best % change (rank order; not a measurement scale)", y = NULL) +
  theme_void(base_size = 10) +
  theme(axis.text.y = element_text(size = 8.5), axis.title.x = element_text(size = 9.5, margin = margin(t = 6)),
        legend.position = "right")

cap <- paste(strwrap(sprintf(paste0(
  "Bars colored by RECIST 1.1 best overall response. %d treatment-arm patient(s) with no post-baseline tumor ",
  "assessment excluded from the waterfall (not shown). %d patient(s) shown lack a Cycle 1 AUC value (gray in PK strip)."
), n_no_tumor, n_missing_pk), width = 150), collapse = "\n")

combined <- (p_bar / p_pk) +
  plot_layout(heights = c(4, 0.6)) +
  plot_annotation(caption = cap,
                   theme = theme(plot.caption = element_text(size = 8, hjust = 0.5, face = "italic")))

ggsave(file.path(OUTPUT_DIR, "waterfall_pk_overlay.png"), combined, width = 13, height = 7.5, dpi = 300)
message("Saved waterfall_pk_overlay.png to ", OUTPUT_DIR)
