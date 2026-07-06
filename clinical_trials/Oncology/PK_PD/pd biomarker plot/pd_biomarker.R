# pd_biomarker.R
# PD Biomarker Plot — ctDNA % change from baseline over time, by best overall response (BOR)
# Standards applied: baseline-normalized (%CHG), mean +/- SEM per group, secondary top axis
# (nominal study day = cycle x 21), consistent RESP_COLORS mapping, N per group in legend.

suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
})

SCRIPT_DIR <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)))
if (length(SCRIPT_DIR) == 0 || SCRIPT_DIR == "") SCRIPT_DIR <- getwd()
source(file.path(SCRIPT_DIR, "..", "pkpd_common.R"))

DATA_DIR   <- file.path(SCRIPT_DIR, "..", "Data", "V1")
OUTPUT_DIR <- file.path(SCRIPT_DIR, "Out")
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

adbm <- read.csv(file.path(DATA_DIR, "ADBM.csv"), stringsAsFactors = FALSE)

CYCLE_DAYS <- 21  # confirmed exact: ADTN = AVISITN * 21

ct <- adbm |> filter(PARAMCD == "CTDNA", ARM == "TREATMENT", AVISITN >= 0)

resp_order <- c("CR", "PR", "SD", "PD")

summ <- ct |>
  filter(BESTRSPC %in% resp_order) |>
  group_by(BESTRSPC, AVISITN) |>
  summarise(mean_pchg = mean(PCHG, na.rm = TRUE),
            sem = sd(PCHG, na.rm = TRUE) / sqrt(sum(!is.na(PCHG))),
            n = sum(!is.na(PCHG)), .groups = "drop") |>
  mutate(BESTRSPC = factor(BESTRSPC, levels = resp_order))

baseline_n <- summ |> filter(AVISITN == 0) |> select(BESTRSPC, n0 = n)
legend_labels <- setNames(
  paste0(resp_order, " (baseline n=", baseline_n$n0[match(resp_order, baseline_n$BESTRSPC)], ")"),
  resp_order
)

visits <- sort(unique(summ$AVISITN))

p <- ggplot(summ, aes(x = AVISITN, y = mean_pchg, color = BESTRSPC, fill = BESTRSPC)) +
  geom_ribbon(aes(ymin = mean_pchg - sem, ymax = mean_pchg + sem), alpha = 0.18, color = NA) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  geom_hline(yintercept = 0, linetype = "dashed", linewidth = 0.5, color = "black", alpha = 0.6) +
  scale_color_manual(values = RESP_COLORS, labels = legend_labels, name = "Best Overall Response") +
  scale_fill_manual(values = RESP_COLORS, labels = legend_labels, name = "Best Overall Response") +
  scale_x_continuous(breaks = visits,
                      sec.axis = sec_axis(~ . * CYCLE_DAYS, name = "Nominal study day post first dose (cycle \u00d7 21)",
                                          breaks = visits * CYCLE_DAYS)) +
  labs(x = "Study cycle (AVISITN)", y = "ctDNA, % change from baseline (mean \u00b1 SEM)",
       title = "PD Biomarker Plot \u2014 ctDNA Dynamics by Best Overall Response",
       subtitle = "Vizatinib (ONCVIZ-001), Treatment arm, PD-evaluable population") +
  theme_pkpd() +
  theme(legend.position = c(0.22, 0.82), legend.background = element_rect(fill = alpha("white", 0.7), color = NA))

cap <- paste(strwrap(paste0(
  "Baseline-normalized (% change from baseline); mean \u00b1 SEM shown, not SD, per group PD summary convention. ",
  "N per group = baseline N shown in legend; N decreases at later cycles due to progression-driven discontinuation ",
  "(attrition, not shown per curve to avoid clutter \u2014 see source data for per-timepoint N). ",
  "Note on time axis: PD sampling here is longitudinal across cycles (days-to-weeks scale, top axis), whereas the ",
  "PK concentration-time plot is intensive single-dose sampling within Cycle 1 Day 1 (hours scale). The two plots ",
  "are intentionally on different time bases and are not meant to be read off the same axis; the top axis above ",
  "converts cycle to nominal day for cross-study orientation only."
), width = 150), collapse = "\n")

p <- p + labs(caption = cap)

ggsave(file.path(OUTPUT_DIR, "pd_biomarker.png"), p, width = 9, height = 7.2, dpi = 300)
message("Saved pd_biomarker.png to ", OUTPUT_DIR)
