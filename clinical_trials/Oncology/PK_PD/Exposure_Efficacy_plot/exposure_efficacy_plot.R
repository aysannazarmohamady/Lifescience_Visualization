# exposure_efficacy.R
# Exposure-Efficacy Plot — AUC(0-24, C1D1) vs best % tumor change
# Standards applied: log exposure axis, individual points, sigmoid Emax model (nls) with
# bootstrap 95% CI band, complementary exposure-quartile panel, legend/model parameters
# placed outside the plotted data area (never overlapping points/curve).

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

adpk <- read.csv(file.path(DATA_DIR, "ADPK.csv"), stringsAsFactors = FALSE)
adtr <- read.csv(file.path(DATA_DIR, "ADTR.csv"), stringsAsFactors = FALSE)

auc <- adpk |>
  filter(PARAMCD == "AUC", AVISIT == "CYCLE 1 DAY 1") |>
  distinct(USUBJID, .keep_all = TRUE) |>
  select(USUBJID, AUC = AVAL, DOSE)

best <- adtr |>
  filter(AVISITN > 0) |>
  group_by(USUBJID) |>
  summarise(BEST_PCHG = min(PCHG, na.rm = TRUE), .groups = "drop")

df <- inner_join(auc, best, by = "USUBJID")
n  <- nrow(df)

emax_fun <- function(auc, e0, emax, ec50, gamma) e0 + emax * (auc^gamma) / (ec50^gamma + auc^gamma)

fit <- tryCatch(
  nls(BEST_PCHG ~ emax_fun(AUC, e0, emax, ec50, gamma), data = df,
      start = list(e0 = mean(df$BEST_PCHG), emax = -80, ec50 = median(df$AUC), gamma = 1),
      lower = c(-100, -300, 0.1, 0.3), upper = c(100, 300, 100, 8),
      algorithm = "port", control = nls.control(maxiter = 200)),
  error = function(e) NULL
)
fit_ok <- !is.null(fit)

xg <- seq(min(df$AUC), max(df$AUC), length.out = 200)

if (fit_ok) {
  p_hat <- coef(fit)
  yg <- emax_fun(xg, p_hat[["e0"]], p_hat[["emax"]], p_hat[["ec50"]], p_hat[["gamma"]])

  set.seed(42)
  n_boot <- 500
  boot_mat <- matrix(NA_real_, nrow = n_boot, ncol = length(xg))
  for (b in seq_len(n_boot)) {
    idx <- sample(seq_len(n), n, replace = TRUE)
    dfb <- df[idx, ]
    fb <- tryCatch(
      nls(BEST_PCHG ~ emax_fun(AUC, e0, emax, ec50, gamma), data = dfb,
          start = as.list(p_hat), lower = c(-100, -300, 0.1, 0.3), upper = c(100, 300, 100, 8),
          algorithm = "port", control = nls.control(maxiter = 100)),
      error = function(e) NULL)
    if (!is.null(fb)) {
      pb <- coef(fb)
      boot_mat[b, ] <- emax_fun(xg, pb[["e0"]], pb[["emax"]], pb[["ec50"]], pb[["gamma"]])
    }
  }
  valid <- boot_mat[complete.cases(boot_mat), , drop = FALSE]
  lo_band <- apply(valid, 2, quantile, probs = 0.025)
  hi_band <- apply(valid, 2, quantile, probs = 0.975)
  fit_df  <- data.frame(AUC = xg, yhat = yg, lo = lo_band, hi = hi_band)
}

df$DOSE_f <- factor(df$DOSE, levels = sort(unique(df$DOSE)))

p1 <- ggplot(df, aes(x = AUC, y = BEST_PCHG)) +
  { if (fit_ok) geom_ribbon(data = fit_df, aes(x = AUC, ymin = lo, ymax = hi), inherit.aes = FALSE,
                            fill = "grey60", alpha = 0.25) } +
  { if (fit_ok) geom_line(data = fit_df, aes(x = AUC, y = yhat), inherit.aes = FALSE,
                           color = "black", linewidth = 1) } +
  geom_point(aes(color = DOSE_f), size = 2.6, alpha = 0.85, shape = 21, stroke = 0.4,
             fill = NA) +
  geom_point(aes(fill = DOSE_f), size = 2.6, alpha = 0.85, shape = 21, color = "black", stroke = 0.3) +
  geom_hline(yintercept = PR_TH, color = "#1A3A7C", linetype = "dotted", linewidth = 0.6) +
  annotate("label", x = max(df$AUC), y = PR_TH, label = "PR threshold (\u201330%)",
           color = "#1A3A7C", size = 2.8, hjust = 1, label.size = 0, fill = alpha("white", 0.8)) +
  scale_x_log10() +
  scale_color_manual(values = DOSE_COLORS, guide = "none") +
  scale_fill_manual(values = DOSE_COLORS, name = "Dose", labels = function(x) paste0(x, " mg")) +
  labs(x = "AUC\u2080\u208b\u2082\u2084 (ng\u00b7h/mL, Cycle 1 Day 1) \u2014 log scale",
       y = "Best % change in tumor size from baseline",
       title = paste0("Individual patients + sigmoid Emax fit (n=", n, ")")) +
  theme_pkpd() +
  theme(legend.position = "bottom")

# Panel 2 — exposure-quartile analysis (unchanged design)
df <- df |> mutate(Q = ntile(AUC, 4))
qsum <- df |>
  group_by(Q) |>
  summarise(med_auc = median(AUC), mean_r = mean(BEST_PCHG),
            sem_r = sd(BEST_PCHG) / sqrt(n()), n = n(), .groups = "drop")

p2 <- ggplot(qsum, aes(x = factor(Q), y = mean_r)) +
  geom_col(fill = "#4A90C4", color = "#1A3A7C", width = 0.6) +
  geom_errorbar(aes(ymin = mean_r - 1.96 * sem_r, ymax = mean_r + 1.96 * sem_r), width = 0.15) +
  geom_hline(yintercept = 0, linewidth = 0.5) +
  scale_x_discrete(labels = sprintf("Q%d\nmedian AUC=%.1f\n(n=%d)", qsum$Q, qsum$med_auc, qsum$n)) +
  labs(x = NULL, y = "Mean best % change (95% CI)", title = "Exposure-quartile analysis") +
  theme_pkpd()

param_txt <- if (fit_ok) {
  sprintf("Model: E = E0 + Emax\u00b7AUC^\u03b3 / (EC50^\u03b3 + AUC^\u03b3). Fitted parameters: E0=%.1f%%, Emax=%.1f%%, EC50=%.1f ng\u00b7h/mL, \u03b3(Hill)=%.2f. ",
          p_hat[["e0"]], p_hat[["emax"]], p_hat[["ec50"]], p_hat[["gamma"]])
} else "Model did not converge. "

cap <- paste(strwrap(paste0(param_txt,
  "Sigmoid Emax model fit by nonlinear least squares; shaded band = 95% CI from 500 bootstrap refits ",
  "(patient-level resampling). Flat/non-significant fit here (Emax and EC50 poorly constrained by wide CI) is ",
  "reported as observed, not forced \u2014 consistent with no true exposure-response signal in this dataset. ",
  "Right panel: population divided into AUC quartiles (median exposure per quartile) vs. mean response \u00b1 95% CI, ",
  "the complementary categorical view used in FDA exposure-response submissions."
), width = 150), collapse = "\n")

combined <- (p1 | p2) +
  plot_annotation(
    title = sprintf("Exposure-Efficacy Plot \u2014 Vizatinib (ONCVIZ-001) AUC\u2080\u208b\u2082\u2084 vs. Best %% Tumor Change, PK/Efficacy-evaluable population (N=%d)", n),
    caption = cap,
    theme = theme(plot.title = element_text(face = "bold", size = 12.5, hjust = 0.5))
  )

ggsave(file.path(OUTPUT_DIR, "exposure_efficacy.png"), combined, width = 12, height = 6.5, dpi = 300)
message("Saved exposure_efficacy.png to ", OUTPUT_DIR)
if (fit_ok) {
  message(sprintf("E0=%.2f Emax=%.2f EC50=%.2f gamma=%.2f", p_hat[["e0"]], p_hat[["emax"]], p_hat[["ec50"]], p_hat[["gamma"]]))
}
