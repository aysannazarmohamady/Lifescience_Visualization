# =============================================================================
# ONCVIZ-001 -- Kaplan-Meier plot, variant 2 (benchmark second implementation)
#
# How this differs from survivalcurves_plot.R (variant 1):
#   * NO survminer: the product-limit estimate is computed here, step data is
#     drawn with ggplot2 and the number-at-risk table with patchwork. The clean
#     curve is cross-checked against survival::survfit when the script runs;
#   * time axis in months (variant 1: days); endpoint OS, arms in one panel,
#     at-risk table every 6 months;
#   * errors are applied to the plotted step data / table by switches in the
#     runner, so no edited copies of the script are needed.
#
# DRAFT NOTE: drafted with LLM assistance; an independent author must rewrite
# or review it if the paper claims hand-written benchmark scripts.
# Usage (repository root): Rscript <this file>   (draws the clean figure)
# =============================================================================

suppressPackageStartupMessages({ library(dplyr); library(ggplot2); library(patchwork); library(survival) })

KM_ARMS   <- c(CONTROL = "#B2182B", TREATMENT = "#2166AC")
KM_MONTH  <- 30.4375
KM_TICKS  <- seq(0, 48, by = 6)

# Product-limit estimate for one arm. Returns step data (with a t = 0 row) and
# the censoring marks (time, survival at that time).
km_estimate <- function(d) {
  d <- d |> mutate(t = AVAL / KM_MONTH, ev = 1L - CNSR) |> arrange(t)
  tt <- sort(unique(d$t))
  steps <- bind_rows(lapply(tt, function(x) {
    data.frame(t = x, n_risk = sum(d$t >= x), n_event = sum(d$t == x & d$ev == 1),
               n_cens = sum(d$t == x & d$ev == 0))
  }))
  steps$surv <- cumprod(1 - steps$n_event / steps$n_risk)
  steps <- bind_rows(data.frame(t = 0, n_risk = nrow(d), n_event = 0, n_cens = 0, surv = 1), steps)
  marks <- steps |> filter(n_cens > 0) |> select(t, surv)
  list(steps = steps, marks = marks)
}

# corrupt: "none" | "KM01" (swap two adjacent survival values; non-monotonic)
#        | "KM03" (curve shifted down, so it does not start at 1.0)
#        | "KM04" (at-risk table numbers corrupted)
render_km_v2 <- function(adtte, out_png, paramcd = "OS", corrupt = "none", seed = 4) {
  d <- adtte |> filter(PARAMCD == paramcd)
  est <- lapply(names(KM_ARMS), function(a) km_estimate(filter(d, ARM == a)))
  names(est) <- names(KM_ARMS)

  if (corrupt == "none") {   # cross-check against survival::survfit
    for (a in names(est)) {
      da <- filter(d, ARM == a)
      sf <- survfit(Surv(AVAL, 1 - CNSR) ~ 1, data = da)
      ours <- est[[a]]$steps; ours <- ours[ours$t > 0, ]
      stopifnot(max(abs(ours$surv - sf$surv)) < 1e-9)
    }
  }
  for (a in names(est)) {
    s <- est[[a]]$steps
    if (corrupt == "KM01") {
      k <- which(diff(s$surv) != 0); k <- k[ceiling(length(k) / 2)]
      tmp <- s$surv[k]; s$surv[k] <- s$surv[k + 1]; s$surv[k + 1] <- tmp
    }
    if (corrupt == "KM03") s$surv <- pmax(0, s$surv - 0.08)
    est[[a]]$steps <- s
    est[[a]]$marks <- s |> filter(n_cens > 0) |> select(t, surv)
  }

  curve <- bind_rows(lapply(names(est), function(a) mutate(est[[a]]$steps, arm = a)))
  marks <- bind_rows(lapply(names(est), function(a) mutate(est[[a]]$marks, arm = a)))
  tab <- expand.grid(arm = names(KM_ARMS), t = KM_TICKS, stringsAsFactors = FALSE)
  tab$n_risk <- mapply(function(a, tm) sum(filter(d, ARM == a)$AVAL / KM_MONTH >= tm),
                       tab$arm, tab$t)
  if (corrupt == "KM04") {
    set.seed(seed)
    tab$n_risk <- pmax(0, tab$n_risk - sample(3:9, nrow(tab), replace = TRUE))
  }

  p1 <- ggplot(curve, aes(t, surv, colour = arm)) +
    geom_step(linewidth = 0.9) +
    geom_point(data = marks, aes(t, surv), shape = 3, size = 2.2, show.legend = FALSE) +
    scale_colour_manual(values = KM_ARMS, name = NULL) +
    scale_x_continuous(breaks = KM_TICKS, limits = c(0, 48)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.2)) +
    labs(x = NULL, y = "Overall survival probability",
         title = "Kaplan-Meier overall survival (+ = censored)") +
    theme_classic(base_size = 11) + theme(legend.position = "top")
  p2 <- ggplot(tab, aes(t, arm, label = n_risk, colour = arm)) +
    geom_text(size = 3.4) +
    scale_colour_manual(values = KM_ARMS, guide = "none") +
    scale_x_continuous(breaks = KM_TICKS, limits = c(0, 48)) +
    labs(x = "Months since randomization", y = NULL, subtitle = "Number at risk") +
    theme_classic(base_size = 10)
  ggsave(out_png, p1 / p2 + plot_layout(heights = c(4, 1)), width = 9, height = 7,
         dpi = 120, device = "png")

  base <- sub("\\.png$", "", out_png)
  write.csv(curve[, c("arm", "t", "surv")], paste0(base, "_curve.csv"), row.names = FALSE)
  write.csv(marks, paste0(base, "_marks.csv"), row.names = FALSE)
  write.csv(tab, paste0(base, "_risktable.csv"), row.names = FALSE)
  write.csv(d, paste0(base, "_data_used.csv"), row.names = FALSE)
  invisible(curve)
}

if (sys.nframe() == 0) {
  out <- "./Outputs/km_v2"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
  render_km_v2(read.csv("./Data/V1/ADTTE.csv", stringsAsFactors = FALSE),
               file.path(out, "KM_v2_clean.png"))
}
