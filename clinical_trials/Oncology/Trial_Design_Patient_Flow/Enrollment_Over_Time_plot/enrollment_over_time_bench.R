# =============================================================================
# ONCVIZ-001 -- Enrollment over time, benchmark implementation (switch-based)
# Differences from enrollment_over_time_plot.R: overall + per-arm cumulative
# step curves, a vertical line at the study closure (data cut-off, ADSL$CUTDTC),
# a horizontal target-enrollment line (planned N = randomized N, a proxy because
# the source data hold no planned-accrual figure), and arm totals in the legend.
# Errors are switches: monotone (EN01), arm_label_shift (EN02), data injection
# of enrollment dates after closure (EN03, in the runner), target_line (EN04).
# DRAFT NOTE: drafted with LLM assistance; needs independent review if the paper
# claims hand-written benchmark scripts.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

prepare_enrollment_bench <- function(adrand, adsl) {
  rand <- adrand |> filter(RANDFL == "Y") |>
    inner_join(adsl |> select(USUBJID, TRTSDT, ARM, CUTDTC), by = "USUBJID", suffix = c("", ".sl")) |>
    mutate(date = as.Date(TRTSDT)) |> arrange(date, USUBJID)
  rand$cum <- seq_len(nrow(rand))
  rand$arm <- rand$ARM
  list(rand = rand, closure = as.Date(unique(adsl$CUTDTC)[1]), n = nrow(rand))
}

render_enrollment_bench <- function(prep, out_png, monotone = TRUE, arm_label_shift = 0,
                                    target_line = TRUE) {
  r <- prep$rand
  overall <- data.frame(date = r$date, cum = r$cum)
  if (!monotone) {
    k <- floor(nrow(overall) / 2)
    tmp <- overall$cum[k]; overall$cum[k] <- overall$cum[k + 1] - 6; overall$cum[k + 1] <- overall$cum[k + 1]
  }
  arms <- r |> group_by(arm) |> mutate(cum_arm = row_number()) |> ungroup()
  n_arm <- table(r$arm)
  n_shown <- c(TREATMENT = as.integer(n_arm["TREATMENT"]), CONTROL = as.integer(n_arm["CONTROL"]))
  n_shown["CONTROL"] <- n_shown["CONTROL"] + arm_label_shift
  lab <- c(TREATMENT = sprintf("Treatment (n=%d)", n_shown["TREATMENT"]),
           CONTROL   = sprintf("Control (n=%d)", n_shown["CONTROL"]))
  arms$arm_lab <- lab[arms$arm]
  title <- sprintf("Cumulative enrollment: overall n=%d (Treatment %d + Control %d)",
                   prep$n, n_shown["TREATMENT"], n_shown["CONTROL"])

  p <- ggplot() +
    geom_step(data = overall, aes(date, cum), linewidth = 1.1, colour = "#08306B") +
    geom_step(data = arms, aes(date, cum_arm, colour = arm_lab), linewidth = 0.7) +
    geom_vline(xintercept = prep$closure, linetype = "dashed", colour = "#B2182B") +
    annotate("text", x = prep$closure, y = 5, label = "Study closure (data cut-off)",
             angle = 90, vjust = -0.4, hjust = 0, colour = "#B2182B", size = 3) +
    (if (target_line) geom_hline(yintercept = prep$n, linetype = "dotted", colour = "#444444")) +
    (if (target_line) annotate("text", x = min(r$date), y = prep$n, vjust = -0.5, hjust = 0,
                               label = sprintf("Target enrollment (planned N = %d)", prep$n), size = 3)) +
    scale_colour_manual(values = c("#2166AC", "#B2182B"), name = NULL) +
    labs(x = "Randomization date (first-dose date proxy)", y = "Cumulative patients randomized",
         title = title) +
    theme_classic(base_size = 11) + theme(legend.position = "top")
  ggsave(out_png, p, width = 10, height = 6, dpi = 120, device = "png")

  base <- sub("\\.png$", "", out_png)
  write.csv(overall, paste0(base, "_curve.csv"), row.names = FALSE)
  writeLines(c(paste("closure", prep$closure), paste("target_line", target_line),
               paste("title", title)), paste0(base, "_meta.txt"))
  invisible(NULL)
}

if (sys.nframe() == 0) {
  out <- "./Outputs/enrollment_bench"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
  prep <- prepare_enrollment_bench(read.csv("./Data/V1/ADRAND.csv", stringsAsFactors = FALSE),
                                   read.csv("./Data/V1/ADSL.csv", stringsAsFactors = FALSE))
  render_enrollment_bench(prep, file.path(out, "EN_clean.png"))
}
