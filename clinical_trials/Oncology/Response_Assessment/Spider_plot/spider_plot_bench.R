# =============================================================================
# ONCVIZ-001 -- Spider plot, benchmark implementation (switch-based)
# Differences from spider_plot.R: one panel, all patients coloured by arm, time
# in days since baseline, each patient's line stops at the earlier of first
# documented PD (ADRS) and death (ADSL OSDTC with OSCR = 1). ADTR itself keeps
# measurements after PD/death, so truncation happens here.
# Errors are switches / data injections: ref_lines (SP01), baseline offset in
# the data (SP02, runner), truncate (SP03), mixed_units (SP04).
# DRAFT NOTE: drafted with LLM assistance; needs independent review if the paper
# claims hand-written benchmark scripts.
# =============================================================================
suppressPackageStartupMessages({ library(dplyr); library(ggplot2) })

prepare_spider_bench <- function(adtr, adrs, adsl) {
  first_pd <- adrs |> filter(PARAMCD == "OVRLRESP", AVALC == "PD", !is.na(ADT), ADT != "") |>
    group_by(USUBJID) |> summarise(pd_date = min(as.Date(ADT)), .groups = "drop")
  death <- adsl |> filter(OSCR == 1) |> transmute(USUBJID, death_date = as.Date(OSDTC))
  stop_df <- adsl |> select(USUBJID, ARM) |>
    left_join(first_pd, by = "USUBJID") |> left_join(death, by = "USUBJID") |>
    mutate(stop_date = pmin(pd_date, death_date, na.rm = TRUE))
  ser <- adtr |> filter(PARAMCD == "SUMDIAM") |> mutate(date = as.Date(ADT)) |>
    group_by(USUBJID) |> mutate(day = as.numeric(date - min(date))) |> ungroup() |>
    inner_join(stop_df, by = "USUBJID") |>
    select(USUBJID, ARM = ARM.y, date, day, pchg = PCHG, stop_date)
  ser
}

render_spider_bench <- function(ser, out_png, ref_lines = TRUE, truncate = TRUE,
                                mixed_units = FALSE, n_mixed = 10) {
  if (truncate) ser <- ser |> filter(is.na(stop_date) | date <= stop_date)
  ser <- ser |> arrange(USUBJID, day)
  if (mixed_units) {   # some patients' time given in weeks while the axis says days
    ids <- sort(unique(ser$USUBJID))[seq_len(n_mixed)]
    ser$day[ser$USUBJID %in% ids] <- ser$day[ser$USUBJID %in% ids] / 7
  }
  p <- ggplot(ser, aes(day, pchg, group = USUBJID, colour = ARM)) +
    geom_line(linewidth = 0.5, alpha = 0.8) +
    geom_hline(yintercept = 0, colour = "black") +
    (if (ref_lines) geom_hline(yintercept = c(-30, 20), linetype = "dashed", colour = "#555555")) +
    (if (ref_lines) annotate("text", x = max(ser$day), y = c(-30, 20), label = c("-30% (PR)", "+20% (PD)"),
                             hjust = 1, vjust = -0.4, size = 3)) +
    scale_colour_manual(values = c(TREATMENT = "#2166AC", CONTROL = "#B2182B")) +
    labs(x = "Days since baseline", y = "% change in sum of longest diameters",
         title = "Spider plot: change in tumor burden over time (lines end at PD or death)") +
    theme_classic(base_size = 11) + theme(legend.position = "top")
  ggsave(out_png, p, width = 10, height = 6.5, dpi = 120, device = "png")
  base <- sub("\\.png$", "", out_png)
  write.csv(ser[, c("USUBJID", "ARM", "date", "day", "pchg", "stop_date")], paste0(base, "_series.csv"), row.names = FALSE)
  writeLines(c(paste("ref_lines", paste(if (ref_lines) c(-30, 20) else "none", collapse = ",")),
               paste("x_unit", "days")), paste0(base, "_meta.txt"))
  invisible(ser)
}

if (sys.nframe() == 0) {
  out <- "./Outputs/spider_bench"; dir.create(out, showWarnings = FALSE, recursive = TRUE)
  rd <- function(f) read.csv(file.path("./Data/V1", f), stringsAsFactors = FALSE)
  render_spider_bench(prepare_spider_bench(rd("ADTR.csv"), rd("ADRS.csv"), rd("ADSL.csv")),
                      file.path(out, "SP_clean.png"))
}
