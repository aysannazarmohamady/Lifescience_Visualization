# Benchmark PK concentration-time plot (ADPK, CONC): geometric mean (positive values
# only; zeros = below quantification) by nominal time and cycle, log10 y axis.
# Error switches are arguments of render_pk_bench().
# DRAFT NOTE: drafted with LLM assistance.
suppressPackageStartupMessages(library(ggplot2))

prepare_pk_bench <- function(adpk) {
  d <- adpk[adpk$PARAMCD == "CONC" & adpk$NOMTPT > 0, ]
  d$cycle <- paste0("Cycle ", d$AVISITN)
  g <- split(d, list(d$cycle, d$NOMTPT), drop = TRUE)
  do.call(rbind, lapply(g, function(x) {
    pos <- x$AVAL[x$AVAL > 0]
    data.frame(cycle = x$cycle[1], time = x$NOMTPT[1], gmean = exp(mean(log(pos))),
               amean = mean(x$AVAL), n_pos = length(pos), n = nrow(x))
  }))
}

# linear_y: linear axis (log expected); arith: arithmetic mean (incl. zeros) plotted as
# "geometric mean"; index_x: equally spaced sampling index instead of hours;
# swap_labels: cycle labels swapped
render_pk_bench <- function(s, out_png, linear_y = FALSE, arith = FALSE,
                            index_x = FALSE, swap_labels = FALSE) {
  s$y <- if (arith) s$amean else s$gmean
  s$x <- if (index_x) match(s$time, sort(unique(s$time))) else s$time
  s$lab <- s$cycle
  if (swap_labels) s$lab <- ifelse(s$cycle == "Cycle 1", "Cycle 3", "Cycle 1")
  p <- ggplot(s, aes(x, y, colour = lab, group = lab)) + geom_line() + geom_point(size = 2) +
    scale_colour_manual(values = c(`Cycle 1` = "#1F77B4", `Cycle 3` = "#D62728"), name = NULL) +
    labs(x = "Nominal time after dose (h)", y = "Geometric mean concentration (ng/mL)",
         title = "Vizatinib plasma concentration-time profile (n = 62)",
         subtitle = "Zero (below-quantification) values excluded from geometric means") +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  if (!linear_y) p <- p + scale_y_log10()
  ggsave(out_png, p, width = 7.5, height = 5.2, dpi = 150)
  base <- sub("\\.png$", "", out_png)
  write.csv(data.frame(label = s$lab, x = s$x, y = s$y), paste0(base, "_points.csv"), row.names = FALSE)
  writeLines(c(paste("y_scale", if (linear_y) "linear" else "log10"),
               paste("xlab", "hours")), paste0(base, "_meta.txt"))
}
