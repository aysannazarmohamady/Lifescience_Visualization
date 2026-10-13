# Benchmark cumulative incidence plot (ADTTE PFS rows): cumulative incidence of
# progression with death-without-progression as competing event, by arm (cmprsk).
# Error switches are arguments of render_cif_bench().
# DRAFT NOTE: drafted with LLM assistance.
suppressPackageStartupMessages({library(ggplot2); library(cmprsk)})

prepare_cif_bench <- function(adtte) {
  d <- adtte[adtte$PARAMCD == "PFS", ]
  d$status <- ifelse(d$COMPTYPE == "PROGRESSION", 1L, ifelse(d$COMPTYPE == "DEATH_WITHOUT_PROGRESSION", 2L, 0L))
  d$t <- d$AVALM
  d
}

# cuminc()'s stored (time, est) pairs are lagged for step plotting, so evaluate the
# estimate at every distinct failure time with timepoints() (right-continuous steps).
cif_tp <- function(g, key) {
  tt <- c(0, sort(unique(g$t[g$status > 0])))
  fit <- cuminc(g$t, g$status, cencode = 0)
  if (is.null(fit[[key]])) return(list(time = c(0, max(g$t)), est = c(0, 0)))   # no such events in arm
  list(time = tt, est = as.numeric(timepoints(fit, tt)$est[key, ]))
}

cif_curves <- function(d, kind = "cif") {
  out <- list()
  for (a in c("TREATMENT", "CONTROL")) {
    g <- d[d$ARM == a, ]
    if (kind == "cif") {
      ci <- cif_tp(g, "1 1")
    } else if (kind == "death") {
      ci <- cif_tp(g, "1 2")
    } else {                                          # 1 - KM treating competing event as censored
      sf <- survfit_km(g$t, g$status == 1)
      ci <- list(time = sf$time, est = 1 - sf$surv)
    }
    if (kind == "km") { ci$time <- c(0, ci$time); ci$est <- c(0, ci$est) }   # cuminc output already starts at 0
    out[[a]] <- data.frame(arm = a, time = ci$time, est = ci$est)
  }
  do.call(rbind, out)
}
survfit_km <- function(t, ev) {
  u <- sort(unique(t[ev])); s <- 1; surv <- numeric(length(u))
  for (i in seq_along(u)) { n <- sum(t >= u[i]); e <- sum(t == u[i] & ev); s <- s * (1 - e / n); surv[i] <- s }
  list(time = u, surv = surv)
}

# one_minus_km: 1-KM instead of CIF; wrong_event: death CIF labelled as progression;
# days_axis: x in days but axis says months; no_risk: omit number-at-risk table
render_cif_bench <- function(d, out_png, one_minus_km = FALSE, wrong_event = FALSE,
                             days_axis = FALSE, no_risk = FALSE) {
  kind <- if (one_minus_km) "km" else if (wrong_event) "death" else "cif"
  cv <- cif_curves(d, kind)
  if (days_axis) cv$time <- cv$time * 30.4375
  xmax <- 42 * (if (days_axis) 30.4375 else 1)
  brk <- if (days_axis) seq(0, 1260, by = 180) else seq(0, 42, by = 6)
  risk <- do.call(rbind, lapply(c("TREATMENT", "CONTROL"), function(a)
    data.frame(arm = a, t = seq(0, 42, by = 6), n = sapply(seq(0, 42, by = 6), function(tt) sum(d$ARM == a & d$t >= tt)))))
  p <- ggplot(cv, aes(time, est, colour = arm)) +
    geom_step(linewidth = 0.9) +
    scale_colour_manual(values = c(TREATMENT = "#1F77B4", CONTROL = "#D62728"), name = NULL) +
    scale_x_continuous(breaks = brk, labels = brk, limits = c(0, xmax)) +
    scale_y_continuous(breaks = seq(0, 1, 0.25)) + coord_cartesian(ylim = c(0, 1), clip = "off") +
    labs(x = "Time since randomisation (months)",
         y = "Cumulative incidence of progression",
         title = "Cumulative incidence of progression (competing risk: death without progression)") +
    theme_bw(base_size = 11) + theme(legend.position = "top", plot.title = element_text(size = 10))
  if (!no_risk) {
    rr <- risk; rr$x <- rr$t * (if (days_axis) 30.4375 else 1); rr$y <- ifelse(rr$arm == "TREATMENT", -0.27, -0.34)
    p <- p + geom_text(data = rr, aes(x, y, label = n, colour = arm), show.legend = FALSE, size = 3.2) +
      annotate("text", x = -Inf, y = -0.20, label = "Number at risk", hjust = -0.02, size = 3.2) +
      theme(plot.margin = margin(5, 10, 70, 5))
  }
  ggsave(out_png, p, width = 7.5, height = 5.5, dpi = 150)
  base <- sub("\\.png$", "", out_png)
  write.csv(cv, paste0(base, "_curves.csv"), row.names = FALSE)
  write.csv(if (no_risk) risk[0, ] else risk, paste0(base, "_risk.csv"), row.names = FALSE)
  writeLines(c(paste("x_unit", if (days_axis) "days" else "months"), paste("x_label", "months"),
               paste("event_label", "progression")), paste0(base, "_meta.txt"))
}
