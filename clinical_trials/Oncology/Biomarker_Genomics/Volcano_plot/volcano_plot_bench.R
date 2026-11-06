# Benchmark volcano plot of on-treatment biomarker change (ADBM, treatment arm).
# Each point = one biomarker x arm; x = mean over patients of the per-patient mean
# log2(AVAL/BASE) over all post-baseline visits, y = -log10(p) from a one-sample
# t-test of those per-patient values (exploratory, unadjusted; synthetic data has
# too little signal for a multiplicity-adjusted threshold).
# Error switches are arguments of render_volcano_bench().
# DRAFT NOTE: drafted with LLM assistance.
suppressPackageStartupMessages(library(ggplot2))

prepare_volcano_bench <- function(adbm, alpha = 0.05) {
  d <- adbm[adbm$AVISIT != "BASELINE" & !is.na(adbm$AVAL) & !is.na(adbm$BASE) &
              adbm$AVAL > 0 & adbm$BASE > 0, ]
  d$l2 <- log2(d$AVAL / d$BASE)
  d$feat <- paste(substr(d$ARM, 1, 3), d$PARAMCD)
  rows <- lapply(split(d, d$feat), function(g) {
    pat <- tapply(g$l2, g$USUBJID, mean)
    if (length(pat) < 10) return(NULL)
    data.frame(feature = g$feat[1], n = length(pat), log2fc = mean(pat),
               ratio = mean(2^pat), p = t.test(pat)$p.value, stringsAsFactors = FALSE)
  })
  r <- do.call(rbind, rows)
  r$sig <- r$p < alpha
  attr(r, "alpha") <- alpha
  r
}

# thr_alpha: alpha used to draw the threshold line (clean = real alpha)
# ln_p: plot -ln(p) although axis says log10; fc_raw: plot ratio, not log2
# wrong_highlight: colour n non-significant points as significant
render_volcano_bench <- function(r, out_png, thr_alpha = NULL, ln_p = FALSE,
                                 fc_raw = FALSE, wrong_highlight = 0) {
  alpha <- attr(r, "alpha"); if (is.null(thr_alpha)) thr_alpha <- alpha
  r$x <- if (fc_raw) r$ratio else r$log2fc
  r$y <- (if (ln_p) -log(r$p) else -log10(r$p))
  r$hl <- r$sig
  if (wrong_highlight > 0) {
    ns <- which(!r$sig); ns <- ns[order(-r$p[ns])][seq_len(wrong_highlight)]  # least significant
    r$hl[ns] <- TRUE
  }
  thr_y <- -log10(thr_alpha)
  lab <- r[r$hl, ]
  p <- ggplot(r, aes(x, y, colour = hl)) +
    geom_hline(yintercept = thr_y, linetype = "dashed", colour = "grey40") +
    geom_vline(xintercept = if (fc_raw) 1 else 0, colour = "grey70") +
    geom_point(size = 2.2) +
    scale_colour_manual(values = c(`FALSE` = "grey60", `TRUE` = "#C0392B"),
                        labels = c("Not significant", "Significant (p < 0.05)"), name = NULL) +
    labs(x = "log2 fold change (on-treatment / baseline)",
         y = "-log10 (p-value)",
         title = "Biomarker change from baseline by arm (exploratory, unadjusted)",
         subtitle = "Dashed line: p = 0.05") +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  if (nrow(lab) > 0)
    p <- p + geom_text(data = lab, aes(label = feature), size = 2.6, vjust = -0.8,
                       check_overlap = TRUE, show.legend = FALSE)
  ggsave(out_png, p, width = 7.5, height = 5.5, dpi = 150)
  base <- sub("\\.png$", "", out_png)
  write.csv(data.frame(feature = r$feature, x = r$x, y = r$y, highlighted = r$hl),
            paste0(base, "_points.csv"), row.names = FALSE)
  writeLines(c(paste("thr_y", thr_y), paste("xlab", "log2 fold change"),
               paste("ylab", "log10")), paste0(base, "_meta.txt"))
}
