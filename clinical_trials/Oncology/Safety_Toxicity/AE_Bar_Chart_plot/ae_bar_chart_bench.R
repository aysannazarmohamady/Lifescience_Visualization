# Benchmark AE bar chart: % of treatment-arm patients by preferred term, stacked by
# worst CTCAE grade per patient. Error switches are arguments of render_ae_bench().
# DRAFT NOTE: drafted with LLM assistance.
suppressPackageStartupMessages(library(ggplot2))

GRADE_COLS <- c(`1` = "#FEE8C8", `2` = "#FDBB84", `3` = "#E34A33", `4` = "#B30000", `5` = "#4D0000")

prepare_ae_bench <- function(adae, adsl, arm = "TREATMENT", top = 12) {
  n_arm <- sum(adsl$ARM == arm)
  a <- adae[adae$ARM == arm, ]
  w <- aggregate(AETOXGR ~ USUBJID + AEPT, a, max)         # worst grade per patient & term
  tab <- as.data.frame(table(term = w$AEPT, grade = w$AETOXGR), stringsAsFactors = FALSE)
  tab <- tab[tab$Freq > 0, ]
  tot <- tapply(tab$Freq, tab$term, sum)
  keep <- names(sort(tot, decreasing = TRUE))[seq_len(top)]
  tab <- tab[tab$term %in% keep, ]
  tab$n_arm <- n_arm
  tab
}

# bad_colours: grade colours reversed; drop_g45: grade 4-5 not drawn;
# denom: denominator used for percentages (default = arm N);
# shuffle: terms not ordered by frequency
render_ae_bench <- function(tab, out_png, bad_colours = FALSE, drop_g45 = FALSE,
                            denom = NULL, shuffle = FALSE) {
  n_arm <- tab$n_arm[1]; if (is.null(denom)) denom <- n_arm
  cols <- GRADE_COLS; if (bad_colours) cols <- setNames(rev(unname(GRADE_COLS)), names(GRADE_COLS))
  d <- tab; if (drop_g45) d <- d[!d$grade %in% c("4", "5"), ]
  d$pct <- 100 * d$Freq / denom
  tot <- tapply(d$pct, d$term, sum)
  lv <- names(sort(tot))                                    # ascending => largest on top
  if (shuffle) { set.seed(3); lv <- sample(lv) }
  d$term <- factor(d$term, levels = lv)
  p <- ggplot(d, aes(pct, term, fill = grade)) +
    geom_col(width = 0.7, colour = "grey30", linewidth = 0.2) +
    scale_fill_manual(values = cols, name = "Worst CTCAE grade") +
    labs(x = sprintf("Patients with event (%%), N = %d", n_arm), y = NULL,
         title = "Adverse events by preferred term and worst grade, treatment arm") +
    theme_bw(base_size = 12) + theme(legend.position = "bottom")
  ggsave(out_png, p, width = 7.5, height = 5.5, dpi = 150)
  base <- sub("\\.png$", "", out_png)
  o <- data.frame(term = as.character(d$term), grade = d$grade, pct = d$pct,
                  order = as.integer(d$term))
  write.csv(o, paste0(base, "_bars.csv"), row.names = FALSE)
  writeLines(c(paste("n_label", n_arm),
               paste("colours", paste(names(cols), cols, sep = "=", collapse = ","))),
             paste0(base, "_meta.txt"))
}
