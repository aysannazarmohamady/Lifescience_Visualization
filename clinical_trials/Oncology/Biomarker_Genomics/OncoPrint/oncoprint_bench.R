# Benchmark OncoPrint (ADMUT): top-10 genes x all 80 randomised patients.
# Cell colour = alteration category (priority Truncating > Missense > Amplification
# > Deletion > Other). Left labels = % of the 80 patients altered; top bar = number of
# altered genes among the shown genes. Error switches are arguments of render_onco_bench().
# DRAFT NOTE: drafted with LLM assistance.
suppressPackageStartupMessages(library(ggplot2))

ONCO_COLS <- c(Truncating = "#222222", Missense = "#2E7D32", Amplification = "#C62828",
               Deletion = "#1565C0", Other = "#9E9E9E")

categorize <- function(vc) {
  ifelse(vc %in% c("Nonsense_Mutation", "Frame_Shift_Del", "Frame_Shift_Ins", "Splice_Site"), "Truncating",
  ifelse(vc %in% c("Missense_Mutation", "In_Frame_Del"), "Missense",
  ifelse(vc == "Amplification", "Amplification",
  ifelse(vc == "Deletion", "Deletion", "Other"))))
}

prepare_onco_bench <- function(admut, adsl, top = 10) {
  samples <- sort(unique(adsl$USUBJID)); N <- length(samples)
  m <- admut
  m$cat <- categorize(m$VARIANT_CLASS)
  pr <- c(Truncating = 1, Missense = 2, Amplification = 3, Deletion = 4, Other = 5)
  m <- m[order(m$USUBJID, m$HUGO_SYMBOL, pr[m$cat]), ]
  m <- m[!duplicated(m[, c("USUBJID", "HUGO_SYMBOL")]), ]     # one cell per patient-gene
  cnt <- sort(table(m$HUGO_SYMBOL), decreasing = TRUE)
  genes <- names(cnt)[seq_len(top)]
  m <- m[m$HUGO_SYMBOL %in% genes, c("USUBJID", "HUGO_SYMBOL", "cat")]
  names(m) <- c("sample", "gene", "cat")
  list(cells = m, genes = genes, samples = samples, N = N)
}

# Standard oncoprint ordering: genes by frequency; samples sorted lexicographically by
# presence in gene order (memo sort); unaltered samples last.
onco_sample_order <- function(cells, genes, samples) {
  M <- sapply(genes, function(g) as.integer(samples %in% cells$sample[cells$gene == g]))
  rownames(M) <- samples
  o <- do.call(order, c(as.data.frame(-M), list(samples)))
  samples[o]
}

# freq_denom: denominator for left % labels (default N=80); legend_swap: swap two legend
# colours; shuffle: sample columns not sorted; topbar_off: top-bar counts +1 for some samples
render_onco_bench <- function(prep, out_png, freq_denom = NULL, legend_swap = FALSE,
                              shuffle = FALSE, topbar_off = FALSE) {
  cells <- prep$cells; genes <- prep$genes; N <- prep$N
  if (is.null(freq_denom)) freq_denom <- N
  ord <- onco_sample_order(cells, genes, prep$samples)
  if (shuffle) { set.seed(5); ord <- sample(ord) }
  cells$x <- match(cells$sample, ord)
  cells$y <- length(genes) + 1 - match(cells$gene, genes)
  cells$colour <- ONCO_COLS[cells$cat]
  nalt <- as.integer(tapply(rep(1, nrow(cells)), factor(cells$sample, levels = ord), sum))
  nalt[is.na(nalt)] <- 0
  top_shown <- nalt
  if (topbar_off) top_shown[seq(2, length(top_shown), by = 4)] <- top_shown[seq(2, length(top_shown), by = 4)] + 1
  gcount <- table(factor(cells$gene, levels = genes))
  lab <- sprintf("%s  %d%%", genes, round(100 * as.integer(gcount) / freq_denom))
  leg <- ONCO_COLS
  if (legend_swap) { leg["Missense"] <- ONCO_COLS["Amplification"]; leg["Amplification"] <- ONCO_COLS["Missense"] }
  bg <- expand.grid(x = seq_along(ord), y = seq_along(genes))
  topdf <- data.frame(x = seq_along(ord), n = top_shown)
  ng <- length(genes)
  p <- ggplot() +
    geom_tile(data = bg, aes(x, y), fill = "#EEEEEE", colour = "white", linewidth = 0.3) +
    geom_tile(data = cells, aes(x, y, fill = I(colour)), colour = "white", linewidth = 0.3, height = 0.9) +
    geom_col(data = topdf, aes(x, n * 0.4 + ng + 0.8 - 0.5 + 0 * n, group = x), alpha = 0, width = 0.01) +
    geom_rect(data = topdf, aes(xmin = x - 0.4, xmax = x + 0.4, ymin = ng + 0.8, ymax = ng + 0.8 + n * 0.35),
              fill = "grey35") +
    geom_point(data = data.frame(x = 1, y = 1, cat = factor(names(leg), levels = names(leg))),
               aes(x, y, colour = cat), alpha = 0, shape = 15, size = 4) +
    scale_colour_manual(values = leg, name = "Alteration",
                        guide = guide_legend(override.aes = list(alpha = 1))) +
    scale_y_continuous(breaks = seq_along(genes), labels = rev(lab), expand = c(0.01, 0.01)) +
    scale_x_continuous(expand = c(0.005, 0.005)) +
    labs(x = sprintf("Patients (N = %d); top bars = altered genes per patient", N), y = NULL,
         title = "OncoPrint of top-10 altered genes (ADMUT)") +
    theme_minimal(base_size = 11) +
    theme(panel.grid = element_blank(), axis.text.x = element_blank(), axis.ticks.x = element_blank(),
          legend.position = "bottom")
  ggsave(out_png, p, width = 10, height = 5.2, dpi = 150)
  base <- sub("\\.png$", "", out_png)
  write.csv(data.frame(sample = cells$sample, gene = cells$gene, cat = cells$cat,
                       x = cells$x, colour = cells$colour), paste0(base, "_cells.csv"), row.names = FALSE)
  write.csv(data.frame(sample = ord, x = seq_along(ord), top = top_shown), paste0(base, "_topbar.csv"), row.names = FALSE)
  writeLines(c(paste("freq", paste(sprintf("%s=%d", genes, round(100 * as.integer(gcount) / freq_denom)), collapse = ",")),
               paste("gene_order", paste(genes, collapse = ",")),
               paste("legend", paste(names(leg), leg, sep = "=", collapse = ","))),
             paste0(base, "_meta.txt"))
}
