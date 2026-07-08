suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
})

DATA_DIR   <- "./Data/V1"
OUTPUT_DIR <- "./Out"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

adrand <- read.csv(file.path(DATA_DIR, "ADRAND.csv"), stringsAsFactors = FALSE)
adsl   <- read.csv(file.path(DATA_DIR, "ADSL.csv"),   stringsAsFactors = FALSE)

n_screened <- nrow(adrand)
n_failed   <- sum(adrand$SFFL == "Y")
sf_reasons <- adrand %>% filter(SFFL == "Y") %>% count(SFREASN, sort = TRUE)
n_rand     <- sum(adrand$RANDFL == "Y")

n_trt      <- sum(adsl$ARM == "TREATMENT")
n_ctrl     <- sum(adsl$ARM == "CONTROL")
n_trt_recv <- sum(adsl$ARM == "TREATMENT" & !is.na(adsl$TRTSDT) & adsl$TRTSDT != "")
n_ctrl_recv<- sum(adsl$ARM == "CONTROL"   & !is.na(adsl$TRTSDT) & adsl$TRTSDT != "")

disc_trt  <- adsl %>% filter(ARM == "TREATMENT", EOSSTT == "DISCONTINUED")
disc_ctrl <- adsl %>% filter(ARM == "CONTROL",   EOSSTT == "DISCONTINUED")
disc_trt_reasons  <- disc_trt  %>% count(DCSREAS, sort = TRUE)
disc_ctrl_reasons <- disc_ctrl %>% count(DCSREAS, sort = TRUE)

ittfl_trt  <- sum(adsl$ARM == "TREATMENT" & adsl$ITTFL == "Y")
ittfl_ctrl <- sum(adsl$ARM == "CONTROL"   & adsl$ITTFL == "Y")
saffl_trt  <- sum(adsl$ARM == "TREATMENT" & adsl$SAFFL == "Y")
saffl_ctrl <- sum(adsl$ARM == "CONTROL"   & adsl$SAFFL == "Y")

reasons_txt <- function(df, name_col) {
  paste(sprintf("\u2022 %s: %d", df[[name_col]], df$n), collapse = "\n")
}

# ── Box + arrow helper geometry (grid coordinates, 0-100 canvas) ──────────────
CX <- 50
box_df <- function(x, y, w, h, label) data.frame(xmin = x, xmax = x + w, ymin = y, ymax = y + h, label = label)

boxes <- list()
arrows <- list()

# 1) Assessed for eligibility
b1 <- box_df(30, 90, 40, 8, sprintf("Assessed for eligibility\n(n = %d)", n_screened))
boxes[[1]] <- b1

# 2) Excluded -- centered directly below box 1, above box 3
gap_a <- 6; excl_h <- 10
excl_top <- b1$ymin - gap_a
b2 <- box_df(CX - 19, excl_top - excl_h, 38, excl_h,
             sprintf("Excluded (n = %d)\n%s", n_failed, reasons_txt(sf_reasons, "SFREASN")))
boxes[[2]] <- b2
arrows[[1]] <- data.frame(x = CX, y = b1$ymin, xend = CX, yend = b2$ymax)

# 3) Randomized -- centered directly below Excluded
gap <- 6; h2 <- 7
rand_top <- b2$ymin - gap
b3 <- box_df(CX - 22, rand_top - h2, 44, h2, sprintf("Randomized (n = %d)", n_rand))
boxes[[3]] <- b3
arrows[[2]] <- data.frame(x = CX, y = b2$ymin, xend = CX, yend = b3$ymax)

# 4) Allocation
alloc_gap <- 8; h3 <- 10
alloc_top <- b3$ymin - alloc_gap
bL <- box_df(5, alloc_top - h3, 40, h3, sprintf("Allocated to Treatment (n = %d)\nReceived intervention (n = %d)", n_trt, n_trt_recv))
bR <- box_df(55, alloc_top - h3, 40, h3, sprintf("Allocated to Control/Placebo (n = %d)\nReceived intervention (n = %d)", n_ctrl, n_ctrl_recv))
boxes[[4]] <- bL; boxes[[5]] <- bR
arrows[[3]] <- data.frame(x = CX - 11, y = b3$ymin, xend = (bL$xmin + bL$xmax) / 2, yend = bL$ymax)
arrows[[4]] <- data.frame(x = CX + 11, y = b3$ymin, xend = (bR$xmin + bR$xmax) / 2, yend = bR$ymax)

# 5) Discontinued
disc_gap <- 8; h4 <- 14
disc_top <- bL$ymin - disc_gap
bDL <- box_df(5, disc_top - h4, 40, h4, sprintf("Discontinued intervention (n = %d)\n%s", nrow(disc_trt), reasons_txt(disc_trt_reasons, "DCSREAS")))
bDR <- box_df(55, disc_top - h4, 40, h4, sprintf("Discontinued intervention (n = %d)\n%s", nrow(disc_ctrl), reasons_txt(disc_ctrl_reasons, "DCSREAS")))
boxes[[6]] <- bDL; boxes[[7]] <- bDR
arrows[[5]] <- data.frame(x = (bL$xmin + bL$xmax) / 2, y = bL$ymin, xend = (bDL$xmin + bDL$xmax) / 2, yend = bDL$ymax)
arrows[[6]] <- data.frame(x = (bR$xmin + bR$xmax) / 2, y = bR$ymin, xend = (bDR$xmin + bDR$xmax) / 2, yend = bDR$ymax)

# 6) Analysis
an_gap <- 8; h5 <- 10
an_top <- bDL$ymin - an_gap
bAL <- box_df(5, an_top - h5, 40, h5, sprintf("Analyzed, ITT (n = %d)\nAnalyzed, Safety (n = %d)", ittfl_trt, saffl_trt))
bAR <- box_df(55, an_top - h5, 40, h5, sprintf("Analyzed, ITT (n = %d)\nAnalyzed, Safety (n = %d)", ittfl_ctrl, saffl_ctrl))
boxes[[8]] <- bAL; boxes[[9]] <- bAR
arrows[[7]] <- data.frame(x = (bDL$xmin + bDL$xmax) / 2, y = bDL$ymin, xend = (bAL$xmin + bAL$xmax) / 2, yend = bAL$ymax)
arrows[[8]] <- data.frame(x = (bDR$xmin + bDR$xmax) / 2, y = bDR$ymin, xend = (bAR$xmin + bAR$xmax) / 2, yend = bAR$ymax)

boxes_df  <- bind_rows(boxes, .id = "id")
arrows_df <- bind_rows(arrows)

p <- ggplot() +
  geom_rect(data = boxes_df, aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
            fill = "white", color = "black", linewidth = 0.6) +
  geom_text(data = boxes_df, aes(x = (xmin + xmax) / 2, y = (ymin + ymax) / 2, label = label),
            size = 3.1, lineheight = 0.95) +
  geom_segment(data = arrows_df, aes(x = x, y = y, xend = xend, yend = yend),
               arrow = arrow(length = unit(0.12, "inches"), type = "closed"), linewidth = 0.55) +
  coord_cartesian(xlim = c(0, 100), ylim = c(min(bAL$ymin, bAR$ymin) - 6, 100)) +
  theme_void() +
  labs(title = "CONSORT 2010 Participant Flow Diagram \u2014 ONCVIZ-001",
       caption = "Reasons for exclusion, discontinuation, and analysis population flags per protocol-defined definitions (SAFFL/ITTFL).") +
  theme(plot.title = element_text(hjust = 0.5, size = 14, face = "bold", margin = margin(b = 10)),
        plot.caption = element_text(hjust = 0.5, size = 8, face = "italic"))

ggsave(file.path(OUTPUT_DIR, "consort_diagram.png"), p, width = 11, height = 16, dpi = 300)
