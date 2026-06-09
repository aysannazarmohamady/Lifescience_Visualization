library(dplyr)
library(ggplot2)
library(patchwork)
library(grid)
library(gridExtra)

DATA_DIR   <- "./Data/V1"
OUTPUT_DIR <- "./Outputs"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

adsl <- read.csv(file.path(DATA_DIR, "ADSL.csv"), stringsAsFactors = FALSE)
adrs <- read.csv(file.path(DATA_DIR, "ADRS.csv"), stringsAsFactors = FALSE)
adsl$TRTSDT <- as.Date(adsl$TRTSDT)
adsl$TRTEDT <- as.Date(adsl$TRTEDT)
adsl$OSDTC  <- suppressWarnings(as.Date(adsl$OSDTC))

DAY2MO   <- 30.4375
BAR_BLUE <- "#1F4E79"
FU_GRAY  <- "#808080"
RESP_COLORS <- c(CR="#1A9641", PR="#4A90C4", SD="#E8A020", PD="#C0392B", NE="#999999")
ID_COLORS   <- c(CR="#1A9641", PR="#2B6CB0", SD="#B45309", PD="#C0392B", NE="#C0392B")
LM_COLORS   <- c(Y="#C0392B", N="#2B6CB0")
TUMOR_ORDER <- c("NSCLC","BRCA","HCC","CRC","PDAC")

rank_map <- c(CR=1, PR=2, SD=3, PD=4, NE=5)

best_resp <- function(uid) {
  rs <- adrs[adrs$USUBJID == uid & adrs$PARAMCD == "OVRLRESP" &
               adrs$AVALC %in% names(rank_map), ]
  if (nrow(rs) == 0) return("NE")
  rs$AVALC[which.min(rank_map[rs$AVALC])]
}

build_swimmer_data <- function(adsl_sub) {
  rows      <- list()
  resp_rows <- list()

  for (i in seq_len(nrow(adsl_sub))) {
    pt     <- adsl_sub[i, ]
    uid    <- pt$USUBJID
    short  <- sub("^[^-]+-[^-]+-", "", uid)
    trtsdt <- pt$TRTSDT
    trtedt <- pt$TRTEDT
    ong    <- pt$EOSSTT == "ONGOING"
    dc     <- trimws(pt$DCSREAS)
    bor    <- best_resp(uid)
    lm     <- trimws(pt$LIVERMETS)
    oscr   <- as.integer(pt$OSCR)

    trt_mo <- max((as.numeric(trtedt - trtsdt)) / DAY2MO, 0.3)
    fu_mo  <- NA_real_
    osdtc  <- pt$OSDTC
    if (!ong && !is.na(osdtc) && osdtc > trtedt)
      fu_mo <- as.numeric(osdtc - trtsdt) / DAY2MO

    pt_rs <- adrs[adrs$USUBJID == uid & adrs$PARAMCD == "OVRLRESP", ]
    pt_rs$ADT <- suppressWarnings(as.Date(pt_rs$ADT))
    pt_rs     <- pt_rs[!is.na(pt_rs$ADT), ]
    pt_rs$mo  <- as.numeric(pt_rs$ADT - trtsdt) / DAY2MO
    pt_rs     <- pt_rs[pt_rs$mo >= 0 & pt_rs$mo <= trt_mo + 1, ]
    pt_rs     <- pt_rs[order(pt_rs$mo), ]
    seen <- character(0)
    for (j in seq_len(nrow(pt_rs))) {
      rv <- pt_rs$AVALC[j]
      if (rv %in% names(RESP_COLORS) && !rv %in% seen) {
        resp_rows[[length(resp_rows)+1]] <- data.frame(
          uid=uid, mo=pt_rs$mo[j], resp=rv, stringsAsFactors=FALSE)
        seen <- c(seen, rv)
      }
    }

    end_t <- if (ong) "ongoing" else
              if (oscr == 1 || dc == "Death") "death" else
              if (grepl("Consent", dc)) "withdraw" else "disc"

    rows[[i]] <- data.frame(
      uid=uid, short=short, lm=lm, bor=bor,
      ong=ong, end_t=end_t, trt_mo=trt_mo, fu_mo=fu_mo,
      eval=bor!="NE", stringsAsFactors=FALSE)
  }
  list(df=bind_rows(rows),
       resp_df=if (length(resp_rows)>0) bind_rows(resp_rows) else data.frame())
}

sort_swimmer <- function(df) {
  bor_n <- c(CR=0, PR=1, SD=2, PD=3, NE=4)
  df$bor_n <- bor_n[df$bor]
  df$bor_n[is.na(df$bor_n)] <- 4
  ev <- df[df$eval, ]  |> arrange(bor_n, desc(trt_mo))
  ne <- df[!df$eval, ] |> arrange(desc(trt_mo))
  out <- bind_rows(ev, ne)
  out$y <- nrow(out):1
  list(df=out, n_eval=nrow(ev))
}

draw_swimmer <- function(adsl_sub, title_str, fname,
                         fig_w=16, fig_h=9) {

  res        <- build_swimmer_data(adsl_sub)
  df_raw     <- res$df
  resp_df    <- res$resp_df
  sorted     <- sort_swimmer(df_raw)
  df         <- sorted$df
  n_eval     <- sorted$n_eval
  n          <- nrow(df)
  n_ne       <- n - n_eval

  all_ends <- c(df$trt_mo, df$fu_mo[!is.na(df$fu_mo)])
  max_data <- max(all_ends)
  max_mo   <- ceiling((max_data + 2) / 3) * 3
  step     <- if (max_mo > 36) 6 else 3
  mstones  <- seq(step, max_mo, by = step)

  ax_h_in   <- fig_h * 0.815
  y_range   <- n + 2.0
  in_per_u  <- ax_h_in / y_range
  bar_h_pt  <- min(0.62 * in_per_u * 72, 14.0)
  BAR_FRAC  <- bar_h_pt / (in_per_u * 72)
  mk_sz_pt  <- bar_h_pt
  MK_SIZE   <- mk_sz_pt / 2.845   # ggplot size in mm ~ pt/2.845

  df$y <- as.numeric(df$y)

  bar_data <- df |>
    mutate(bar_end = trt_mo,
           fu_end  = ifelse(!is.na(fu_mo) & fu_mo > trt_mo, fu_mo, NA_real_))

  fu_data <- bar_data |> filter(!is.na(fu_end))

  resp_plot_df <- if (nrow(resp_df) > 0)
    left_join(resp_df, df[, c("uid","y")], by="uid") else NULL

  end_sym <- df |>
    mutate(
      x_e   = ifelse(!is.na(fu_mo) & fu_mo > trt_mo, fu_mo, trt_mo),
      sym_x = ifelse(end_t == "ongoing", x_e + 0.55, x_e + 0.35)
    )

  p <- ggplot() +
    geom_hline(yintercept = seq(0.5, n+1.5, 1),
               color = NA, linewidth = 0) +
    theme_classic(base_size = 10) +
    theme(
      panel.background   = element_rect(fill="white", color=NA),
      panel.grid         = element_blank(),
      axis.line          = element_blank(),
      axis.ticks.y       = element_blank(),
      axis.text.y        = element_blank(),
      axis.title.y       = element_blank(),
      axis.text.x        = element_text(size=9.5),
      axis.title.x       = element_blank(),
      plot.title         = element_text(face="bold", size=14,
                                        hjust=0.5, margin=margin(b=12)),
      legend.position    = "none",
      plot.margin        = margin(6, 6, 4, 6)
    )

  for (y_even in df$y[df$y %% 2 == 0]) {
    p <- p + annotate("rect",
      xmin=-Inf, xmax=Inf,
      ymin=y_even-0.5, ymax=y_even+0.5,
      fill="#EFF4FB", alpha=0.55)
  }

  for (mo in mstones) {
    p <- p +
      annotate("segment",
               x=mo, xend=mo, y=0.2, yend=n+2,
               color="#BBBBBB", linewidth=0.8, linetype="dashed") +
      annotate("text",
               x=mo, y=n+1.55,
               label=as.character(mo),
               hjust=0.5, vjust=0, size=2.7,
               fontface="bold.italic", color="#444444")
  }

  p <- p +
    geom_tile(data=bar_data,
              aes(x=trt_mo/2, y=y, width=trt_mo, height=BAR_FRAC),
              fill=BAR_BLUE, alpha=0.93)

  if (nrow(fu_data) > 0) {
    fu_data2 <- fu_data |>
      mutate(fu_w = fu_end - trt_mo,
             fu_cx = trt_mo + fu_w/2)
    p <- p +
      geom_tile(data=fu_data2,
                aes(x=fu_cx, y=y, width=fu_w, height=BAR_FRAC*0.42),
                fill=FU_GRAY, alpha=0.75)
  }

  if (!is.null(resp_plot_df) && nrow(resp_plot_df) > 0) {
    mk_shapes <- c(CR=16, PR=15, SD=17, PD=18)
    for (rv in intersect(names(mk_shapes), unique(resp_plot_df$resp))) {
      sub_r <- resp_plot_df[resp_plot_df$resp == rv, ]
      p <- p + geom_point(data=sub_r,
                           aes(x=mo, y=y),
                           shape=mk_shapes[rv],
                           size=MK_SIZE,
                           color=RESP_COLORS[rv],
                           stroke=0.4)
    }
  }

  ong_df  <- end_sym[end_sym$end_t == "ongoing", ]
  dth_df  <- end_sym[end_sym$end_t == "death",   ]
  disc_df <- end_sym[end_sym$end_t %in% c("withdraw","disc"), ]

  if (nrow(ong_df) > 0)
    p <- p + geom_point(data=ong_df, aes(x=sym_x, y=y),
                        shape=17, size=MK_SIZE*1.05,
                        color=BAR_BLUE)
  if (nrow(dth_df) > 0)
    p <- p + geom_point(data=dth_df, aes(x=sym_x, y=y),
                        shape=4, size=MK_SIZE*0.92,
                        color="black", stroke=1.8)
  if (nrow(disc_df) > 0)
    p <- p + geom_point(data=disc_df, aes(x=sym_x, y=y),
                        shape=4, size=MK_SIZE*0.92,
                        color="#C0392B", stroke=1.8)

  if (n_ne > 0) {
    sep_y <- n_eval + 0.5
    p <- p + geom_hline(yintercept=sep_y, linewidth=1.3, color="black")
  }

  p <- p +
    scale_x_continuous(breaks=mstones,
                       labels=as.character(mstones),
                       limits=c(0, max_mo+0.5),
                       expand=c(0,0)) +
    scale_y_continuous(limits=c(0.2, n+2.0), expand=c(0,0)) +
    labs(title=title_str)

  summary_txt <- sprintf(
    "Summary\nEvaluable       %d\nNon-evaluable   %d\nTotal Enrolled  %d",
    n_eval, n_ne, n)

  p <- p + annotate("label",
    x=max_mo*0.97, y=1.8,
    label=summary_txt,
    hjust=1, vjust=0,
    size=3.4, family="mono",
    fill=BAR_BLUE, color="white",
    label.r=unit(0.25,"lines"),
    label.padding=unit(0.45,"lines"),
    alpha=0.9)

  metadata <- df |>
    mutate(
      id_x     = -0.5,
      recist_x = -1.5,
      lm_x     = -2.6
    )

  p <- p +
    geom_text(data=metadata,
              aes(x=id_x, y=y, label=short,
                  color=bor),
              hjust=1, size=2.4, family="mono") +
    geom_text(data=metadata,
              aes(x=recist_x, y=y, label=bor,
                  color=bor),
              hjust=0.5, size=2.8, fontface="bold") +
    geom_text(data=metadata,
              aes(x=lm_x, y=y, label=lm),
              hjust=0.5, size=2.8, fontface="bold",
              color=ifelse(metadata$lm=="Y", LM_COLORS["Y"], LM_COLORS["N"])) +
    scale_color_manual(values=ID_COLORS) +
    annotate("text", x=-1.5, y=n+1.55,
             label="RECIST\nResponse", hjust=0.5, vjust=0,
             size=2.5, fontface="bold", color="#111111", lineheight=1.2) +
    annotate("text", x=-2.6, y=n+1.55,
             label="Liver\nmets", hjust=0.5, vjust=0,
             size=2.5, fontface="bold", color="#111111", lineheight=1.2)

  legend_df <- data.frame(
    x     = rep(1:9, each=1),
    label = c("On Treatment", "Off Treatment / Follow-Up",
              "Still on Treatment",
              "CR","PR","SD","PD","Death","Discontinuation"),
    type  = c("bar_blue","bar_gray","arrow",
              "CR","PR","SD","PD","death","disc"),
    stringsAsFactors = FALSE
  )

  leg_p <- ggplot() +
    annotate("rect", xmin=0.05, xmax=0.18, ymin=0.62, ymax=0.72,
             fill=BAR_BLUE, alpha=0.93) +
    annotate("text", x=0.21, y=0.67, label="On Treatment",
             hjust=0, size=3.3) +
    annotate("rect", xmin=0.05, xmax=0.18, ymin=0.38, ymax=0.48,
             fill=FU_GRAY, alpha=0.75) +
    annotate("text", x=0.21, y=0.43,
             label="Off Treatment and On Study Follow-Up",
             hjust=0, size=3.3) +
    annotate("point", x=0.115, y=0.19, shape=17, size=3.5, color=BAR_BLUE) +
    annotate("text",  x=0.21,  y=0.19, label="Still on Treatment",
             hjust=0, size=3.3) +
    annotate("point", x=0.62, y=0.87, shape=16, size=3.5,
             color=RESP_COLORS["CR"]) +
    annotate("text",  x=0.66, y=0.87, label="CR", hjust=0, size=3.3) +
    annotate("point", x=0.72, y=0.87, shape=15, size=3.2,
             color=RESP_COLORS["PR"]) +
    annotate("text",  x=0.76, y=0.87, label="PR", hjust=0, size=3.3) +
    annotate("point", x=0.82, y=0.87, shape=17, size=3.2,
             color=RESP_COLORS["SD"]) +
    annotate("text",  x=0.86, y=0.87, label="SD", hjust=0, size=3.3) +
    annotate("point", x=0.92, y=0.87, shape=18, size=3.5,
             color=RESP_COLORS["PD"]) +
    annotate("text",  x=0.96, y=0.87, label="PD", hjust=0, size=3.3) +
    annotate("point", x=0.62, y=0.62, shape=4, size=3.5,
             color="black", stroke=1.8) +
    annotate("text",  x=0.66, y=0.62, label="Death", hjust=0, size=3.3) +
    annotate("point", x=0.62, y=0.38, shape=4, size=3.5,
             color="#C0392B", stroke=1.8) +
    annotate("text",  x=0.66, y=0.38, label="Discontinuation",
             hjust=0, size=3.3) +
    xlim(0, 1.15) + ylim(0, 1) +
    theme_void() +
    theme(plot.margin=margin(2,2,2,2))

  png(fname, width=fig_w, height=fig_h, units="in", res=120)
  print(
    p + inset_element(leg_p,
                      left=0.78, right=1.02,
                      bottom=0.10, top=0.90,
                      align_to="plot")
  )
  dev.off()
  message("  Saved: ", basename(fname))
}

trt <- adsl[adsl$ARM == "TREATMENT", ]

draw_swimmer(trt,
  "Swimmer Plot: ONCVIZ-001 \u00b7 Vizatinib 300mg QD \u2014 All Patients  \u00b7  Treatment Arm",
  file.path(OUTPUT_DIR, "swimmer_all_treatment.png"),
  fig_w=16, fig_h=11)

for (tumor in TUMOR_ORDER) {
  sub <- trt[trt$TUMORTYPE == tumor, ]
  draw_swimmer(sub,
    sprintf("Swimmer Plot: ONCVIZ-001 \u00b7 Vizatinib 300mg QD \u2014 All Patients  \u00b7  %s", tumor),
    file.path(OUTPUT_DIR, sprintf("swimmer_%s.png", tolower(tumor))),
    fig_w=16, fig_h=9)
}

message("All swimmer plots saved to: ", OUTPUT_DIR)
