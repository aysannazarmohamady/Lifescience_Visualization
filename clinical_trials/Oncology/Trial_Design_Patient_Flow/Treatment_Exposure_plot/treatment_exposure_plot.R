library(dplyr); library(ggplot2)

DATA_DIR   <- "./Data/V1"
OUTPUT_DIR <- "./clinical_trials/Oncology/Trial_Design_Patient_Flow/Treatment_Exposure_plot/output"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

adsl <- read.csv(file.path(DATA_DIR, "ADSL.csv"), stringsAsFactors = FALSE)
adrs <- read.csv(file.path(DATA_DIR, "ADRS.csv"), stringsAsFactors = FALSE)

DAY2MO <- 30.44

trt <- adsl %>%
  filter(ARM == "TREATMENT") %>%
  transmute(USUBJID, TRTDURD, EOSSTT, BESTRSPC,
            dur_mo = TRTDURD / DAY2MO,
            ongoing = EOSSTT == "ONGOING")

ovrl <- adrs %>%
  filter(PARAMCD == "OVRLRESP") %>%
  mutate(ADTN = suppressWarnings(as.numeric(ADTN)),
         mo = ADTN / DAY2MO) %>%
  filter(!is.na(mo))

first_resp <- ovrl %>%
  filter(AVALC %in% c("CR","PR")) %>%
  group_by(USUBJID) %>%
  slice_min(mo, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(USUBJID, resp_mo = mo)

first_pd <- ovrl %>%
  filter(AVALC == "PD") %>%
  group_by(USUBJID) %>%
  slice_min(mo, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  select(USUBJID, pd_mo = mo)

trt <- trt %>%
  left_join(first_resp, by = "USUBJID") %>%
  left_join(first_pd, by = "USUBJID") %>%
  arrange(dur_mo) %>%
  mutate(y = row_number())

RESP_COLORS <- c(CR = "#08306b", PR = "#6baed6", SD = "#e8a020", PD = "#c0392b", NE = "#999999")

p <- ggplot(trt) +
  geom_segment(aes(x = 0, xend = dur_mo, y = y, yend = y, color = BESTRSPC),
               linewidth = 3, lineend = "butt") +
  geom_point(data = filter(trt, !is.na(resp_mo)),
             aes(x = resp_mo, y = y), shape = 8, size = 2.2, color = "#d4af00", stroke = 1) +
  geom_point(data = filter(trt, !is.na(pd_mo)),
             aes(x = pd_mo, y = y), shape = 4, size = 2.2, color = "black", stroke = 1) +
  geom_text(data = filter(trt, ongoing),
            aes(x = dur_mo, y = y, label = "\u25b6"), size = 2.6, color = "black", hjust = -0.1) +
  scale_color_manual(values = RESP_COLORS, name = "Best overall response") +
  scale_x_continuous(name = "Months on treatment", expand = expansion(mult = c(0, 0.08))) +
  scale_y_continuous(name = NULL, breaks = NULL) +
  theme_classic(base_size = 10) +
  theme(panel.background = element_rect(fill = "white", color = NA)) +
  labs(
    title = "Treatment Exposure Plot (Swimmer Plot) \u2014 Vizatinib (ONCVIZ-001)",
    subtitle = sprintf("Treatment arm, n = %d \u00b7 sorted by ascending exposure duration \u00b7 \u2605 first CR/PR \u00b7 \u2715 first PD \u00b7 \u25b6 ongoing", nrow(trt)),
    caption = "Duration = TRTDURD / 30.44. Response markers use earliest ADRS assessment per category (no independent confirmation-window check)."
  )

ggsave(file.path(OUTPUT_DIR, "treatment_exposure_plot.png"), p,
       width = 9, height = max(4, 0.14 * nrow(trt)), dpi = 300, bg = "white")
message("Saved: treatment_exposure_plot.png")
