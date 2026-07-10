suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(lubridate)
})

DATA_DIR   <- "./Data/V1"
OUTPUT_DIR <- "./Out"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

adrand <- read.csv(file.path(DATA_DIR, "ADRAND.csv"), stringsAsFactors = FALSE)
adsl   <- read.csv(file.path(DATA_DIR, "ADSL.csv"),   stringsAsFactors = FALSE) %>% select(USUBJID, TRTSDT)

rand <- adrand %>% filter(RANDFL == "Y") %>%
  left_join(adsl, by = "USUBJID") %>%
  mutate(TRTSDT = as.Date(TRTSDT)) %>%
  arrange(TRTSDT) %>%
  mutate(cum = row_number())

monthly <- rand %>% mutate(month = floor_date(TRTSDT, "month")) %>% count(month)

scale_factor <- max(rand$cum) / max(monthly$n)

p <- ggplot() +
  geom_area(data = rand, aes(x = TRTSDT, y = cum), fill = "#4393C3", alpha = 0.25) +
  geom_step(data = rand, aes(x = TRTSDT, y = cum), color = "#08306B", linewidth = 1) +
  geom_col(data = monthly, aes(x = month, y = n * scale_factor), fill = "#B2182B", alpha = 0.35, width = 20) +
  scale_y_continuous(name = "Cumulative Patients Randomized",
                      sec.axis = sec_axis(~ . / scale_factor, name = "Patients enrolled per month")) +
  labs(x = "Calendar Date",
       title = "Patient Enrollment Over Time \u2014 ONCVIZ-001",
       caption = "Enrollment date proxied by first dose date (TRTSDT); randomized, ITT population. No planned-accrual curve available in source data.") +
  theme_classic(base_size = 12) +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 13),
        plot.caption = element_text(size = 8, face = "italic", hjust = 0.5),
        axis.title.y.left = element_text(color = "#08306B"),
        axis.title.y.right = element_text(color = "#B2182B"))

ggsave(file.path(OUTPUT_DIR, "enrollment_over_time.png"), p, width = 10, height = 6, dpi = 300)
