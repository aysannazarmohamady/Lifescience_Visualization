suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

# ── Paths ─────────────────────────────────────────────────────────────────────
DATA_DIR   <- "./Data/V1"
OUTPUT_DIR <- "./Out"
dir.create(OUTPUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ── Load data ─────────────────────────────────────────────────────────────────
adae <- read.csv(file.path(DATA_DIR, "ADAE.csv"), stringsAsFactors = FALSE)
adsl <- read.csv(file.path(DATA_DIR, "ADSL.csv"), stringsAsFactors = FALSE)

n_arm <- adsl %>% count(ARM) %>% { setNames(.$n, .$ARM) }

# ── Palette ───────────────────────────────────────────────────────────────────
C_TRT  <- "#2166AC"
C_CTRL <- "#B2182B"

# ── Any-grade incidence, patient-level, threshold >=5% in either arm ──────────
inc <- adae %>%
  distinct(AEPT, ARM, USUBJID) %>%
  count(AEPT, ARM, name = "n_pt") %>%
  complete(AEPT, ARM, fill = list(n_pt = 0)) %>%
  mutate(pct = 100 * n_pt / n_arm[ARM])

keep_pt <- inc %>%
  group_by(AEPT) %>%
  summarise(keep = any(pct >= 5)) %>%
  filter(keep) %>%
  pull(AEPT)

ord <- inc %>% filter(ARM == "TREATMENT", AEPT %in% keep_pt) %>% arrange(pct) %>% pull(AEPT)
inc <- inc %>% filter(AEPT %in% keep_pt) %>% mutate(AEPT = factor(AEPT, levels = ord))

# ── Grade >=3 incidence, same PT ordering ──────────────────────────────────────
inc_g3 <- adae %>%
  filter(AETOXGR >= 3) %>%
  distinct(AEPT, ARM, USUBJID) %>%
  count(AEPT, ARM, name = "n_pt") %>%
  complete(AEPT = keep_pt, ARM = names(n_arm), fill = list(n_pt = 0)) %>%
  mutate(pct = 100 * n_pt / n_arm[ARM],
         AEPT = factor(AEPT, levels = ord))

ctcae_v <- unique(na.omit(adae$CTCAE_V))[1]

arm_labels <- c(TREATMENT = paste0("Treatment (N=", n_arm["TREATMENT"], ")"),
                 CONTROL   = paste0("Control (N=", n_arm["CONTROL"], ")"))

p1 <- ggplot(inc, aes(x = AEPT, y = pct, fill = ARM)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_manual(values = c(TREATMENT = C_TRT, CONTROL = C_CTRL), labels = arm_labels) +
  coord_flip() +
  scale_x_discrete(limits = rev(ord)) +
  scale_y_reverse() +
  labs(x = NULL, y = "Patients with AE, any grade (%)", fill = NULL,
       title = "Any-Grade Adverse Events (\u2265 5% in either arm)") +
  theme_classic(base_size = 11) +
  theme(plot.title = element_text(size = 10))

p2 <- ggplot(inc_g3, aes(x = AEPT, y = pct, fill = ARM)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.65) +
  scale_fill_manual(values = c(TREATMENT = C_TRT, CONTROL = C_CTRL), labels = arm_labels) +
  coord_flip() +
  scale_x_discrete(limits = rev(ord)) +
  labs(x = NULL, y = "Patients with AE, Grade \u2265 3 (%)", fill = NULL,
       title = "Grade \u2265 3 Adverse Events (CTCAE)") +
  theme_classic(base_size = 11) +
  theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(),
        legend.position = "none", plot.title = element_text(size = 10))

combined <- (p1 | p2) +
  plot_layout(guides = "collect") +
  plot_annotation(
    title = "Adverse Events by Preferred Term and Treatment Arm",
    caption = sprintf("Safety population. Sorted by incidence in Treatment arm. CTCAE v%s; MedDRA-style Preferred Term (AEPT). Threshold: \u2265 5%% incidence in any arm.", ctcae_v),
    theme = theme(plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
                  plot.caption = element_text(size = 8, hjust = 0.5, face = "italic"))
  ) &
  theme(legend.position = "top")

ggsave(file.path(OUTPUT_DIR, "ae_bar_chart.png"), combined, width = 12, height = 7, dpi = 300)
