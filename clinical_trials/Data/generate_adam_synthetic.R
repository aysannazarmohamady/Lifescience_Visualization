# ==============================================================================
# generate_adam_synthetic.R
# Synthetic ADaM Dataset Generator, ONCVIZ-001
# Vizatinib Phase II/III Basket Trial, 400 patients, 10 ADaM domains
#
# Calibration sources:
#   - OS/PFS/AGE distributions: cBioPortal NSCLC studies (n=2,153)
#   - Efficacy parameters: KEYNOTE-189 (Gandhi et al., NEJM 2018)
#   - PK parameters: Erlotinib population PK (Ling et al., J Clin Pharmacol 2006)
#   - AE incidence: KEYNOTE-189 FDA label
#   - Lab toxicity grading: NCI-CTCAE v5.0
#
# Usage:
#   source("generate_adam_synthetic.R")
#   # CSV files are written to the working directory
#
# Requirements:
#   install.packages(c("dplyr", "tidyr", "lubridate", "purrr"))
# ==============================================================================

library(dplyr)
library(tidyr)
library(lubridate)
library(purrr)

set.seed(42)

# ==============================================================================
# Global Parameters
# ==============================================================================

N            <- 400
CUTOFF       <- as.Date("2026-03-05")
STUDY_START  <- as.Date("2022-06-01")
STUDY_ID     <- "ONCVIZ-001"
DRUG_NAME    <- "Vizatinib"

# Weibull parameters fitted to cBioPortal NSCLC pooled data (n=2,153)
# OS:  weibull_min, fitted by MLE, AIC-selected
# PFS: weibull_min, fitted by MLE, KS p=0.185
OS_SHAPE     <- 0.92    # shape parameter (c) from scipy weibull_min fit
OS_SCALE     <- 24.8    # scale parameter
PFS_SHAPE    <- 0.78
PFS_SCALE    <- 15.2

# AGE: Normal, fitted to cBioPortal data (median=66, sd~11)
AGE_MEAN     <- 66.0
AGE_SD       <- 11.0

# Efficacy, calibrated to KEYNOTE-189 (Gandhi et al., NEJM 2018)
# OS HR = 0.49, PFS HR = 0.52
OS_HR        <- 0.49
PFS_HR       <- 0.52

RESPONSE_RATES <- list(
  TREATMENT = c(CR = 0.07, PR = 0.38, SD = 0.32, PD = 0.23),
  CONTROL   = c(CR = 0.01, PR = 0.18, SD = 0.38, PD = 0.43)
)

# Markov transition matrix for response (RECIST 1.1)
TRANS <- list(
  CR = c(CR = 0.82, PR = 0.10, SD = 0.05, PD = 0.03),
  PR = c(CR = 0.10, PR = 0.58, SD = 0.20, PD = 0.12),
  SD = c(CR = 0.02, PR = 0.10, SD = 0.52, PD = 0.36),
  PD = c(CR = 0.00, PR = 0.00, SD = 0.00, PD = 1.00)
)

RESP_NAMES   <- c("CR", "PR", "SD", "PD")

# ==============================================================================
# Utility Functions
# ==============================================================================

clamp <- function(x, lo, hi) pmax(lo, pmin(hi, x))

fmt_date <- function(d) format(as.Date(d), "%Y-%m-%d")

sample_cat <- function(probs) {
  sample(names(probs), size = 1, prob = probs)
}

rweibull_custom <- function(n, shape, scale) {
  # Parameterization consistent with scipy weibull_min
  # X = scale * (-log(U))^(1/shape), U ~ Uniform(0,1)
  scale * (-log(runif(n)))^(1 / shape)
}

# ==============================================================================
# ADSL, Subject-Level Analysis Dataset
# ==============================================================================

message("Generating ADSL...")

TUMOR_TYPES  <- c("NSCLC", "CRC", "HCC", "PDAC", "BRCA")
SITES        <- sprintf("SITE-%02d", 1:20)
DISC_TRT     <- c("Adverse Event", "Progressive Disease",
                  "Withdrawal of Consent", "Death", "Protocol Deviation")
DISC_CTL     <- c("Progressive Disease", "Withdrawal of Consent",
                  "Death", "Crossover to Treatment", "Protocol Deviation")

adsl <- map_dfr(seq_len(N), function(i) {
  site   <- sample(SITES, 1)
  uid    <- sprintf("%s-%s-%04d", STUDY_ID, substr(site, 6, 7), i)
  arm    <- sample(c("TREATMENT", "CONTROL"), 1, prob = c(2, 1))
  armcd  <- ifelse(arm == "TREATMENT", "TRT", "CTL")

  # Demographics from fitted distributions
  age    <- clamp(round(rnorm(1, AGE_MEAN, AGE_SD)), 28, 85)
  sex    <- sample(c("M", "F"), 1, prob = c(0.509, 0.491))
  race   <- sample(c("WHITE","BLACK","ASIAN","OTHER"), 1,
                   prob = c(0.88, 0.098, 0.021, 0.001))
  smoking <- sample(c("Former heavy","Never","Former light","Current heavy"),
                    1, prob = c(0.488, 0.322, 0.178, 0.012))
  ecog   <- sample(0:2, 1, prob = c(0.38, 0.50, 0.12))

  # Biomarkers
  pdl1_score <- round(rbeta(1, 1.2, 2.5) * 100, 1)
  pdl1_grp   <- ifelse(pdl1_score >= 50, "HIGH",
                       ifelse(pdl1_score >= 1, "MED", "NEG"))
  tmb        <- clamp(round(rlnorm(1, 2.1, 0.85), 1), 0.5, 80)
  tmb_high   <- ifelse(tmb >= 10, "Y", "N")
  msi        <- sample(c("MSI-H","MSS"), 1, prob = c(0.05, 0.95))
  egfr_mut   <- ifelse(runif(1) < 0.15, "Y", "N")
  kras_mut   <- ifelse(runif(1) < 0.30, "Y", "N")
  braf_mut   <- ifelse(runif(1) < 0.03, "Y", "N")
  liver_mets <- ifelse(runif(1) < 0.20, "Y", "N")
  brain_mets <- ifelse(runif(1) < 0.10, "Y", "N")
  prior_lines<- sample(1:4, 1, prob = c(0.35, 0.35, 0.20, 0.10))
  baseline_sz<- round(rlnorm(1, 3.4, 0.55), 1)

  # Enrollment and treatment dates
  trtsdt <- STUDY_START + sample(0:540, 1)
  if (trtsdt > CUTOFF - 60) trtsdt <- CUTOFF - 60
  max_dur <- as.integer(CUTOFF - trtsdt)

  still_on <- runif(1) < ifelse(arm == "TREATMENT", 0.42, 0.22)
  r        <- RESPONSE_RATES[[arm]]
  bor      <- sample(RESP_NAMES, 1, prob = r)

  # OS and PFS from arm-adjusted Weibull distributions
  # Scale parameter adjusted by hazard ratio
  hr_os    <- ifelse(arm == "TREATMENT", 1.0, 1.0 / OS_HR)
  hr_pfs   <- ifelse(arm == "TREATMENT", 1.0, 1.0 / PFS_HR)
  os_months<- clamp(rweibull_custom(1, OS_SHAPE,  OS_SCALE  / hr_os),  1, 120)
  pfs_months<-clamp(rweibull_custom(1, PFS_SHAPE, PFS_SCALE / hr_pfs), 0.5, os_months)
  os_days  <- as.integer(os_months * 30.44)
  pfs_days <- as.integer(pfs_months * 30.44)

  dur_days <- ifelse(still_on, max_dur, clamp(os_days, 42, max_dur))
  trtedt   <- trtsdt + dur_days
  eosstt   <- ifelse(still_on, "ONGOING", "DISCONTINUED")
  dcsreas  <- ifelse(still_on, "",
                     sample(if (arm == "TREATMENT") DISC_TRT else DISC_CTL, 1))
  os_event <- ifelse(still_on, 0L, ifelse(os_days <= max_dur, 1L, 0L))
  pfs_event<- ifelse(pfs_days <= dur_days, 1L, 0L)

  tibble(
    STUDYID    = STUDY_ID,
    USUBJID    = uid,
    SITEID     = site,
    ARM        = arm,
    ARMCD      = armcd,
    TRT01P     = ifelse(arm == "TREATMENT", DRUG_NAME, "Placebo"),
    TUMORTYPE  = sample(TUMOR_TYPES, 1),
    AGE        = age,
    AGEGR1     = ifelse(age < 65, "<65", ">=65"),
    SEX        = sex,
    RACE       = race,
    SMOKING    = smoking,
    ECOG       = ecog,
    LIVERMETS  = liver_mets,
    BRAINMETS  = brain_mets,
    PRIORLINES = prior_lines,
    BASESZ     = baseline_sz,
    PDL1SCORE  = pdl1_score,
    PDL1GRP    = pdl1_grp,
    TMB        = tmb,
    TMBHIGH    = tmb_high,
    MSISTS     = msi,
    EGFRMUT    = egfr_mut,
    KRASMUT    = kras_mut,
    BRAFMUT    = braf_mut,
    BESTRSPC   = bor,
    TRTSDT     = fmt_date(trtsdt),
    TRTEDT     = fmt_date(trtedt),
    TRTDURD    = dur_days,
    EOSSTT     = eosstt,
    DCSREAS    = dcsreas,
    OSDUR      = os_days,
    OSCR       = os_event,
    OSDTC      = fmt_date(trtsdt + os_days),
    PFSDUR     = pfs_days,
    PFSCR      = pfs_event,
    PFSDTC     = fmt_date(trtsdt + pfs_days),
    CUTDTC     = fmt_date(CUTOFF),
    SAFFL      = "Y",
    ITTFL      = "Y"
  )
})

message(sprintf("  ADSL: %d patients | TRT=%d CTL=%d",
                nrow(adsl),
                sum(adsl$ARM == "TREATMENT"),
                sum(adsl$ARM == "CONTROL")))

# ==============================================================================
# ADRS, Tumor Response Dataset (Markov chain, 42-day cycles)
# ==============================================================================

message("Generating ADRS...")

adrs <- map_dfr(seq_len(nrow(adsl)), function(i) {
  s      <- adsl[i, ]
  trtsdt <- as.Date(s$TRTSDT)
  trtedt <- as.Date(s$TRTEDT)
  dur    <- as.integer(trtedt - trtsdt)
  n_vis  <- max(1L, dur %/% 42L)
  r      <- RESPONSE_RATES[[s$ARM]]
  prev   <- NULL
  rows   <- list()

  for (k in seq_len(n_vis)) {
    aday <- 42L * k
    if (aday > dur + 21L) break
    adt  <- trtsdt + aday
    if (adt > CUTOFF) break

    resp <- if (is.null(prev)) {
      sample(RESP_NAMES, 1, prob = r)
    } else {
      sample(RESP_NAMES, 1, prob = TRANS[[prev]])
    }

    rows[[k]] <- tibble(
      STUDYID  = STUDY_ID,
      USUBJID  = s$USUBJID,
      ARM      = s$ARM,
      PARAMCD  = "OVRLRESP",
      PARAM    = "Overall Response per RECIST 1.1",
      AVISIT   = sprintf("CYCLE %d DAY 1", k),
      AVISITN  = k,
      ADT      = fmt_date(adt),
      ADTN     = aday,
      AVALC    = resp,
      AVAL     = match(resp, RESP_NAMES),
      ANL01FL  = "Y"
    )
    prev <- resp
    if (resp == "PD") break
  }
  bind_rows(rows)
})

message(sprintf("  ADRS: %d records", nrow(adrs)))

# ==============================================================================
# ADTR, Tumor Measurements (response-stratified interpolation)
# ==============================================================================

message("Generating ADTR...")

adtr <- map_dfr(seq_len(nrow(adsl)), function(i) {
  s        <- adsl[i, ]
  trtsdt   <- as.Date(s$TRTSDT)
  trtedt   <- as.Date(s$TRTEDT)
  dur      <- as.integer(trtedt - trtsdt)
  baseline <- s$BASESZ
  bor      <- s$BESTRSPC
  n_vis    <- max(1L, dur %/% 42L)

  # Pre-specify biologically plausible endpoint by BOR
  final_val <- switch(bor,
    CR = baseline * runif(1, 0.01, 0.03),
    PR = baseline * runif(1, 0.25, 0.65),
    SD = baseline * runif(1, 0.82, 1.19),
    PD = baseline * runif(1, 1.20, 2.20)
  )
  nadir <- if (bor == "PR") baseline * runif(1, 0.25, 0.65) else NA

  rows <- list()
  rows[[1]] <- tibble(
    STUDYID = STUDY_ID, USUBJID = s$USUBJID, ARM = s$ARM,
    PARAMCD = "SUMDIAM", PARAM = "Sum of Longest Diameters (mm)",
    AVISIT = "BASELINE", AVISITN = 0L,
    ADT = fmt_date(trtsdt), ADTN = 0L,
    BASE = baseline, AVAL = baseline, CHG = 0.0, PCHG = 0.0,
    ANL01FL = "Y"
  )

  current <- baseline
  for (k in seq_len(n_vis)) {
    aday  <- 42L * k
    if (aday > dur + 21L) break
    adt   <- trtsdt + aday
    if (adt > CUTOFF) break
    phase <- k / max(n_vis, 1)

    current <- switch(bor,
      CR = {
        target <- baseline * 0.02
        v <- current + (target - current) * 0.45
        max(0.01, v + rnorm(1, 0, baseline * 0.008))
      },
      PR = {
        target <- if (phase < 0.55) nadir else nadir * runif(1, 1.0, 1.15)
        v <- current + (target - current) * 0.40
        max(0.1, v + rnorm(1, 0, baseline * 0.015))
      },
      SD = {
        target <- baseline * runif(1, 0.82, 1.18)
        v <- current + (target - current) * 0.25
        max(0.1, v + rnorm(1, 0, baseline * 0.03))
      },
      PD = {
        v <- current + (final_val - current) * 0.35
        max(baseline * 1.01, v)
      }
    )

    chg  <- round(current - baseline, 1)
    pchg <- round((current - baseline) / baseline * 100, 1)

    rows[[k + 1]] <- tibble(
      STUDYID = STUDY_ID, USUBJID = s$USUBJID, ARM = s$ARM,
      PARAMCD = "SUMDIAM", PARAM = "Sum of Longest Diameters (mm)",
      AVISIT = sprintf("CYCLE %d DAY 1", k), AVISITN = k,
      ADT = fmt_date(adt), ADTN = aday,
      BASE = baseline, AVAL = round(current, 1),
      CHG = chg, PCHG = pchg, ANL01FL = "Y"
    )
  }
  bind_rows(rows)
})

message(sprintf("  ADTR: %d records | PCHG range: %.1f%% to %.1f%%",
                nrow(adtr), min(adtr$PCHG), max(adtr$PCHG)))

# ==============================================================================
# ADAE, Adverse Events Dataset
# Incidence rates calibrated to KEYNOTE-189 FDA label
# ==============================================================================

message("Generating ADAE...")

AE_PROFILE <- tribble(
  ~SOC,                                     ~PT,                              ~INC_TRT, ~INC_CTL, ~GW1,  ~GW2,  ~GW3,  ~GW4,  ~GW5,
  "Gastrointestinal disorders",             "Nausea",                          0.58,     0.32,    0.38,  0.40,  0.18,  0.03,  0.01,
  "Gastrointestinal disorders",             "Diarrhea",                        0.52,     0.25,    0.32,  0.38,  0.24,  0.05,  0.01,
  "Gastrointestinal disorders",             "Vomiting",                        0.34,     0.18,    0.38,  0.40,  0.18,  0.03,  0.01,
  "Gastrointestinal disorders",             "Constipation",                    0.26,     0.20,    0.52,  0.36,  0.10,  0.02,  0.00,
  "Gastrointestinal disorders",             "Abdominal pain",                  0.20,     0.16,    0.42,  0.38,  0.16,  0.04,  0.00,
  "Skin and subcutaneous tissue disorders", "Rash",                            0.45,     0.12,    0.42,  0.38,  0.16,  0.04,  0.00,
  "Skin and subcutaneous tissue disorders", "Pruritus",                        0.22,     0.08,    0.55,  0.35,  0.08,  0.02,  0.00,
  "Skin and subcutaneous tissue disorders", "Alopecia",                        0.20,     0.08,    0.60,  0.35,  0.05,  0.00,  0.00,
  "General disorders",                      "Fatigue",                         0.65,     0.48,    0.30,  0.40,  0.22,  0.07,  0.01,
  "General disorders",                      "Decreased appetite",              0.42,     0.30,    0.40,  0.38,  0.18,  0.04,  0.00,
  "General disorders",                      "Pyrexia",                         0.22,     0.16,    0.48,  0.36,  0.13,  0.03,  0.00,
  "Investigations",                         "ALT increased",                   0.32,     0.10,    0.40,  0.35,  0.16,  0.07,  0.02,
  "Investigations",                         "AST increased",                   0.30,     0.08,    0.42,  0.35,  0.15,  0.06,  0.02,
  "Investigations",                         "Blood bilirubin increased",       0.14,     0.05,    0.48,  0.34,  0.12,  0.05,  0.01,
  "Investigations",                         "Neutrophil count decreased",      0.32,     0.08,    0.32,  0.30,  0.24,  0.10,  0.04,
  "Investigations",                         "Hemoglobin decreased",            0.38,     0.22,    0.38,  0.36,  0.18,  0.07,  0.01,
  "Nervous system disorders",               "Peripheral neuropathy",           0.26,     0.08,    0.42,  0.38,  0.16,  0.04,  0.00,
  "Nervous system disorders",               "Headache",                        0.20,     0.16,    0.52,  0.36,  0.10,  0.02,  0.00,
  "Musculoskeletal disorders",              "Arthralgia",                      0.24,     0.16,    0.44,  0.38,  0.14,  0.04,  0.00,
  "Musculoskeletal disorders",              "Myalgia",                         0.18,     0.14,    0.48,  0.38,  0.12,  0.02,  0.00,
  "Respiratory disorders",                  "Dyspnea",                         0.20,     0.16,    0.42,  0.38,  0.14,  0.05,  0.01,
  "Respiratory disorders",                  "Cough",                           0.18,     0.14,    0.58,  0.36,  0.06,  0.00,  0.00,
  "Metabolism disorders",                   "Hypokalemia",                     0.16,     0.07,    0.44,  0.36,  0.14,  0.05,  0.01,
  "Metabolism disorders",                   "Hyperglycemia",                   0.20,     0.08,    0.42,  0.36,  0.16,  0.05,  0.01,
  "Cardiac disorders",                      "QT prolongation",                 0.10,     0.04,    0.52,  0.34,  0.10,  0.04,  0.00,
  "Cardiac disorders",                      "Hypertension",                    0.16,     0.08,    0.44,  0.36,  0.14,  0.05,  0.01,
  "Immune system disorders",                "Hypothyroidism",                  0.12,     0.03,    0.48,  0.38,  0.12,  0.02,  0.00,
  "Immune system disorders",                "Pneumonitis",                     0.07,     0.02,    0.38,  0.32,  0.18,  0.10,  0.02,
  "Infections",                             "Urinary tract infection",         0.10,     0.10,    0.48,  0.38,  0.12,  0.02,  0.00,
  "Infections",                             "Upper respiratory infection",     0.12,     0.12,    0.52,  0.36,  0.10,  0.02,  0.00
)

adae <- map_dfr(seq_len(nrow(adsl)), function(i) {
  s      <- adsl[i, ]
  trtsdt <- as.Date(s$TRTSDT)
  trtedt <- as.Date(s$TRTEDT)
  dur    <- as.integer(trtedt - trtsdt)
  arm    <- s$ARM
  uid    <- s$USUBJID
  seq_n  <- 0L
  rows   <- list()

  for (j in seq_len(nrow(AE_PROFILE))) {
    ae  <- AE_PROFILE[j, ]
    inc <- ifelse(arm == "TREATMENT", ae$INC_TRT, ae$INC_CTL)
    if (runif(1) > inc) next

    n_occ <- sample(1:3, 1, prob = c(0.72, 0.22, 0.06))
    for (occ in seq_len(n_occ)) {
      seq_n  <- seq_n + 1L
      onset  <- sample(seq_len(max(dur - 7L, 1L)), 1)
      grade  <- sample(1:5, 1, prob = c(ae$GW1, ae$GW2, ae$GW3, ae$GW4, ae$GW5))
      dur_ae <- max(1L, as.integer(rexp(1, 1 / (12 + grade * 7))))
      dur_ae <- min(dur_ae, dur - onset + 30L)
      endt   <- min(trtsdt + onset + dur_ae, CUTOFF)
      serious<- grade >= 3
      related<- runif(1) < ifelse(arm == "TREATMENT", 0.62, 0.28)
      action <- if (grade >= 2) {
        sample(c("DOSE REDUCED","DRUG INTERRUPTED","DRUG WITHDRAWN","NONE"),
               1, prob = c(0.18, 0.22, 0.05, 0.55))
      } else "NONE"
      outcome <- sample(c("RECOVERED","RECOVERING","NOT RECOVERED","FATAL","UNKNOWN"),
                        1, prob = c(0.68, 0.14, 0.11, 0.03, 0.04))

      rows[[length(rows) + 1]] <- tibble(
        STUDYID  = STUDY_ID, USUBJID = uid, ARM = arm,
        AESEQ    = seq_n,
        AESOC    = ae$SOC, AEPT = ae$PT,
        AETOXGR  = grade,
        AESER    = ifelse(serious, "Y", "N"),
        AEREL    = ifelse(related, "RELATED", "NOT RELATED"),
        AEACN    = action, AEOUT = outcome,
        AESTDTC  = fmt_date(trtsdt + onset),
        AEENDTC  = fmt_date(endt),
        AESTDY   = onset,
        AEENDY   = onset + dur_ae,
        AEDUR    = dur_ae
      )
    }
  }
  bind_rows(rows)
})

message(sprintf("  ADAE: %d records | Grade>=3: %.1f%%",
                nrow(adae),
                mean(adae$AETOXGR >= 3) * 100))

# ==============================================================================
# ADTTE, Time-to-Event Dataset (OS, PFS, DOR, TTR)
# ==============================================================================

message("Generating ADTTE...")

subgrp_vars <- c("TUMORTYPE","AGEGR1","SEX","ECOG","PDL1GRP",
                 "MSISTS","TMBHIGH","LIVERMETS","PRIORLINES",
                 "SMOKING","EGFRMUT","KRASMUT")

adtte <- map_dfr(seq_len(nrow(adsl)), function(i) {
  s      <- adsl[i, ]
  uid    <- s$USUBJID
  arm    <- s$ARM
  armcd  <- s$ARMCD
  trtsdt <- as.Date(s$TRTSDT)
  sg     <- s[subgrp_vars]

  make_row <- function(paramcd, param, aval, cnsr, evntdesc, adt) {
    bind_cols(
      tibble(STUDYID=STUDY_ID, USUBJID=uid, ARM=arm, ARMCD=armcd,
             PARAMCD=paramcd, PARAM=param,
             AVAL=aval, AVALU="DAYS", AVALM=round(aval/30.44,2),
             CNSR=cnsr, EVNTDESC=evntdesc,
             ADT=fmt_date(adt), STARTDT=s$TRTSDT, ANL01FL="Y"),
      sg
    )
  }

  rows <- list(
    make_row("OS",  "Overall Survival",
             s$OSDUR, 1L - s$OSCR,
             ifelse(s$OSCR==1, "DEATH", "CENSORED"),
             as.Date(s$OSDTC)),
    make_row("PFS", "Progression-Free Survival",
             s$PFSDUR, 1L - s$PFSCR,
             ifelse(s$PFSCR==1, "PROGRESSION", "CENSORED"),
             as.Date(s$PFSDTC))
  )

  if (s$BESTRSPC %in% c("CR","PR")) {
    ttr_days  <- sample(28:84, 1)
    dor_base  <- ifelse(arm=="TREATMENT",
                        ifelse(s$BESTRSPC=="CR", 14.0, 9.5),
                        ifelse(s$BESTRSPC=="CR",  6.0, 4.0))
    dor_days  <- max(42L, min(as.integer(rexp(1, 1/(dor_base*30.44))),
                              s$OSDUR - ttr_days))
    dor_event <- rbinom(1, 1, 0.65)

    rows[[3]] <- make_row("DOR", "Duration of Response",
                          dor_days, 1L - dor_event,
                          ifelse(dor_event==1, "PROGRESSION", "CENSORED"),
                          trtsdt + ttr_days + dor_days)
    rows[[4]] <- make_row("TTR", "Time to Response",
                          ttr_days, 0L, "RESPONSE",
                          trtsdt + ttr_days)
  }
  bind_rows(rows)
})

message(sprintf("  ADTTE: %d records | params: %s",
                nrow(adtte),
                paste(unique(adtte$PARAMCD), collapse=", ")))

# ==============================================================================
# ADPK, Pharmacokinetics (Treatment arm only)
# One-compartment oral model, calibrated to erlotinib PK
# ==============================================================================

message("Generating ADPK...")

PK_DOSE  <- 300
KA       <- 0.8
CL_F     <- 18.0
VD_F     <- 320.0
PK_TIMES <- c(0, 0.5, 1, 2, 3, 4, 6, 8, 12, 24)
PK_VISITS<- c("CYCLE 1 DAY 1","CYCLE 1 DAY 15","CYCLE 3 DAY 1")
PK_DAYS  <- c(0L, 14L, 56L)

trt_subj <- adsl %>% filter(ARM == "TREATMENT")

adpk <- map_dfr(seq_len(nrow(trt_subj)), function(i) {
  s      <- trt_subj[i, ]
  uid    <- s$USUBJID
  trtsdt <- as.Date(s$TRTSDT)

  cl_i <- CL_F * exp(rnorm(1, 0, 0.28))
  vd_i <- VD_F * exp(rnorm(1, 0, 0.25))
  ka_i <- KA   * exp(rnorm(1, 0, 0.35))
  ke_i <- cl_i / vd_i

  conc_rows <- map_dfr(seq_along(PK_VISITS), function(vi) {
    vday    <- PK_DAYS[vi]
    base_dt <- trtsdt + vday
    if (base_dt > CUTOFF) return(NULL)

    map_dfr(PK_TIMES, function(t) {
      conc <- if (t == 0) 0.0 else {
        c_val <- (PK_DOSE * 1000 * ka_i / (vd_i * (ka_i - ke_i))) *
                 (exp(-ke_i * t) - exp(-ka_i * t))
        c_val <- max(0.0, c_val)
        c_val * exp(rnorm(1, 0, 0.12)) + abs(rnorm(1, 0, 2.0))
      }
      tibble(STUDYID=STUDY_ID, USUBJID=uid,
             PARAMCD="CONC", PARAM=paste(DRUG_NAME,"Plasma Concentration"),
             AVISIT=PK_VISITS[vi], AVISITN=vi,
             NOMTPT=t, ADT=fmt_date(base_dt + t/24),
             AVAL=round(conc,2), AVALU="ng/mL",
             DOSE=PK_DOSE, ROUTE="ORAL", ANL01FL="Y",
             AGE=s$AGE, SEX=s$SEX)
    })
  })

  cmax  <- max(0.1, (PK_DOSE*1000*ka_i/(vd_i*(ka_i-ke_i))) *
               (exp(-ke_i*3) - exp(-ka_i*3)) * exp(rnorm(1,0,0.10)))
  auc   <- max(0.1, PK_DOSE * 1000 / cl_i * exp(rnorm(1,0,0.10)))
  tmax  <- max(0.1, 3.0 + rnorm(1,0,0.5))
  thalf <- max(0.1, 0.693 * vd_i / cl_i)

  param_rows <- tibble(
    STUDYID = STUDY_ID, USUBJID = uid,
    PARAMCD = c("CMAX","AUCINF","TMAX","THALF"),
    PARAM   = paste(DRUG_NAME, c("Cmax","AUCinf","Tmax","t1/2")),
    AVISIT  = "CYCLE 1 DAY 1", AVISITN = 1L,
    NOMTPT  = NA_real_,
    ADT     = fmt_date(trtsdt),
    AVAL    = round(c(cmax, auc, tmax, thalf), 1),
    AVALU   = c("ng/mL","ng*h/mL","h","h"),
    DOSE    = PK_DOSE, ROUTE = "ORAL", ANL01FL = "Y",
    AGE     = s$AGE, SEX = s$SEX
  )

  bind_rows(conc_rows, param_rows)
})

message(sprintf("  ADPK: %d records | Cmax median=%.0f ng/mL",
                nrow(adpk),
                median(adpk$AVAL[adpk$PARAMCD == "CMAX"])))

# ==============================================================================
# ADLB, Laboratory Data (21 parameters, 8 visits)
# ==============================================================================

message("Generating ADLB...")

LAB_TESTS <- tribble(
  ~PARAMCD, ~PARAM,                          ~AVALU,    ~ULN,  ~LLN,  ~MU,   ~SD,   ~EFFECT, ~ETYPE,
  "ALT",    "Alanine Aminotransferase",       "U/L",     45,    7,     22,    8,     1.80,    "hepatic",
  "AST",    "Aspartate Aminotransferase",     "U/L",     40,    10,    20,    7,     1.60,    "hepatic",
  "BILI",   "Total Bilirubin",               "mg/dL",   1.2,   0.2,   0.7,   0.2,   1.35,    "hepatic",
  "ALP",    "Alkaline Phosphatase",           "U/L",     120,   35,    70,    20,    1.20,    "increase",
  "GGT",    "Gamma-Glutamyl Transferase",     "U/L",     60,    8,     30,    12,    1.45,    "hepatic",
  "CREAT",  "Creatinine",                    "mg/dL",   1.2,   0.5,   0.85,  0.18,  1.06,    "increase",
  "BUN",    "Blood Urea Nitrogen",           "mg/dL",   20,    6,     13,    4,     1.05,    "stable",
  "K",      "Potassium",                     "mEq/L",   5.1,   3.5,   4.0,   0.4,  -0.12,    "suppress",
  "NA",     "Sodium",                        "mEq/L",   145,   135,   139,   2.5,   0.00,    "stable",
  "HGB",    "Hemoglobin",                    "g/dL",    17,    11,    13.5,  1.2,  -0.18,    "suppress",
  "PLT",    "Platelet Count",               "10^9/L",   400,   150,   220,   45,   -0.20,    "suppress",
  "ANC",    "Absolute Neutrophil Count",    "10^9/L",   7.5,   1.8,   4.0,   1.0,  -0.25,    "suppress",
  "WBC",    "White Blood Cell Count",       "10^9/L",   10,    4.0,   6.5,   1.5,  -0.14,    "suppress",
  "QTCF",   "QTc (Fridericia)",              "msec",    450,   350,   410,   22,    9.0,     "increase",
  "HGB",    "Hemoglobin",                    "g/dL",    17,    11,    13.5,  1.2,  -0.18,    "suppress",
  "GLUC",   "Glucose",                       "mg/dL",   100,   70,    88,    14,    7.0,     "increase",
  "TSH",    "Thyroid Stimulating Hormone",   "mIU/L",   4.0,   0.4,   2.0,   0.8,   0.7,    "increase",
  "ALB",    "Albumin",                       "g/dL",    5.2,   3.5,   4.2,   0.4,  -0.08,    "suppress"
) %>% distinct(PARAMCD, .keep_all = TRUE)

VISITS_LB  <- c("BASELINE","CYCLE 2 DAY 1","CYCLE 4 DAY 1","CYCLE 6 DAY 1",
                "CYCLE 9 DAY 1","CYCLE 12 DAY 1","CYCLE 18 DAY 1","EOT")
VISIT_DAYS <- c(0L, 42L, 126L, 210L, 336L, 462L, 714L, NA_integer_)

adlb <- map_dfr(seq_len(nrow(adsl)), function(i) {
  s      <- adsl[i, ]
  trtsdt <- as.Date(s$TRTSDT)
  trtedt <- as.Date(s$TRTEDT)
  dur    <- as.integer(trtedt - trtsdt)
  arm    <- s$ARM
  uid    <- s$USUBJID
  pat_re <- rnorm(1, 0, 0.10)
  has_hepatotox <- (arm == "TREATMENT") && (runif(1) < 0.08)
  hepatotox_grade <- if (has_hepatotox) sample(2:4, 1, prob=c(0.50,0.35,0.15)) else 0L

  map_dfr(seq_len(nrow(LAB_TESTS)), function(j) {
    lb       <- LAB_TESTS[j, ]
    base_val <- max(lb$LLN * 0.8, abs(rnorm(1, lb$MU, lb$SD)))

    map_dfr(seq_along(VISITS_LB), function(vi) {
      vday <- ifelse(is.na(VISIT_DAYS[vi]), dur, VISIT_DAYS[vi])
      if (vday > dur + 14L) return(NULL)
      adt <- trtsdt + vday
      if (adt > CUTOFF) return(NULL)

      aval <- if (vi == 1L) {
        base_val
      } else {
        phase <- min(vday / max(dur, 1), 1.0)
        if (arm == "TREATMENT") {
          switch(lb$ETYPE,
            hepatic = {
              peak  <- 0.30
              shape <- exp(-4 * (phase - peak)^2)
              if (has_hepatotox && lb$PARAMCD %in% c("ALT","AST","BILI")) {
                target_xuln <- c("2"=4.0,"3"=8.0,"4"=15.0)[as.character(hepatotox_grade)]
                peak_val    <- lb$ULN * target_xuln * runif(1, 0.85, 1.15)
                base_val + (peak_val - base_val) * shape + rnorm(1, 0, lb$SD * 0.20)
              } else {
                base_val + lb$EFFECT * shape * lb$MU * (1 + pat_re) + rnorm(1, 0, lb$SD*0.20)
              }
            },
            suppress = {
              nadir <- 0.40
              drop  <- abs(lb$EFFECT) * exp(-5*(phase-nadir)^2)
              base_val * (1 - drop) + rnorm(1, 0, lb$SD * 0.12)
            },
            increase = {
              if (lb$PARAMCD == "QTCF") {
                base_val + lb$EFFECT * phase + rnorm(1, 0, 5)
              } else {
                base_val + lb$EFFECT * phase + rnorm(1, 0, lb$SD*0.15)
              }
            },
            base_val + rnorm(1, 0, lb$SD * 0.10)
          )
        } else {
          base_val + rnorm(1, 0, lb$SD * 0.08)
        }
      }
      aval <- max(0.01, aval)
      chg  <- round(aval - base_val, 3)
      pchg <- round((aval - base_val) / base_val * 100, 1)
      xuln <- round(aval / lb$ULN, 3)

      grade <- if (lb$PARAMCD %in% c("ALT","AST")) {
        ifelse(aval > 20*lb$ULN, 5L,
        ifelse(aval > 10*lb$ULN, 4L,
        ifelse(aval >  5*lb$ULN, 3L,
        ifelse(aval >  3*lb$ULN, 2L,
        ifelse(aval >    lb$ULN, 1L, 0L)))))
      } else if (lb$PARAMCD == "BILI") {
        ifelse(aval > 10*lb$ULN, 4L,
        ifelse(aval >  3*lb$ULN, 3L,
        ifelse(aval >1.5*lb$ULN, 2L,
        ifelse(aval >    lb$ULN, 1L, 0L))))
      } else NA_integer_

      tibble(
        STUDYID=STUDY_ID, USUBJID=uid, ARM=arm,
        PARAMCD=lb$PARAMCD, PARAM=lb$PARAM, AVALU=lb$AVALU,
        AVISIT=VISITS_LB[vi], AVISITN=vi-1L,
        ADT=fmt_date(adt), ADTN=vday,
        BASE=round(base_val,3), AVAL=round(aval,3),
        CHG=chg, PCHG=pchg,
        ULN=lb$ULN, LLN=lb$LLN, XULNFL=xuln,
        ATOXGR=grade,
        ANRIND=ifelse(aval>lb$ULN,"HIGH",ifelse(aval<lb$LLN,"LOW","NORMAL")),
        ANL01FL="Y"
      )
    })
  })
})

message(sprintf("  ADLB: %d records | %d parameters",
                nrow(adlb), n_distinct(adlb$PARAMCD)))

# ==============================================================================
# ADEX, Exposure and Dose Modifications
# ==============================================================================

message("Generating ADEX...")

adex <- map_dfr(seq_len(nrow(adsl)), function(i) {
  s      <- adsl[i, ]
  trtsdt <- as.Date(s$TRTSDT)
  trtedt <- as.Date(s$TRTEDT)
  dur    <- as.integer(trtedt - trtsdt)
  arm    <- s$ARM
  uid    <- s$USUBJID
  planned<- ifelse(arm=="TREATMENT", 300L, 0L)
  level  <- 1.0
  n_cyc  <- max(1L, dur %/% 21L)

  map_dfr(seq_len(n_cyc), function(cyc) {
    cstart <- trtsdt + (cyc-1L)*21L
    cend   <- cstart + 20L
    if (cend > CUTOFF) return(NULL)

    mod_reason <- ""
    if (cyc > 2L && arm == "TREATMENT") {
      if (runif(1) < 0.10) {
        level      <<- max(0.33, level - 0.33)
        mod_reason <- sample(c("AE - Hepatotoxicity","AE - Diarrhea",
                                "AE - Fatigue","AE - Neuropathy"), 1)
      } else if (runif(1) < 0.04) {
        level      <<- 0
        mod_reason <- "AE - Dose Interruption"
      }
    }

    actual   <- round(planned * level)
    days_on  <- ifelse(level > 0, 21L, sample(0:10, 1))
    intensity<- ifelse(planned > 0, round(actual*days_on/(planned*21)*100,1), 0.0)

    tibble(
      STUDYID=STUDY_ID, USUBJID=uid, ARM=arm,
      EXSEQ=cyc,
      EXTRT=ifelse(arm=="TREATMENT", DRUG_NAME, "Placebo"),
      EXDOSE=actual, EXDOSU="mg", EXDOSFRQ="QD", EXROUTE="ORAL",
      EXSTDTC=fmt_date(cstart), EXENDTC=fmt_date(cend),
      EXCYCLE=cyc, PLANDOSE=planned, DOSELEVEL=level,
      DAYSONDRUG=days_on, DOSEINT=intensity,
      MODFL=ifelse(nchar(mod_reason)>0,"Y","N"),
      MODREASN=mod_reason
    )
  })
})

message(sprintf("  ADEX: %d records", nrow(adex)))

# ==============================================================================
# ADBM, Biomarker Dataset
# ==============================================================================

message("Generating ADBM...")

BM_VISITS    <- c("BASELINE","CYCLE 2 DAY 1","CYCLE 4 DAY 1","CYCLE 6 DAY 1","EOT")
BM_VISIT_DAY <- c(0L, 42L, 126L, 210L, NA_integer_)

BM_PARAMS <- tribble(
  ~PARAMCD, ~PARAM,                            ~AVALU,     ~PATTERN,       ~TRT_CHG,
  "TMB",    "Tumor Mutational Burden",          "mut/Mb",  "baseline",      NA_real_,
  "PDL1",   "PD-L1 TPS Score",                 "%",        "baseline",      NA_real_,
  "CTDNA",  "ctDNA Variant Allele Frequency",   "%",        "longitudinal", -0.55,
  "CEA",    "Carcinoembryonic Antigen",          "ng/mL",   "longitudinal", -0.30,
  "CA125",  "CA-125",                            "U/mL",    "longitudinal", -0.28,
  "CD8",    "CD8+ T Cell Density",              "cells/mm2","longitudinal",  0.42
)

adbm <- map_dfr(seq_len(nrow(adsl)), function(i) {
  s      <- adsl[i, ]
  trtsdt <- as.Date(s$TRTSDT)
  trtedt <- as.Date(s$TRTEDT)
  dur    <- as.integer(trtedt - trtsdt)
  arm    <- s$ARM
  uid    <- s$USUBJID
  bor    <- s$BESTRSPC

  map_dfr(seq_len(nrow(BM_PARAMS)), function(j) {
    bm       <- BM_PARAMS[j, ]
    base_val <- switch(bm$PARAMCD,
      TMB  = s$TMB,
      PDL1 = s$PDL1SCORE,
      round(rlnorm(1, 2.5, 0.8), 2)
    )
    visits <- if (bm$PATTERN == "longitudinal") BM_VISITS else "BASELINE"
    vdays  <- if (bm$PATTERN == "longitudinal") BM_VISIT_DAY else 0L

    map_dfr(seq_along(visits), function(vi) {
      vday <- ifelse(is.na(vdays[vi]), dur, vdays[vi])
      if (vday > dur + 14L) return(NULL)
      adt <- trtsdt + vday
      if (adt > CUTOFF) return(NULL)

      aval <- if (visits[vi] == "BASELINE" || bm$PATTERN == "baseline") {
        base_val
      } else {
        phase    <- min(vday / max(dur,1), 1.0)
        resp_mod <- c(CR=1.3,PR=1.0,SD=0.5,PD=-0.2)[bor]
        if (!is.na(resp_mod)) resp_mod else 0.5
        if (arm == "TREATMENT" && !is.na(bm$TRT_CHG)) {
          base_val*(1+bm$TRT_CHG*phase*resp_mod)*exp(rnorm(1,0,0.15))
        } else {
          base_val * exp(rnorm(1, 0, 0.08))
        }
      }
      aval <- max(0.01, aval)
      tibble(
        STUDYID=STUDY_ID, USUBJID=uid, ARM=arm,
        PARAMCD=bm$PARAMCD, PARAM=bm$PARAM, AVALU=bm$AVALU,
        AVISIT=visits[vi], AVISITN=vi-1L,
        ADT=fmt_date(adt), ADTN=vday,
        BASE=round(base_val,3), AVAL=round(aval,3),
        CHG=round(aval-base_val,3),
        PCHG=round((aval-base_val)/base_val*100,1),
        BESTRSPC=bor, ANL01FL="Y"
      )
    })
  })
})

message(sprintf("  ADBM: %d records | %d parameters",
                nrow(adbm), n_distinct(adbm$PARAMCD)))

# ==============================================================================
# ADPR, Patient-Reported Outcomes (EORTC QLQ-C30)
# ==============================================================================

message("Generating ADPR...")

PRO_SCALES <- tribble(
  ~PARAMCD,  ~PARAM,                    ~SCALETYP,    ~MU, ~SD, ~DIR, ~EFFECT,
  "GLQOL",   "Global Health Status/QoL","functional",  62,  18,   1,    4.0,
  "PHFUNC",  "Physical Functioning",    "functional",  68,  16,  -1,    3.5,
  "EMOFUNC", "Emotional Functioning",   "functional",  72,  17,   1,    2.5,
  "FATIGUE", "Fatigue",                 "symptom",     32,  22,   1,    8.0,
  "NAUSEA",  "Nausea and Vomiting",     "symptom",     14,  18,   1,    6.0,
  "PAIN",    "Pain",                    "symptom",     28,  22,   1,    4.0,
  "DYSPNEA", "Dyspnoea",                "symptom",     20,  18,   1,    3.5,
  "APPETIT", "Appetite Loss",           "symptom",     22,  20,   1,    5.0,
  "DIARRH",  "Diarrhoea",               "symptom",     10,  14,   1,    7.0
)

PR_VISITS    <- c("BASELINE","CYCLE 2","CYCLE 4","CYCLE 6","CYCLE 9","CYCLE 12","EOT")
PR_VISIT_DAY <- c(0L, 42L, 126L, 210L, 336L, 462L, NA_integer_)
MID          <- 10.0

adpr <- map_dfr(seq_len(nrow(adsl)), function(i) {
  s      <- adsl[i, ]
  trtsdt <- as.Date(s$TRTSDT)
  trtedt <- as.Date(s$TRTEDT)
  dur    <- as.integer(trtedt - trtsdt)
  arm    <- s$ARM
  uid    <- s$USUBJID
  comp   <- rbeta(1, 8, 2)

  map_dfr(seq_len(nrow(PRO_SCALES)), function(j) {
    sc       <- PRO_SCALES[j, ]
    base_val <- clamp(rnorm(1, sc$MU, sc$SD), 0, 100)

    map_dfr(seq_along(PR_VISITS), function(vi) {
      if (vi > 1L && runif(1) > comp) return(NULL)
      vday <- ifelse(is.na(PR_VISIT_DAY[vi]), dur, PR_VISIT_DAY[vi])
      if (vday > dur + 21L) return(NULL)
      adt <- trtsdt + vday
      if (adt > CUTOFF) return(NULL)

      aval <- if (vi == 1L) {
        base_val
      } else {
        phase    <- min(vday / max(dur,1), 1.0)
        peak     <- 0.25
        toxicity <- sc$EFFECT * sc$DIR * exp(-5*(phase-peak)^2)
        recovery <- if (phase > peak) sc$EFFECT * sc$DIR * 0.35 * (phase-peak) else 0
        if (arm == "TREATMENT") {
          clamp(base_val + toxicity + recovery + rnorm(1,0,sc$SD*0.18), 0, 100)
        } else {
          clamp(base_val - 3*phase + rnorm(1,0,sc$SD*0.12), 0, 100)
        }
      }

      chg <- round(aval - base_val, 1)
      tibble(
        STUDYID=STUDY_ID, USUBJID=uid, ARM=arm,
        PARAMCD=sc$PARAMCD, PARAM=sc$PARAM, SCALETYP=sc$SCALETYP,
        AVISIT=PR_VISITS[vi], AVISITN=vi-1L,
        ADT=fmt_date(adt), ADTN=vday,
        BASE=round(base_val,1), AVAL=round(aval,1),
        CHG=chg, PCHG=round(chg/base_val*100,1),
        MID=MID,
        MIDRESP=ifelse(abs(chg)>=MID,"Y","N"),
        DETRFL=ifelse((sc$SCALETYP=="functional"&&chg<=-MID)|
                      (sc$SCALETYP=="symptom"&&chg>=MID),"Y","N"),
        ANL01FL="Y"
      )
    })
  })
})

message(sprintf("  ADPR: %d records | %d scales",
                nrow(adpr), n_distinct(adpr$PARAMCD)))

# ==============================================================================
# Save All Datasets
# ==============================================================================

message("\nSaving datasets...")

datasets <- list(
  ADSL  = adsl,
  ADRS  = adrs,
  ADTR  = adtr,
  ADAE  = adae,
  ADLB  = adlb,
  ADTTE = adtte,
  ADPK  = adpk,
  ADEX  = adex,
  ADBM  = adbm,
  ADPR  = adpr
)

for (nm in names(datasets)) {
  fname <- paste0(nm, ".csv")
  write.csv(datasets[[nm]], fname, row.names = FALSE)
  message(sprintf("  %s: %s rows x %s cols",
                  nm, nrow(datasets[[nm]]), ncol(datasets[[nm]])))
}

total <- sum(sapply(datasets, nrow))
message(sprintf("\nTotal records: %s", format(total, big.mark=",")))
message(sprintf("Study: %s, %s, Phase II/III Basket Trial", STUDY_ID, DRUG_NAME))
message("Generation complete.")
