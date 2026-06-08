suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(purrr)
})

generate_adam_oncviz001 <- function(
    output_dir  = "./data_v2",
    seed        = 42L,
    n_phase1    = 20L,
    n_phase2_trt = 40L,
    n_phase2_ctl = 20L,
    verbose     = TRUE
) {
  set.seed(seed)

  N           <- n_phase1 + n_phase2_trt + n_phase2_ctl
  CUTOFF      <- as.Date("2026-03-05")
  STUDY_START <- as.Date("2022-06-01")
  STUDY_ID    <- "ONCVIZ-001"
  DRUG_NAME   <- "Vizatinib"
  RP2D        <- 300L
  DOSE_LEVELS <- c(100L, 200L, 300L, 400L)

  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

  log  <- function(...) if (verbose) message(sprintf(...))
  clamp     <- function(v, lo, hi) pmax(lo, pmin(hi, v))
  fmt_date  <- function(d) format(as.Date(d), "%Y-%m-%d")
  wc        <- function(opts, wts) sample(opts, 1L, prob = wts / sum(wts))
  rwei      <- function(sh, sc) sc * (-log(max(runif(1), 1e-9)))^(1 / sh)
  mk        <- function(p) ifelse(runif(1) < p, "Y", "N")

  TUMOR_TYPES <- c("NSCLC", "CRC", "HCC", "PDAC", "BRCA")
  SITES       <- sprintf("SITE-%02d", 1:15)
  RESP_NAMES  <- c("CR", "PR", "SD", "PD")

  RR <- list(
    NSCLC = list(TRT = c(0.07,0.38,0.32,0.23), CTL = c(0.01,0.18,0.38,0.43)),
    CRC   = list(TRT = c(0.02,0.10,0.33,0.55), CTL = c(0.00,0.03,0.27,0.70)),
    HCC   = list(TRT = c(0.04,0.22,0.38,0.36), CTL = c(0.01,0.08,0.35,0.56)),
    PDAC  = list(TRT = c(0.01,0.09,0.32,0.58), CTL = c(0.00,0.02,0.25,0.73)),
    BRCA  = list(TRT = c(0.06,0.25,0.37,0.32), CTL = c(0.01,0.12,0.40,0.47))
  )

  SP <- list(
    NSCLC = list(os_sh=0.92,os_sc=24.8,pfs_sh=0.78,pfs_sc=15.2,hr=0.49),
    CRC   = list(os_sh=0.88,os_sc=14.0,pfs_sh=0.82,pfs_sc= 7.8,hr=0.72),
    HCC   = list(os_sh=0.90,os_sc=19.2,pfs_sh=0.80,pfs_sc=10.6,hr=0.58),
    PDAC  = list(os_sh=0.95,os_sc= 8.5,pfs_sh=0.90,pfs_sc= 4.2,hr=0.76),
    BRCA  = list(os_sh=0.85,os_sc=40.0,pfs_sh=0.80,pfs_sc=18.0,hr=0.61)
  )

  TRANS <- list(
    CR = c(CR=0.82,PR=0.10,SD=0.05,PD=0.03),
    PR = c(CR=0.10,PR=0.58,SD=0.20,PD=0.12),
    SD = c(CR=0.02,PR=0.10,SD=0.52,PD=0.36),
    PD = c(CR=0.00,PR=0.00,SD=0.00,PD=1.00)
  )

  MP <- list(
    NSCLC = c(0.46,0.30,0.15,0.17,0.12,0.03,0.08,0.07,0.06,0.12,0.08,0.02,0.07,0.05,0.04),
    CRC   = c(0.60,0.45,0.02,0.03,0.05,0.10,0.03,0.20,0.08,0.05,0.05,0.32,0.12,0.06,0.02),
    HCC   = c(0.28,0.03,0.01,0.02,0.04,0.02,0.05,0.05,0.04,0.08,0.04,0.03,0.18,0.05,0.02),
    PDAC  = c(0.72,0.92,0.02,0.04,0.06,0.03,0.04,0.03,0.03,0.30,0.06,0.22,0.08,0.04,0.02),
    BRCA  = c(0.38,0.03,0.02,0.05,0.03,0.02,0.02,0.35,0.22,0.15,0.15,0.02,0.14,0.04,0.02)
  )

  GENES <- c("TP53","KRAS","EGFR","STK11","KEAP1","BRAF",
             "MET","PIK3CA","PTEN","CDKN2A","RB1","SMAD4",
             "ARID1A","NF1","RET")

  GENE_CHR <- c(TP53="17p13.1",KRAS="12p12.1",EGFR="7p11.2",STK11="19p13.3",
                KEAP1="19p13.2",BRAF="7q34",MET="7q31.2",PIK3CA="3q26.32",
                PTEN="10q23.31",CDKN2A="9p21.3",RB1="13q14.2",SMAD4="18q21.2",
                ARID1A="1p36.11",NF1="17q11.2",RET="10q11.21")

  GENE_DOM <- c(TP53="DNA-binding domain",KRAS="GTPase domain",EGFR="Kinase domain",
                STK11="Kinase domain",KEAP1="Kelch domain",BRAF="Kinase domain",
                MET="Kinase domain",PIK3CA="PI3K domain",PTEN="Phosphatase domain",
                CDKN2A="Ankyrin repeat",RB1="Pocket domain",SMAD4="MH2 domain",
                ARID1A="ARID domain",NF1="RasGAP domain",RET="Kinase domain")

  VC_OPTIONS <- c("Missense_Mutation","Nonsense_Mutation","Frame_Shift_Del",
                  "Frame_Shift_Ins","Splice_Site","Amplification",
                  "Deletion","Fusion","In_Frame_Del")
  VC_WEIGHTS <- c(0.55,0.17,0.12,0.05,0.05,0.03,0.01,0.01,0.01)

  SIG_PROFS <- list(
    NSCLC = list(list(sig="SBS4",desc="Tobacco smoking",w=0.35),
                 list(sig="SBS2",desc="APOBEC activity",w=0.20),
                 list(sig="SBS13",desc="APOBEC activity",w=0.15),
                 list(sig="SBS1",desc="Spontaneous deamination",w=0.18),
                 list(sig="SBS5",desc="Unknown clock-like",w=0.12)),
    CRC   = list(list(sig="SBS1",desc="Spontaneous deamination",w=0.30),
                 list(sig="SBS5",desc="Unknown clock-like",w=0.25),
                 list(sig="SBS15",desc="Defective MMR",w=0.20),
                 list(sig="SBS6",desc="Defective MMR",w=0.15),
                 list(sig="SBS44",desc="Defective MMR",w=0.10)),
    HCC   = list(list(sig="SBS4",desc="Tobacco smoking",w=0.25),
                 list(sig="SBS22",desc="Aristolochic acid",w=0.20),
                 list(sig="SBS24",desc="Aflatoxin exposure",w=0.20),
                 list(sig="SBS1",desc="Spontaneous deamination",w=0.20),
                 list(sig="SBS5",desc="Unknown clock-like",w=0.15)),
    PDAC  = list(list(sig="SBS1",desc="Spontaneous deamination",w=0.35),
                 list(sig="SBS5",desc="Unknown clock-like",w=0.30),
                 list(sig="SBS2",desc="APOBEC activity",w=0.15),
                 list(sig="SBS3",desc="HR deficiency",w=0.12),
                 list(sig="SBS18",desc="ROS damage",w=0.08)),
    BRCA  = list(list(sig="SBS3",desc="HR deficiency / BRCA1/2",w=0.30),
                 list(sig="SBS2",desc="APOBEC activity",w=0.25),
                 list(sig="SBS13",desc="APOBEC activity",w=0.20),
                 list(sig="SBS1",desc="Spontaneous deamination",w=0.15),
                 list(sig="SBS5",desc="Unknown clock-like",w=0.10))
  )

  AE_PROF <- tribble(
    ~SOC, ~PT, ~INC_TRT, ~INC_CTL, ~GW1, ~GW2, ~GW3, ~GW4, ~GW5,
    "Gastrointestinal disorders","Nausea",              0.58,0.32, 0.38,0.40,0.18,0.03,0.01,
    "Gastrointestinal disorders","Diarrhea",            0.52,0.25, 0.32,0.38,0.24,0.05,0.01,
    "Gastrointestinal disorders","Vomiting",            0.34,0.18, 0.38,0.40,0.18,0.03,0.01,
    "Gastrointestinal disorders","Constipation",        0.26,0.20, 0.52,0.36,0.10,0.02,0.00,
    "Gastrointestinal disorders","Abdominal pain",      0.20,0.16, 0.42,0.38,0.16,0.04,0.00,
    "Gastrointestinal disorders","Mucositis",           0.16,0.04, 0.38,0.40,0.18,0.04,0.00,
    "Skin and subcutaneous tissue disorders","Rash",    0.45,0.12, 0.42,0.38,0.16,0.04,0.00,
    "Skin and subcutaneous tissue disorders","Pruritus",0.22,0.08, 0.55,0.35,0.08,0.02,0.00,
    "Skin and subcutaneous tissue disorders","Palmar-plantar erythrodysaesthesia",0.22,0.06,0.40,0.38,0.18,0.04,0.00,
    "General disorders","Fatigue",                     0.65,0.48, 0.30,0.40,0.22,0.07,0.01,
    "General disorders","Decreased appetite",          0.42,0.30, 0.40,0.38,0.18,0.04,0.00,
    "General disorders","Pyrexia",                     0.22,0.16, 0.48,0.36,0.13,0.03,0.00,
    "General disorders","Peripheral edema",            0.18,0.10, 0.52,0.36,0.10,0.02,0.00,
    "Investigations","ALT increased",                  0.32,0.10, 0.40,0.35,0.16,0.07,0.02,
    "Investigations","AST increased",                  0.30,0.08, 0.42,0.35,0.15,0.06,0.02,
    "Investigations","Platelet count decreased",       0.26,0.07, 0.40,0.34,0.18,0.07,0.01,
    "Investigations","Neutrophil count decreased",     0.32,0.08, 0.32,0.30,0.24,0.10,0.04,
    "Investigations","Hemoglobin decreased",           0.38,0.22, 0.38,0.36,0.18,0.07,0.01,
    "Nervous system disorders","Peripheral neuropathy",0.26,0.08, 0.42,0.38,0.16,0.04,0.00,
    "Nervous system disorders","Headache",             0.20,0.16, 0.52,0.36,0.10,0.02,0.00,
    "Musculoskeletal disorders","Arthralgia",           0.24,0.16, 0.44,0.38,0.14,0.04,0.00,
    "Musculoskeletal disorders","Myalgia",              0.18,0.14, 0.48,0.38,0.12,0.02,0.00,
    "Respiratory disorders","Dyspnea",                 0.20,0.16, 0.42,0.38,0.14,0.05,0.01,
    "Respiratory disorders","Pneumonitis",             0.07,0.02, 0.38,0.32,0.18,0.10,0.02,
    "Cardiac disorders","QT prolongation",             0.10,0.04, 0.52,0.34,0.10,0.04,0.00,
    "Immune system disorders","Hypothyroidism",        0.12,0.03, 0.48,0.38,0.12,0.02,0.00
  )

  LAB_TESTS <- tribble(
    ~PARAMCD, ~PARAM,                             ~AVALU,    ~ULN, ~LLN, ~MU,   ~SD,
    "ALT",   "Alanine Aminotransferase",          "U/L",      45,   7,   25,   12,
    "AST",   "Aspartate Aminotransferase",        "U/L",      40,  10,   22,   10,
    "BILI",  "Total Bilirubin",                   "mg/dL",   1.2, 0.2,  0.6,  0.3,
    "ALP",   "Alkaline Phosphatase",              "U/L",     147,  44,   80,   30,
    "GGT",   "Gamma-Glutamyltransferase",         "U/L",      61,   8,   25,   15,
    "CREAT", "Creatinine",                        "mg/dL",   1.2, 0.6,  0.85, 0.2,
    "BUN",   "Blood Urea Nitrogen",               "mg/dL",    20,   7,   13,    4,
    "HGB",   "Haemoglobin",                       "g/dL",     17,  12,   14,  1.5,
    "PLT",   "Platelets",                         "10^9/L",  400, 150,  220,   60,
    "ANC",   "Absolute Neutrophil Count",         "10^9/L",  7.7, 1.8,  3.5,  1.5,
    "WBC",   "White Blood Cells",                 "10^9/L",   11,   4,   6.5,   2,
    "LYMPH", "Lymphocytes",                       "10^9/L",  4.8,   1,   2.0,  0.8,
    "GLUC",  "Glucose",                           "mg/dL",   100,  70,   85,   12,
    "ALB",   "Albumin",                           "g/dL",    5.0, 3.5,  4.2,  0.4,
    "QTCF",  "QTcF",                              "msec",    450, 350,  400,   25,
    "TSH",   "Thyroid Stimulating Hormone",       "mIU/L",   4.0, 0.4,  1.8,  0.8,
    "CD8",   "CD8+ T Cells",                      "%",        40,  10,   22,    8,
    "CD4",   "CD4+ T Cells",                      "%",        60,  20,   38,   10,
    "NK",    "NK Cells",                          "%",        25,   5,   12,    5,
    "TREG",  "Regulatory T Cells",               "%",         8,   1,    3,  1.5
  )

  LB_VISITS  <- c(0L, 42L, 126L, 210L, 336L, 462L, 714L)
  LB_VLABELS <- c("BASELINE","CYCLE 2 DAY 1","CYCLE 4 DAY 1","CYCLE 6 DAY 1",
                  "CYCLE 9 DAY 1","CYCLE 12 DAY 1","CYCLE 18 DAY 1")

  BM_PARAMS <- tribble(
    ~PARAMCD, ~PARAM,                          ~AVALU,     ~MU,  ~SD,
    "TMB",   "Tumor Mutational Burden",        "mut/Mb",    12,    8,
    "PDL1",  "PD-L1 TPS Score",               "%",         30,   25,
    "MSI",   "MSI Sensor Score",              "score",      1.2,   2,
    "CTDNA", "ctDNA Variant Allele Frequency", "%",         30,   20,
    "CEA",   "Carcinoembryonic Antigen",       "ng/mL",      5,   10,
    "CA125", "CA-125",                         "U/mL",      20,   40,
    "IFNg",  "IFN-gamma",                      "pg/mL",     15,   10
  )
  BM_VISITS  <- c(0L, 21L, 42L, 84L, 126L, 168L, 210L)
  BM_VLABELS <- sprintf("CYCLE %d DAY 1", c(1,2,3,5,7,9,11))
  BM_VLABELS[1] <- "BASELINE"

  PRO_PARAMS <- tribble(
    ~PARAMCD, ~PARAM,                   ~SCALETYP,
    "GLQOL",   "Global Health Status",   "functional",
    "PHFUNC",  "Physical Functioning",   "functional",
    "EMOFUNC", "Emotional Functioning",  "functional",
    "FATIGUE", "Fatigue",                "symptom",
    "NAUSEA",  "Nausea/Vomiting",        "symptom",
    "PAIN",    "Pain",                   "symptom",
    "DYSPNEA", "Dyspnoea",               "symptom",
    "APPETIT", "Appetite Loss",          "symptom"
  )
  PRO_MU <- c(65,70,68,35,20,28,25,30)
  PRO_SD <- c(15,18,16,20,18,22,20,22)
  AA     <- strsplit("ACDEFGHIKLMNPQRSTVWY","")[[1]]

  log("Generating ADRAND...")

  n_screened    <- as.integer(N * 1.16)
  n_screen_fail <- n_screened - N
  sf_reasons <- c("Does not meet inclusion criteria",
                  "Laboratory values outside protocol limits",
                  "Prior treatment exclusion","Patient declined",
                  "ECOG performance status >2","Active CNS metastases")

  adrand <- map_dfr(seq_len(n_screened), function(i) {
    site <- sample(SITES, 1L)
    scr_date <- STUDY_START + sample(0:540, 1L)
    if (i <= n_screen_fail) {
      tibble(STUDYID=STUDY_ID, SCRSEQ=i,
             SCRID=sprintf("%s-SCR-%04d",STUDY_ID,i), USUBJID="",
             SITEID=site, SCRDTC=fmt_date(scr_date),
             RANDFL="N", SFFL="Y", SFREASN=sample(sf_reasons,1L),
             ARM="", SAFFL="N", ITTFL="N")
    } else {
      idx <- i - n_screen_fail
      uid <- sprintf("%s-%s-%04d", STUDY_ID, substr(site,6,7), idx)
      arm <- if (idx <= n_phase1) "TREATMENT" else wc(c("TREATMENT","CONTROL"), c(2,1))
      tibble(STUDYID=STUDY_ID, SCRSEQ=i,
             SCRID=sprintf("%s-SCR-%04d",STUDY_ID,i), USUBJID=uid,
             SITEID=site, SCRDTC=fmt_date(scr_date),
             RANDFL="Y", SFFL="N", SFREASN="",
             ARM=arm, SAFFL="Y", ITTFL="Y")
    }
  })

  log("  ADRAND: %d rows | enrolled=%d sf=%d", nrow(adrand), N, n_screen_fail)

  log("Generating ADSL...")

  enrolled <- adrand %>% filter(RANDFL=="Y") %>% arrange(SCRSEQ)

  adsl_raw <- map_dfr(seq_len(N), function(i) {
    rr    <- enrolled[i,]
    uid   <- rr$USUBJID
    arm   <- rr$ARM
    site  <- rr$SITEID
    armcd <- ifelse(arm=="TREATMENT","TRT","CTL")

    dose_mg <- if (arm == "TREATMENT") {
      if (i <= n_phase1) DOSE_LEVELS[min((i-1L) %/% 5L + 1L, 4L)]
      else RP2D
    } else 0L

    tt      <- sample(TUMOR_TYPES, 1L)
    age     <- clamp(round(abs(rnorm(1,63,11))), 28L, 85L)
    sex     <- wc(c("M","F"), c(0.51,0.49))
    race    <- wc(c("WHITE","BLACK","ASIAN","OTHER"), c(0.88,0.098,0.021,0.001))
    smoking <- wc(c("Former heavy","Never","Former light","Current heavy"), c(0.488,0.322,0.178,0.012))
    ecog    <- wc(0:2, c(0.38,0.50,0.12))
    bmi     <- round(clamp(rnorm(1,25.8,4.5), 16, 42), 1)
    t_st    <- wc(c("T1","T2","T3","T4"), c(0.22,0.38,0.28,0.12))
    n_st    <- wc(c("N0","N1","N2","N3"), c(0.28,0.30,0.32,0.10))
    m_st    <- wc(c("M0","M1"), c(0.65,0.35))
    stage   <- if (m_st=="M1") "STAGE IV" else
               if (n_st %in% c("N2","N3")) wc(c("STAGE IIIA","STAGE IIIB"), c(0.65,0.35)) else
               if (n_st=="N1") wc(c("STAGE IIA","STAGE IIB"), c(0.55,0.45)) else
               wc(c("STAGE IA","STAGE IB"), c(0.55,0.45))

    pdl1  <- round(rbeta(1,1.2,2.5)*100, 1)
    pdl1g <- ifelse(pdl1>=50,"HIGH",ifelse(pdl1>=1,"MED","NEG"))
    tmb   <- clamp(round(rlnorm(1,2.1,0.85),1), 0.5, 80)
    msi   <- clamp(round(rexp(1,1/1.2),3), 0, 15)

    sp      <- SP[[tt]]
    os_sc   <- if (arm=="TREATMENT") sp$os_sc  else sp$os_sc  * sp$hr
    pfs_sc  <- if (arm=="TREATMENT") sp$pfs_sc else sp$pfs_sc * sp$hr
    trtsdt  <- STUDY_START + sample(0:540, 1L)
    if (trtsdt > CUTOFF - 60) trtsdt <- CUTOFF - 60
    max_dur <- as.integer(CUTOFF - trtsdt)

    still_on <- runif(1) < ifelse(arm=="TREATMENT", 0.38, 0.18)
    bor  <- sample(RESP_NAMES, 1L, prob=RR[[tt]][[ifelse(arm=="TREATMENT","TRT","CTL")]])
    os_m <- clamp(rwei(sp$os_sh, os_sc), 1, 120)
    pfs_m<- clamp(rwei(sp$pfs_sh, pfs_sc), 0.5, os_m)
    os_d <- as.integer(os_m * 30.44)
    pfs_d<- as.integer(pfs_m * 30.44)
    dur  <- if (still_on) max_dur else clamp(os_d, 42L, max_dur)
    trtedt <- trtsdt + dur

    os_ev  <- if (still_on) 0L else ifelse(os_d <= max_dur, 1L, 0L)
    ncd    <- (!still_on) && os_ev==1 && runif(1) < 0.05
    pfs_ev <- if (ncd || still_on) 0L else ifelse(pfs_d <= dur, 1L, 0L)
    comp_event <- if (ncd) 1L else if (pfs_ev==1L) 1L else if (os_ev==1L) 1L else 0L
    comp_type  <- if (ncd) "DEATH_WITHOUT_PROGRESSION" else
                  if (pfs_ev==1L) "PROGRESSION" else
                  if (os_ev==1L)  "DEATH_WITHOUT_PROGRESSION" else "CENSORED"

    dcsreas_trt <- c("Adverse Event","Progressive Disease","Withdrawal of Consent","Death","Protocol Deviation")
    dcsreas_ctl <- c("Progressive Disease","Withdrawal of Consent","Death","Crossover to Treatment","Protocol Deviation")
    dcsreas <- if (still_on) "" else sample(if(arm=="TREATMENT") dcsreas_trt else dcsreas_ctl, 1L)

    mp <- MP[[tt]]

    tibble(
      STUDYID=STUDY_ID, USUBJID=uid, SITEID=site,
      ARM=arm, ARMCD=armcd,
      TRT01P=ifelse(arm=="TREATMENT", paste(DRUG_NAME,dose_mg,"mg"), "Placebo"),
      PHASE=ifelse(i<=n_phase1,"I","II"), DOSELEVEL=dose_mg,
      TUMORTYPE=tt, AGE=age, AGEGR1=ifelse(age>=65,">=65","<65"),
      SEX=sex, RACE=race, SMOKING=smoking, ECOG=ecog, BMI=bmi,
      T_STAGE=t_st, N_STAGE=n_st, M_STAGE=m_st, STAGE=stage,
      LIVERMETS=mk(0.25), BRAINMETS=mk(0.08),
      PRIORSURG=mk(0.45), PRIORRAD=mk(0.38),
      PRIORLINES=wc(0:3, c(0.15,0.45,0.30,0.10)),
      BASESZ=round(rlnorm(1,3.5,0.6),1),
      PDL1SCORE=pdl1, PDL1GRP=pdl1g,
      TMB=tmb, TMBHIGH=ifelse(tmb>=10,"Y","N"),
      MSI_SCORE=msi, MSISTS=ifelse(msi>=3.5,"MSI-H","MSS"),
      TP53MUT =mk(mp[1]),  KRASMUT =mk(mp[2]),  EGFRMUT =mk(mp[3]),
      STK11MUT=mk(mp[4]),  KEAP1MUT=mk(mp[5]),  BRAFMUT =mk(mp[6]),
      METAMP  =mk(mp[7]),  PIK3CAMT=mk(mp[8]),  PTENMUT =mk(mp[9]),
      CDKN2AMT=mk(mp[10]), RB1MUT  =mk(mp[11]), SMAD4MUT=mk(mp[12]),
      ARID1AMT=mk(mp[13]), NF1MUT  =mk(mp[14]), RETFUS  =mk(mp[15]),
      BESTRSPC=bor,
      TRTSDT=fmt_date(trtsdt), TRTEDT=fmt_date(trtedt),
      TRTDURD=dur, EOSSTT=ifelse(still_on,"ONGOING","DISCONTINUED"), DCSREAS=dcsreas,
      OSDUR=os_d, OSCR=os_ev, OSDTC=fmt_date(trtsdt+os_d),
      PFSDUR=pfs_d, PFSCR=pfs_ev, PFSDTC=fmt_date(trtsdt+pfs_d),
      COMPEVENT=comp_event, COMPTYPE=comp_type,
      CUTDTC=fmt_date(CUTOFF), SAFFL="Y", ITTFL="Y",
      PPROTFL=ifelse(runif(1)>0.12,"Y","N"),
      .trtsdt=trtsdt, .dur=dur, .os_d=os_d, .pfs_d=pfs_d,
      .bor=bor, .still_on=still_on, .arm_key=ifelse(arm=="TREATMENT","TRT","CTL")
    )
  })

  log("  ADSL: %d | TRT=%d CTL=%d", nrow(adsl_raw),
      sum(adsl_raw$ARM=="TREATMENT"), sum(adsl_raw$ARM=="CONTROL"))

  log("Generating ADTR (SLD)...")

  adtr <- map_dfr(seq_len(nrow(adsl_raw)), function(i) {
    s      <- adsl_raw[i,]
    bl     <- s$BASESZ
    bor    <- s$.bor
    trtsdt <- s$.trtsdt
    dur    <- s$.dur
    n_vis  <- max(2L, min(dur %/% 42L, 12L))

    final_pchg <- switch(bor,
      CR = -100 + runif(1)*5,
      PR = runif(1, -80, -30),
      SD = runif(1, -29, 19),
      PD = runif(1, 20, 120)
    )
    nadir_vis <- min(3L, n_vis - 1L)

    rows <- list()
    rows[[1]] <- tibble(
      STUDYID=STUDY_ID, USUBJID=s$USUBJID, ARM=s$ARM, TUMORTYPE=s$TUMORTYPE,
      PARAMCD="SUMDIAM", PARAM="Sum of Longest Diameters (mm)",
      AVISIT="BASELINE", AVISITN=0L,
      ADT=fmt_date(trtsdt), ADTN=0L,
      BASE=bl, AVAL=bl, CHG=0.0, PCHG=0.0, ANL01FL="Y"
    )
    for (k in seq_len(n_vis)) {
      aday <- 42L * k
      adt  <- trtsdt + aday
      if (adt > CUTOFF) break
      t    <- k / max(n_vis, 1)
      pchg <- if (bor %in% c("CR","PR","SD")) {
        if (k <= nadir_vis)
          final_pchg * (k / nadir_vis) + rnorm(1,0,4)
        else {
          rec <- (k - nadir_vis) / max(n_vis - nadir_vis, 1)
          final_pchg + (final_pchg * ifelse(bor=="SD",1.05,0.95) - final_pchg)*rec + rnorm(1,0,4)
        }
      } else {
        final_pchg * t + rnorm(1,0,5)
      }
      pchg <- round(clamp(pchg, -100, 300), 1)
      aval <- round(max(0, bl*(1+pchg/100)), 1)
      rows[[k+1L]] <- tibble(
        STUDYID=STUDY_ID, USUBJID=s$USUBJID, ARM=s$ARM, TUMORTYPE=s$TUMORTYPE,
        PARAMCD="SUMDIAM", PARAM="Sum of Longest Diameters (mm)",
        AVISIT=sprintf("CYCLE %d DAY 1",k), AVISITN=k,
        ADT=fmt_date(adt), ADTN=aday,
        BASE=bl, AVAL=aval, CHG=round(aval-bl,1), PCHG=pchg, ANL01FL="Y"
      )
    }
    bind_rows(rows)
  })

  log("  ADTR: %d rows | median visits/pt=%.0f",
      nrow(adtr), median(table(adtr$USUBJID)))

  log("Generating ADRS (longitudinal, aligned with ADTR)...")

  resp_rank <- c(CR=1,PR=2,SD=3,PD=4)

  adrs <- map_dfr(seq_len(nrow(adsl_raw)), function(i) {
    s      <- adsl_raw[i,]
    uid    <- s$USUBJID
    trtsdt <- s$.trtsdt
    visits <- adtr %>% filter(USUBJID==uid) %>% pull(AVISITN)
    prev   <- s$.bor

    map_dfr(visits, function(v) {
      adt    <- trtsdt + v*42L
      avisit <- if (v==0L) "BASELINE" else sprintf("CYCLE %d DAY 1",v)
      resp   <- if (v==0L) "NE" else {
        r <- sample(RESP_NAMES, 1L, prob=TRANS[[prev]])
        prev <<- r
        r
      }
      tibble(
        STUDYID=STUDY_ID, USUBJID=uid, ARM=s$ARM, TUMORTYPE=s$TUMORTYPE,
        PARAMCD="OVRLRESP", PARAM="Overall Response per RECIST 1.1",
        AVISIT=avisit, AVISITN=v,
        ADT=fmt_date(adt), ADTN=v*42L,
        AVALC=resp, AVAL=ifelse(resp=="NE",5L,resp_rank[resp]),
        BICR_CONF=ifelse(runif(1)>0.15,"Y","N"), ANL01FL="Y"
      )
    })
  })

  bor_from_adrs <- adrs %>%
    filter(AVALC %in% RESP_NAMES) %>%
    mutate(rk=resp_rank[AVALC]) %>%
    group_by(USUBJID) %>%
    slice_min(rk, n=1L, with_ties=FALSE) %>%
    ungroup() %>%
    select(USUBJID, BESTRSPC_ADJ=AVALC)

  adsl_raw <- adsl_raw %>%
    left_join(bor_from_adrs, by="USUBJID") %>%
    mutate(BESTRSPC=coalesce(BESTRSPC_ADJ, BESTRSPC),
           .bor=coalesce(BESTRSPC_ADJ, .bor)) %>%
    select(-BESTRSPC_ADJ)

  adsl <- adsl_raw %>% select(-starts_with("."))

  log("  ADRS: %d rows | median visits/pt=%.0f | BOR reconciled",
      nrow(adrs), median(table(adrs$USUBJID)))

  responders <- adsl_raw %>% filter(BESTRSPC %in% c("CR","PR")) %>% pull(USUBJID)

  log("Generating ADTTE...")

  SGVARS <- c("TUMORTYPE","AGEGR1","SEX","ECOG","PDL1GRP","MSISTS",
              "TMBHIGH","LIVERMETS","PRIORLINES","EGFRMUT","KRASMUT",
              "TP53MUT","STK11MUT","STAGE")

  adtte <- map_dfr(seq_len(nrow(adsl_raw)), function(i) {
    s      <- adsl_raw[i,]
    uid    <- s$USUBJID
    arm    <- s$ARM
    armcd  <- s$ARMCD
    trtsdt <- s$.trtsdt
    os_d   <- s$.os_d
    pfs_d  <- s$.pfs_d
    bor    <- s$BESTRSPC
    dur    <- s$.dur
    max_dur<- as.integer(CUTOFF - trtsdt)
    os_ev  <- s$OSCR
    pfs_ev <- s$PFSCR
    ncd    <- s$COMPTYPE == "DEATH_WITHOUT_PROGRESSION" && !s$.still_on
    sg     <- s[SGVARS]

    mrow <- function(paramcd, param, aval, cnsr, evnt, adt) {
      bind_cols(
        tibble(
          STUDYID=STUDY_ID, USUBJID=uid, ARM=arm, ARMCD=armcd,
          PARAMCD=paramcd, PARAM=param,
          AVAL=as.integer(aval), AVALU="DAYS",
          AVALM=round(aval/30.44,2),
          CNSR=cnsr, EVNTDESC=evnt,
          ADT=fmt_date(adt), STARTDT=s$TRTSDT,
          LM6MFL =ifelse(aval>=182,"Y","N"),
          LM12MFL=ifelse(aval>=365,"Y","N"),
          LM24MFL=ifelse(aval>=730,"Y","N"),
          COMPEVENT=s$COMPEVENT, COMPTYPE=s$COMPTYPE,
          ANL01FL="Y"
        ), sg
      )
    }

    rows <- list(
      mrow("OS","Overall Survival",
           min(os_d,max_dur), 1L-os_ev,
           ifelse(os_ev==1L,"DEATH","CENSORED"),
           trtsdt+min(os_d,max_dur)),
      mrow("PFS","Progression-Free Survival",
           min(pfs_d,dur), 1L-pfs_ev,
           ifelse(pfs_ev==1L,"PROGRESSION","CENSORED"),
           trtsdt+min(pfs_d,dur)),
      mrow("TTP","Time to Progression",
           min(pfs_d,dur), ifelse(ncd||!pfs_ev,1L,0L),
           ifelse(pfs_ev&&!ncd,"PROGRESSION","CENSORED"),
           trtsdt+min(pfs_d,dur)),
      mrow("EFS","Event-Free Survival",
           min(pfs_d,dur), 1L-pfs_ev,
           ifelse(pfs_ev==1L,"EVENT","CENSORED"),
           trtsdt+min(pfs_d,dur))
    )

    if (uid %in% responders) {
      ttr_d <- clamp(as.integer(rexp(1,1/45)), 14L, pfs_d)
      dor_d <- clamp(pfs_d - ttr_d, 21L, dur)
      dor_ev<- ifelse(pfs_ev==1L,1L,0L)
      resp_date <- trtsdt + ttr_d
      rows[[5]] <- bind_cols(
        tibble(
          STUDYID=STUDY_ID, USUBJID=uid, ARM=arm, ARMCD=armcd,
          PARAMCD="DOR", PARAM="Duration of Response",
          AVAL=dor_d, AVALU="DAYS", AVALM=round(dor_d/30.44,2),
          CNSR=1L-dor_ev, EVNTDESC=ifelse(dor_ev==1L,"PROGRESSION","CENSORED"),
          ADT=fmt_date(resp_date+dor_d), STARTDT=fmt_date(resp_date),
          LM6MFL =ifelse(dor_d>=182,"Y","N"),
          LM12MFL=ifelse(dor_d>=365,"Y","N"),
          LM24MFL=ifelse(dor_d>=730,"Y","N"),
          COMPEVENT=s$COMPEVENT, COMPTYPE=s$COMPTYPE, ANL01FL="Y"
        ), sg
      )
      rows[[6]] <- bind_cols(
        tibble(
          STUDYID=STUDY_ID, USUBJID=uid, ARM=arm, ARMCD=armcd,
          PARAMCD="TTR", PARAM="Time to Response",
          AVAL=ttr_d, AVALU="DAYS", AVALM=round(ttr_d/30.44,2),
          CNSR=0L, EVNTDESC="RESPONSE",
          ADT=fmt_date(trtsdt+ttr_d), STARTDT=s$TRTSDT,
          LM6MFL=ifelse(ttr_d>=182,"Y","N"),
          LM12MFL=ifelse(ttr_d>=365,"Y","N"),
          LM24MFL=ifelse(ttr_d>=730,"Y","N"),
          COMPEVENT=s$COMPEVENT, COMPTYPE=s$COMPTYPE, ANL01FL="Y"
        ), sg
      )
      if (bor == "CR") {
        dfs_start <- trtsdt + ttr_d
        dfs_d     <- clamp(os_d - ttr_d, 30L, max_dur)
        rows[[7]] <- bind_cols(
          tibble(
            STUDYID=STUDY_ID, USUBJID=uid, ARM=arm, ARMCD=armcd,
            PARAMCD="DFS", PARAM="Disease-Free Survival",
            AVAL=dfs_d, AVALU="DAYS", AVALM=round(dfs_d/30.44,2),
            CNSR=1L-os_ev, EVNTDESC=ifelse(os_ev==1L,"RELAPSE/DEATH","CENSORED"),
            ADT=fmt_date(dfs_start+dfs_d), STARTDT=fmt_date(dfs_start),
            LM6MFL =ifelse(dfs_d>=182,"Y","N"),
            LM12MFL=ifelse(dfs_d>=365,"Y","N"),
            LM24MFL=ifelse(dfs_d>=730,"Y","N"),
            COMPEVENT=s$COMPEVENT, COMPTYPE=s$COMPTYPE, ANL01FL="Y"
          ), sg
        )
      }
    }

    bind_rows(rows)
  })

  log("  ADTTE: %d rows | params: %s",
      nrow(adtte), paste(sort(unique(adtte$PARAMCD)), collapse=","))

  log("Generating ADAE...")

  adae <- map_dfr(seq_len(nrow(adsl_raw)), function(i) {
    s      <- adsl_raw[i,]
    uid    <- s$USUBJID
    arm    <- s$ARM
    trtsdt <- s$.trtsdt
    dur    <- max(1L, s$.dur)
    seq_n  <- 0L
    rows   <- list()

    for (j in seq_len(nrow(AE_PROF))) {
      ae  <- AE_PROF[j,]
      inc <- clamp(ifelse(arm=="TREATMENT", ae$INC_TRT, ae$INC_CTL), 0, 0.95)
      if (runif(1) > inc) next
      seq_n <- seq_n + 1L
      onset <- sample(seq_len(dur), 1L)
      grade <- sample(1:5, 1L, prob=c(ae$GW1,ae$GW2,ae$GW3,ae$GW4,ae$GW5))
      dur_ae <- max(1L, as.integer(rexp(1, 1/(10+grade*5))))
      dur_ae <- min(dur_ae, dur - onset + 30L)
      action <- if (grade>=3) sample(c("DOSE REDUCED","DRUG INTERRUPTED","DRUG WITHDRAWN","NONE"),
                                     1L, prob=c(0.22,0.25,0.08,0.45)) else "NONE"
      outcome <- sample(c("RECOVERED","RECOVERING","NOT RECOVERED","FATAL"),
                        1L, prob=if(grade>=4) c(0.60,0.18,0.18,0.04) else c(0.80,0.12,0.07,0.01))
      rows[[seq_n]] <- tibble(
        STUDYID=STUDY_ID, USUBJID=uid, ARM=arm, TUMORTYPE=s$TUMORTYPE,
        AESEQ=seq_n, AESOC=ae$SOC, AEPT=ae$PT, AETOXGR=grade,
        AESER=ifelse(grade>=3,"Y","N"),
        AEREL=ifelse(runif(1)<ifelse(arm=="TREATMENT",0.60,0.25),"RELATED","NOT RELATED"),
        AEACN=action, AEOUT=outcome,
        AESTDTC=fmt_date(trtsdt+onset),
        AEENDTC=fmt_date(min(trtsdt+onset+dur_ae, CUTOFF)),
        AESTDY=onset, AEENDY=onset+dur_ae, AEDUR=dur_ae,
        DOSELEVEL=s$DOSELEVEL, CTCAE_V="5.0"
      )
    }
    if (length(rows)==0L) return(NULL)
    bind_rows(rows)
  })

  log("  ADAE: %d rows | Grade>=3: %d", nrow(adae), sum(adae$AETOXGR>=3))

  log("Generating ADEX...")

  adex <- map_dfr(seq_len(nrow(adsl_raw)), function(i) {
    s      <- adsl_raw[i,]
    uid    <- s$USUBJID
    arm    <- s$ARM
    trtsdt <- s$.trtsdt
    dur    <- max(1L, s$.dur)
    planned<- s$DOSELEVEL
    level  <- 1.0
    cum    <- 0.0
    rows   <- list()

    for (cyc in seq_len(max(1L, dur %/% 21L))) {
      cs <- trtsdt + (cyc-1L)*21L
      ce <- cs + 20L
      if (ce > CUTOFF) break
      mod_reason <- ""
      if (cyc > 2L && arm=="TREATMENT" && planned > 0L) {
        if (runif(1) < 0.08) {
          level <- max(0.33, level - 0.33)
          mod_reason <- sample(c("AE - Hepatotoxicity","AE - Diarrhea","AE - Fatigue"), 1L)
        } else if (runif(1) < 0.03) {
          level <- 0.0
          mod_reason <- "AE - Dose Interruption"
        }
      }
      actual   <- round(planned * level)
      days_on  <- if (level > 0) 21L else sample(0:10, 1L)
      doseint  <- if (planned > 0L) round(actual * days_on / (planned * 21), 4) else 0.0
      cum      <- cum + actual * days_on
      rows[[cyc]] <- tibble(
        STUDYID=STUDY_ID, USUBJID=uid, ARM=arm, TUMORTYPE=s$TUMORTYPE,
        EXSEQ=cyc, EXTRT=ifelse(arm=="TREATMENT",DRUG_NAME,"Placebo"),
        EXDOSE=actual, EXDOSU="mg", EXDOSFRQ="QD", EXROUTE="ORAL",
        EXSTDTC=fmt_date(cs), EXENDTC=fmt_date(ce),
        EXCYCLE=cyc, PLANDOSE=planned,
        DOSELEVEL=s$DOSELEVEL, DAYSONDRUG=days_on,
        DOSEINT=doseint, CUMDOSE=round(cum,0),
        MODFL=ifelse(nchar(mod_reason)>0,"Y","N"), MODREASN=mod_reason
      )
    }
    bind_rows(rows)
  })

  log("  ADEX: %d rows | DOSEINT range=%.2f–%.2f",
      nrow(adex), min(adex$DOSEINT), max(adex$DOSEINT))

  log("Generating ADLB...")

  adlb <- map_dfr(seq_len(nrow(adsl_raw)), function(i) {
    s      <- adsl_raw[i,]
    uid    <- s$USUBJID
    trtsdt <- s$.trtsdt
    dur    <- s$.dur

    map_dfr(seq_len(nrow(LAB_TESTS)), function(j) {
      lb   <- LAB_TESTS[j,]
      base <- max(lb$LLN*0.5, abs(rnorm(1, lb$MU, lb$SD)))

      map_dfr(seq_along(LB_VISITS), function(vi) {
        vday <- LB_VISITS[vi]
        if (vday > dur + 14L) return(NULL)
        adt <- trtsdt + vday
        if (adt > CUTOFF) return(NULL)
        drift <- rnorm(1, 0, lb$SD*0.12) * vday / max(dur,1) * 100
        aval  <- max(0.01, base + drift + rnorm(1,0,lb$SD*0.10))
        chg   <- round(aval - base, 3)
        pchg  <- round(chg / base * 100, 1)
        xuln  <- round(aval / lb$ULN, 3)
        atxg  <- if (aval > 5*lb$ULN) 3L else if (aval > 3*lb$ULN) 2L else
                 if (aval > lb$ULN)   1L else 0L
        tibble(
          STUDYID=STUDY_ID, USUBJID=uid, ARM=s$ARM, TUMORTYPE=s$TUMORTYPE,
          PARAMCD=lb$PARAMCD, PARAM=lb$PARAM, AVALU=lb$AVALU,
          AVISIT=LB_VLABELS[vi], AVISITN=vi-1L,
          ADT=fmt_date(adt), ADTN=vday,
          BASE=round(base,3), AVAL=round(aval,3), CHG=chg, PCHG=pchg,
          ULN=lb$ULN, LLN=lb$LLN, XULNFL=xuln, ATOXGR=atxg,
          ANRIND=ifelse(aval>lb$ULN,"HIGH",ifelse(aval<lb$LLN,"LOW","NORMAL")),
          ANL01FL="Y"
        )
      })
    })
  })

  log("  ADLB: %d rows | %d params", nrow(adlb), n_distinct(adlb$PARAMCD))

  log("Generating ADBM...")

  adbm <- map_dfr(seq_len(nrow(adsl_raw)), function(i) {
    s      <- adsl_raw[i,]
    uid    <- s$USUBJID
    trtsdt <- s$.trtsdt
    dur    <- s$.dur
    bor    <- s$BESTRSPC

    map_dfr(seq_len(nrow(BM_PARAMS)), function(j) {
      bm   <- BM_PARAMS[j,]
      base <- clamp(abs(rnorm(1, bm$MU, bm$SD)), 0, bm$MU*4)

      map_dfr(seq_along(BM_VISITS), function(vi) {
        vday <- BM_VISITS[vi]
        if (vday > dur + 14L) return(NULL)
        adt <- trtsdt + vday
        if (adt > CUTOFF) return(NULL)
        decay <- if (bm$PARAMCD == "CTDNA") {
          if (bor %in% c("CR","PR")) exp(-0.4*(vday/30.44)) else
          if (bor == "PD") 1 + 0.15*(vday/30.44) else 1
        } else 1
        aval <- clamp(base*decay + rnorm(1,0,bm$SD*0.15), 0, bm$MU*5)
        chg  <- round(aval - base, 3)
        pchg <- round(chg / max(base, 0.001) * 100, 1)
        tibble(
          STUDYID=STUDY_ID, USUBJID=uid, ARM=s$ARM, TUMORTYPE=s$TUMORTYPE,
          PARAMCD=bm$PARAMCD, PARAM=bm$PARAM, AVALU=bm$AVALU,
          AVISIT=BM_VLABELS[vi], AVISITN=vi-1L,
          ADT=fmt_date(adt), ADTN=vday,
          BASE=round(base,3), AVAL=round(aval,3), CHG=chg, PCHG=pchg,
          BESTRSPC=bor, ANL01FL="Y"
        )
      })
    })
  })

  log("  ADBM: %d rows | %d params", nrow(adbm), n_distinct(adbm$PARAMCD))

  log("Generating ADPK...")

  PK_TIMES  <- c(0, 0.5, 1, 2, 3, 4, 6, 8, 12, 24)
  PK_VISITS <- c("CYCLE 1 DAY 1","CYCLE 1 DAY 15","CYCLE 3 DAY 1")
  PK_DAYS   <- c(0L, 14L, 56L)
  KA_POP    <- 1.5; CL_POP <- 15.0; VD_POP <- 100.0

  adpk <- map_dfr(adsl_raw %>% filter(ARM=="TREATMENT") %>% seq_len(nrow(.)), function(i) {
    s      <- (adsl_raw %>% filter(ARM=="TREATMENT"))[i,]
    uid    <- s$USUBJID
    trtsdt <- s$.trtsdt
    cl_i   <- CL_POP * exp(rnorm(1,0,0.35))
    vd_i   <- VD_POP * exp(rnorm(1,0,0.28))
    ka_i   <- KA_POP * exp(rnorm(1,0,0.30))
    ke_i   <- cl_i / vd_i

    conc_rows <- map_dfr(seq_along(PK_VISITS), function(vi) {
      bdt <- trtsdt + PK_DAYS[vi]
      if (bdt > CUTOFF) return(NULL)
      map_dfr(PK_TIMES, function(t) {
        conc <- if (t==0) 0 else {
          raw <- (RP2D * ka_i / (vd_i*(ka_i-ke_i))) * (exp(-ke_i*t) - exp(-ka_i*t))
          max(0, raw) * exp(rnorm(1,0,0.14))
        }
        tibble(
          STUDYID=STUDY_ID, USUBJID=uid, ARM="TREATMENT", TUMORTYPE=s$TUMORTYPE,
          PARAMCD="CONC", PARAM=paste(DRUG_NAME,"Plasma Concentration"),
          AVISIT=PK_VISITS[vi], AVISITN=vi, NOMTPT=t,
          ADT=fmt_date(bdt), AVAL=round(conc,3), AVALU="ng/mL",
          DOSE=RP2D, ROUTE="ORAL", ANL01FL="Y", AGE=s$AGE, SEX=s$SEX
        )
      })
    })

    cmax   <- max(0.1, (RP2D * ka_i/(vd_i*(ka_i-ke_i))) *
                 (exp(-ke_i*3)-exp(-ka_i*3)) * exp(rnorm(1,0,0.10)))
    auc    <- max(0.1, RP2D / cl_i * exp(rnorm(1,0,0.10)))
    aucinf <- auc + (conc_rows %>% filter(NOMTPT==max(PK_TIMES)) %>%
                     pull(AVAL) %>% mean(na.rm=TRUE)) / ke_i
    thalf  <- 0.693 * vd_i / cl_i
    trough <- conc_rows %>% filter(NOMTPT==max(PK_TIMES)) %>%
              pull(AVAL) %>% mean(na.rm=TRUE)

    param_rows <- tibble(
      STUDYID=STUDY_ID, USUBJID=uid, ARM="TREATMENT", TUMORTYPE=s$TUMORTYPE,
      PARAMCD=c("CMAX","AUC","AUCINF","TMAX","THALF","TROUGH"),
      PARAM=paste(DRUG_NAME,c("Cmax","AUClast","AUCinf","Tmax","t1/2","Trough Conc")),
      AVISIT="CYCLE 1 DAY 1", AVISITN=1L, NOMTPT=NA_real_,
      ADT=fmt_date(trtsdt),
      AVAL=round(c(cmax,auc,aucinf,3+rnorm(1,0,0.4),thalf,trough),3),
      AVALU=c("ng/mL","ng*h/mL","ng*h/mL","h","h","ng/mL"),
      DOSE=RP2D, ROUTE="ORAL", ANL01FL="Y", AGE=s$AGE, SEX=s$SEX
    )
    bind_rows(conc_rows, param_rows)
  })

  log("  ADPK: %d rows | params: %s", nrow(adpk),
      paste(sort(unique(adpk$PARAMCD)),collapse=","))

  log("Generating ADMUT...")

  admut <- map_dfr(seq_len(nrow(adsl_raw)), function(i) {
    s   <- adsl_raw[i,]
    uid <- s$USUBJID
    tt  <- s$TUMORTYPE
    mp  <- MP[[tt]]
    if (runif(1) > 0.92) return(NULL)

    map_dfr(seq_along(GENES), function(gi) {
      if (runif(1) > mp[gi]) return(NULL)
      gene <- GENES[gi]
      vc   <- sample(VC_OPTIONS, 1L, prob=VC_WEIGHTS)
      pos  <- sample(1:1200, 1L)
      ref  <- sample(AA, 1L)
      alt  <- sample(AA[AA!=ref], 1L)
      hgvsp <- if (grepl("Missense",vc)) sprintf("p.%s%d%s",ref,pos,alt) else
               if (grepl("Frame",vc))    sprintf("p.%s%dfs",ref,pos) else
               if (grepl("Nonsense",vc)) sprintf("p.%s%d*",ref,pos) else
               sprintf("%s %s",gene,vc)
      depth <- sample(80:400,1L)
      vaf   <- clamp(runif(1,0.05,0.55)+rnorm(1,0,0.03), 0.01, 0.85)
      alt_c <- round(vaf*depth)
      tibble(
        STUDYID=STUDY_ID, USUBJID=uid, ARM=s$ARM, TUMORTYPE=tt,
        HUGO_SYMBOL=gene, CHROMOSOME=GENE_CHR[gene],
        VARIANT_CLASS=vc, VARIANT_TYPE=sub("_Mutation","",vc),
        HGVSP=hgvsp, PROTEIN_POS=pos, PROTEIN_DOMAIN=GENE_DOM[gene],
        REF_ALLELE=ref, ALT_ALLELE=alt,
        VAF=round(vaf,3), T_DEPTH=depth,
        T_ALT_COUNT=alt_c, T_REF_COUNT=depth-alt_c,
        CLONAL=mk(0.65), IMPACT=wc(c("HIGH","MODERATE","LOW"),c(0.35,0.50,0.15)),
        BESTRSPC=s$BESTRSPC, TMB=s$TMB, TMBHIGH=s$TMBHIGH, ANL01FL="Y"
      )
    })
  })

  log("  ADMUT: %d rows | %d genes | %d patients",
      nrow(admut), n_distinct(admut$HUGO_SYMBOL), n_distinct(admut$USUBJID))

  log("Generating ADSIG...")

  adsig <- map_dfr(seq_len(nrow(adsl_raw)), function(i) {
    s    <- adsl_raw[i,]
    uid  <- s$USUBJID
    tt   <- s$TUMORTYPE
    sigs <- SIG_PROFS[[tt]]
    n_muts <- max(10L, as.integer(s$TMB * 50))
    raw_w  <- sapply(sigs, function(sg) max(0.01, sg$w + rnorm(1,0,0.05)))
    norm_w <- raw_w / sum(raw_w)
    map_dfr(seq_along(sigs), function(si) {
      sg <- sigs[[si]]
      tibble(
        STUDYID=STUDY_ID, USUBJID=uid, ARM=s$ARM, TUMORTYPE=tt,
        SIG_NAME=sg$sig, SIG_DESC=sg$desc,
        SIG_WEIGHT=round(norm_w[si],4),
        N_MUTS=round(n_muts*norm_w[si]),
        TMB=s$TMB, TMBHIGH=s$TMBHIGH, MSISTS=s$MSISTS,
        BESTRSPC=s$BESTRSPC, ANL01FL="Y"
      )
    })
  })

  log("  ADSIG: %d rows | %d signatures", nrow(adsig), n_distinct(adsig$SIG_NAME))

  log("Generating ADPR...")

  adpr <- map_dfr(seq_len(nrow(adsl_raw)), function(i) {
    s      <- adsl_raw[i,]
    uid    <- s$USUBJID
    trtsdt <- s$.trtsdt
    dur    <- s$.dur
    bor    <- s$BESTRSPC
    map_dfr(seq_len(nrow(PRO_PARAMS)), function(j) {
      pp   <- PRO_PARAMS[j,]
      base <- clamp(rnorm(1, PRO_MU[j], PRO_SD[j]), 0, 100)
      is_sx <- pp$SCALETYP == "symptom"
      map_dfr(seq(0L, min(168L, dur), by=42L), function(vday) {
        adt <- trtsdt + vday
        if (adt > CUTOFF) return(NULL)
        delta <- if (bor %in% c("CR","PR") && vday>0)
                   ifelse(is_sx,-5,5)*min(vday/120,1) else
                   ifelse(is_sx,3,-3)*min(vday/168,1)
        aval  <- clamp(round(base+delta+rnorm(1,0,PRO_SD[j]*0.18),1), 0, 100)
        chg   <- round(aval-base,1)
        pchg  <- round(chg/max(base,0.001)*100,1)
        vis   <- if (vday==0L) "BASELINE" else sprintf("CYCLE %d DAY 1", vday%/%42L)
        tibble(
          STUDYID=STUDY_ID, USUBJID=uid, ARM=s$ARM, TUMORTYPE=s$TUMORTYPE,
          PARAMCD=pp$PARAMCD, PARAM=pp$PARAM, SCALETYP=pp$SCALETYP,
          AVISIT=vis, AVISITN=vday%/%42L,
          ADT=fmt_date(adt), ADTN=vday,
          BASE=round(base,1), AVAL=aval, CHG=chg, PCHG=pchg,
          MID=10.0,
          MIDRESP=ifelse(abs(chg)>=10,"Y","N"),
          DETRFL=ifelse((is_sx&&chg>=10)||(!is_sx&&chg<=-10),"Y","N"),
          ANL01FL="Y"
        )
      })
    })
  })

  log("  ADPR: %d rows | %d params", nrow(adpr), n_distinct(adpr$PARAMCD))

  log("\nSaving CSVs to: %s", normalizePath(output_dir, mustWork=FALSE))

  datasets <- list(
    ADRAND=adrand, ADSL=adsl,   ADRS=adrs,  ADTR=adtr,
    ADTTE =adtte,  ADAE=adae,   ADEX=adex,  ADLB=adlb,
    ADBM  =adbm,   ADPK=adpk,  ADMUT=admut, ADSIG=adsig,
    ADPR  =adpr
  )

  total <- 0L
  for (nm in names(datasets)) {
    path <- file.path(output_dir, paste0(nm, ".csv"))
    write.csv(datasets[[nm]], path, row.names=FALSE)
    log("  %-8s %6d rows x %d cols", nm, nrow(datasets[[nm]]), ncol(datasets[[nm]]))
    total <- total + nrow(datasets[[nm]])
  }

  zip_path <- paste0(output_dir, ".zip")
  csv_files <- file.path(output_dir, paste0(names(datasets),".csv"))
  zip(zipfile=zip_path, files=csv_files, flags="-j")

  log("\n============================================================")
  log("  DONE")
  log("  Total rows : %s", format(total, big.mark=","))
  log("  Patients   : %d (Phase I=%d TRT, Phase II=%d TRT + %d CTL)",
      N, n_phase1, n_phase2_trt, n_phase2_ctl)
  log("  Responders : %d (CR+PR)", length(responders))
  log("  Study      : %s | %s", STUDY_ID, DRUG_NAME)
  log("  Seed       : %d", seed)
  log("  CSV dir    : %s", normalizePath(output_dir, mustWork=FALSE))
  log("  Zip file   : %s", normalizePath(zip_path,   mustWork=FALSE))
  log("============================================================")

  invisible(datasets)
}

datasets <- generate_adam_oncviz001(
  output_dir   = "./data_v2",
  seed         = 42L,
  n_phase1     = 20L,
  n_phase2_trt = 40L,
  n_phase2_ctl = 20L,
  verbose      = TRUE
)
