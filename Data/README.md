# ONCVIZ-001 · ADaM v1
### A Tumor-Stratified Synthetic Clinical Trial Dataset for Oncology Visualization

```
Study     ONCVIZ-001  ·  Vizatinib 300 mg QD vs Placebo
Design    Phase I/II Open-Label Dose-Escalation + Randomized Basket Trial
Patients  80  (Phase I = 20 TRT  ·  Phase II = 40 TRT + 20 CTL  ·  2:1 ratio)
Records   26,723 across 13 ADaM domains
Seed      42  ·  fully reproducible
Cutoff    March 5, 2026
```

---

## Study Design

This dataset models a two-part oncology basket trial across five solid tumor histologies.

**Phase I — Dose Escalation (n=20, Treatment only)**
Four dose levels (100 → 200 → 300 → 400 mg QD), five patients per cohort.
Enables: DLT plots, Dose Escalation (3+3 / BOIN), Exposure-Response, Trough Level.

**Phase II — Randomized Expansion (n=60)**
RP2D = 300 mg QD vs Placebo, 2:1 allocation (40 TRT : 20 CTL).
Enables: all response assessment, survival, subgroup, biomarker, and QoL charts.

---

## Why this dataset

There is a real gap in the oncology data visualization ecosystem. Most published visualization
catalogs use toy datasets too simple to reveal anything interesting, or stitch together figures
from different studies with incompatible assumptions. Reproducibility suffers, and readers cannot
verify whether a spider plot was generated from the same patients as the forest plot next to it.

Real patient-level RCT data is almost never publicly available. The pharmaverse `admiral` test
datasets are designed for software testing, not visual demonstration, and lack PK profiles, dose
modification histories, mutation records, and immune biomarkers. cBioPortal and TCGA provide rich
genomic data but no longitudinal tumor measurements, no QoL, no pharmacokinetics, and no
cross-domain consistency.

Synthetic data is the right answer — provided it is done carefully. The goal is a dataset that
behaves statistically like real oncology data, respects biological relationships between domains,
and is anchored to published benchmarks so that every number can be traced to a source.

A core design decision was **tumor stratification**: each of the five histologies has its own
calibrated ORR, survival parameters, mutation prevalence table, toxicity modifier, and mutational
signature profile. A KRAS mutation in PDAC (93%) means something completely different from one
in BRCA (2%). Most synthetic datasets use a single parameter set across all histologies, which
makes per-tumor subgroup analysis meaningless. This one does not.

---

## Dataset Inventory

| Domain   | Description                                   |   Rows | Cols | Notes |
|----------|-----------------------------------------------|-------:|-----:|-------|
| ADSL     | Subject-level analysis dataset                |     80 |   65 | PHASE, DOSELEVEL, COMPTYPE |
| ADRS     | Tumor response per RECIST 1.1                 |    769 |   14 | Longitudinal — one row per visit |
| ADTR     | Sum of longest diameters (mm)                 |    769 |   15 | Aligned with ADRS visit structure |
| ADAE     | Adverse events · MedDRA / CTCAE v5            |    752 |   19 | DOSELEVEL embedded |
| ADLB     | Laboratory + immune cell parameters · 20 tests| 11,880 |   21 | CD4/CD8/NK/TREG included |
| ADTTE    | Time-to-event · 7 endpoints                   |    377 |   33 | OS/PFS/EFS/TTP/DOR/TTR/DFS |
| ADPK     | Pharmacokinetics · 1-compartment model        |  1,984 |   17 | Treatment arm only · 10 timepoints |
| ADEX     | Dose exposure and modifications               |  2,361 |   20 | DOSEINT on 0–1 scale (RDI) |
| ADBM     | Biomarkers · ctDNA · immune panel             |  3,321 |   17 | Serial ctDNA with BOR-linked decay |
| ADPR     | Patient-reported outcomes · EORTC QLQ-C30     |  3,888 |   19 | 8 scales · MID=10 |
| ADMUT    | Somatic mutations · 15 genes                  |    187 |   23 | CNV records included |
| ADRAND   | Screening and randomization                   |     92 |   12 | 92 screened → 80 enrolled |
| ADSIG    | Mutational signatures · SBS                   |    263 |   13 | 12 unique SBS signatures |
| **Total**|                                               |**26,723**| | |

---

## Calibration Strategy

### Literature-based calibration

Each tumor type was anchored to a specific published trial or database:

| Histology / Domain | Primary source | Parameters borrowed |
|--------------------|---------------|---------------------|
| NSCLC | KEYNOTE-189 (Gandhi et al., *NEJM* 2018) | ORR, OS/PFS HR, AE incidence |
| HCC | IMbrave150 (Finn et al., *NEJM* 2020) | ORR, OS HR, PFS HR |
| CRC | KEYNOTE-177 (André et al., *NEJM* 2020) | ORR, survival |
| BRCA | OlympiAD (Robson et al., *NEJM* 2017) | ORR, survival |
| PDAC | NAPOLI-1 / PRODIGE-4 | ORR, survival |
| All histologies | TCGA PanCancer Atlas (cBioPortal) | Mutation prevalences |
| All histologies | COSMIC v3.3 SBS catalogue | Signature profiles |
| PK | 1-compartment oral model (Ling et al., 2006) | CL/F=15, Vd/F=100, Ka=1.5 |
| Safety grading | NCI CTCAE v5.0 | All lab threshold values |
| QoL | EORTC QLQ-C30 manual (3rd ed.) | Baseline norms, MID=10 |

---

## Tumor-Stratified Parameters

### Response and survival (Treatment arm)

| Histology | n TRT | n CTL | ORR TRT | ORR CTL | OS median TRT |
|-----------|------:|------:|--------:|--------:|--------------:|
| NSCLC     |    23 |     2 |     52% |      0% |        18.3 m |
| BRCA      |     8 |     5 |     38% |     20% |        16.7 m |
| HCC       |    10 |     2 |     10% |      — |         4.0 m |
| CRC       |     9 |     7 |     44% |     14% |        13.5 m |
| PDAC      |    12 |     2 |     17% |      0% |         3.5 m |
| **All**   |**62** |**18** | **33%** |         |               |

BOR distribution (all 80 patients): CR=16 (20%) · PR=10 (12%) · SD=16 (20%) · PD=38 (48%)
Responders (CR+PR): **26 patients**

### Somatic mutation prevalence · TCGA PanCancer Atlas

| Gene    | NSCLC | CRC | HCC | PDAC | BRCA | Key biology |
|---------|------:|----:|----:|-----:|-----:|-------------|
| TP53    |   46% | 62% | 28% |  72% |  38% | Universal tumour suppressor |
| KRAS    |   30% | 45% |  3% | **93%** | 3% | Dominant driver in PDAC |
| EGFR    |   15% |  2% |  1% |   2% |   2% | Targetable in NSCLC |
| PIK3CA  |    7% | 20% |  5% |   3% |  35% | PI3K pathway, high in BRCA |
| SMAD4   |    2% | 32% |  3% |  22% |   2% | TGF-β · CRC and PDAC |
| CDKN2A  |   12% |  5% |  8% |  30% |  15% | Cell cycle regulation |
| STK11   |   17% |  3% |  2% |   4% |   5% | LKB1 pathway · NSCLC |

### Toxicity modifiers by histology

KEYNOTE-189 base incidence rates multiplied by the factors below:

| SOC category     | NSCLC | HCC    | CRC    | BRCA   | PDAC   |
|------------------|------:|-------:|-------:|-------:|-------:|
| Hepatic          |  1.0× | **2.2×** | 0.8× | 0.9×  | 1.2×   |
| Gastrointestinal |  1.0× | 1.1×   | **1.8×** | 1.2× | **2.0×** |
| Dermatologic     |  1.0× | 0.9×   | 1.2×   | 1.3×   | 0.7×   |
| Hematologic      |  1.0× | 0.9×   | 1.1×   | **1.4×** | 1.3× |

### Mutational signatures · COSMIC v3.3 SBS

| Histology | Dominant signatures | Weight | Biology |
|-----------|--------------------:|-------:|---------|
| NSCLC     | SBS4                |    35% | Tobacco carcinogens |
| NSCLC     | SBS2 + SBS13        |    35% | APOBEC cytidine deaminase |
| CRC       | SBS15 + SBS6 + SBS44|    45% | Defective mismatch repair |
| HCC       | SBS22 + SBS24       |    40% | Aristolochic acid / aflatoxin |
| PDAC      | SBS1 + SBS5         |    65% | Age-related clock-like mutagenesis |
| BRCA      | SBS3                |    30% | HR deficiency / BRCA1–2 loss |
| BRCA      | SBS2 + SBS13        |    45% | APOBEC activity |

12 unique SBS signatures · weights sum to 1.000 per patient (validated)

---

## Domain Design Notes

### ADRS — Longitudinal response

One row per patient per assessment cycle (every 42 days), generated by a tumour-stratified
Markov chain. Median 13 visits per patient. This enables Spider Plots (≥3 on-study assessments
required), BOR derivation from ADRS, and response trajectory analyses.

BOR in ADSL (`BESTRSPC`) is derived from ADRS post-generation using the best observed response
across all visits, ensuring full cross-domain consistency. ADRS and ADSL BOR agree for all
80 patients (0 mismatches).

### ADTTE — Seven time-to-event endpoints

```
OS    Overall Survival              80 patients   events=44   censored=36
PFS   Progression-Free Survival     80 patients   events=77   censored=3
EFS   Event-Free Survival           80 patients   (≈ PFS for solid tumours)
TTP   Time to Progression           80 patients   (censors non-progression deaths)
DOR   Duration of Response          26 patients   responders only
TTR   Time to Response              26 patients   responders only
DFS   Disease-Free Survival          5 patients   CR patients only
```

Each record embeds 15 subgroup variables from ADSL — forest plots require no joins.
Three landmark flags pre-computed: `LM6MFL`, `LM12MFL`, `LM24MFL`.

Competing events in `COMPTYPE`:
```
PROGRESSION                  76 patients
CENSORED                      3 patients
DEATH_WITHOUT_PROGRESSION     1 patient   non-cancer death
```

### ADEX — Dose intensity

`DOSEINT` is stored as a proportion (0.0–1.0) — the correct scale for Relative Dose Intensity
(RDI) plots. Mean RDI across treatment cycles: **96.5%**, consistent with an early-phase
study with limited long-term dose modifications.

Dose modification logic per cycle (after cycle 2, treatment arm only):
```
Dose reduction     8% probability  →  level drops one step (300→200→100 mg)
Dose interruption  3% probability  →  one-cycle hold then automatic resumption
```

### ADLB — Immune cell markers

Four immune cell parameters alongside standard chemistry and haematology:
`CD4`, `CD8`, `NK`, `TREG` — enabling Immune Cell Infiltration Heatmaps and longitudinal
immune panel plots without a separate immunophenotyping dataset.

### ADBM — Serial ctDNA

`PARAMCD = "CTDNA"` measured at seven timepoints per patient with response-linked kinetics:
exponential decay for CR/PR, linear increase for PD. Enables ctDNA Dynamics plots and
response correlation analyses.

### ADPK — Pharmacokinetics

Treatment arm only (62 patients). Two PK cycles (Cycle 1 Day 1, Cycle 3 Day 1).
10 nominal timepoints: 0, 0.5, 1, 2, 3, 4, 6, 8, 12, 24 hours.
Summary parameters: `CMAX`, `AUC`, `AUCINF`, `TMAX`, `THALF`, `TROUGH`.

```
Cmax   median=1.9 ng/mL   range=0.5–3.4 ng/mL
AUCinf median=13.5 ng·h/mL
```

*Note: PK values are generated from a simplified 1-compartment oral model for structural
validity (shape, variability, dose proportionality). Absolute concentrations are not
calibrated to a specific drug and should not be interpreted as clinically meaningful.*

### ADRAND — Screening

```
Screened     92
Randomized   80
Screen-fail  12  (6 reason categories)
```

Enables a complete CONSORT flow diagram with reason-level breakdown.

---

## Data Validation

A dedicated validation script (`validate_oncology_data.py`) checks all 13 ADaM datasets
against chart-specific requirements before any visualization is produced. Running this
script against the released files confirms the following:

### Programmatic validation results

| Check | Result | Benchmark / Requirement |
|-------|--------|-------------------------|
| OS ≥ PFS, all patients | **0 violations** | Required |
| ADSL ↔ ADRAND exact ID match | **80 / 80** | Required |
| ADSL BOR consistent with ADRS | **0 mismatches** | Required |
| ADRS median visits per patient | **13** | ≥3 required for Spider Plot |
| ADTTE parameters present | **OS/PFS/EFS/TTP/DOR/TTR/DFS** | Required |
| DOSEINT range | **0.000 – 1.000** | Required for RDI plot |
| AVAL (survival) all non-negative | **0 violations** | Required |
| ADMUT gene count | **15 genes** | Required |
| ADMUT patients with sequencing | **72 / 80 (90%)** | Required |
| ADSIG weight sum per patient | **1.000 (all 80)** | Required |
| ADSIG unique signatures | **12** | ≥3 required for signature plot |
| Immune markers in ADLB | **CD4/CD8/NK/TREG** | Required |
| ctDNA timepoints per patient | **7** | ≥4 required for dynamics plot |
| Cross-domain orphan IDs | **0 in any domain** | Required |
| Grade ≥ 3 AE rate (treatment) | **19.5%** | 15–20% · KEYNOTE-189 |
| KRAS prevalence in PDAC | **93%** | ~90% · TCGA |
| TP53 prevalence in CRC | **62%** | ~60–65% · TCGA |
| Competing event categories | **3 present** | Required for CIF plots |

### How to run validation

```python
# Google Colab — run after uploading CSVs to /content/data/
# Change DATA_DIR at top of script if needed

# Expected output:
# ✅ PASS=10  ⚠️ WARN=0  ❌ FAIL=0
```

```bash
# Local
python validate_oncology_data.py
```

The validation script checks each dataset independently, then performs cross-dataset
consistency checks (BOR reconciliation, patient ID matching, BASESZ alignment between
ADSL and ADTR), and finally maps all 70+ chart types to their required columns to flag
any gaps before chart generation begins.

---

## Visualization Coverage

### Fully supported

| Category | Chart types | Primary domain |
|----------|-------------|----------------|
| Response Assessment | Waterfall · Spider · Swimmer · BOR · SLD over time · Tumor Burden | ADTR · ADRS |
| Survival | KM (OS/PFS/EFS/TTP/DOR/TTR/DFS) · Landmark · CIF · RMST | ADTTE |
| Genomics | OncoPrint · VAF · TMB · MSI · ctDNA dynamics · Mutational Signature | ADMUT · ADSIG · ADBM |
| Safety | AE bar · Toxicity heatmap · DLT · Dose Escalation · Exposure-Response | ADAE · ADLB |
| PK/PD | Concentration-time · Trough · Exposure-efficacy · Waterfall+PK overlay | ADPK · ADEX |
| Biomarker | Forest plot (15 subgroups) · Immune cell panel · PD-L1 by response | ADTTE · ADBM · ADLB |
| QoL / PRO | PRO trajectories · MID responder rate · Deterioration-free | ADPR |
| Trial Design | CONSORT · Enrollment curve · Dose intensity · RDI | ADRAND · ADEX |

### Not supported (structural incompatibility)

```
CNV / Circos / Manhattan / Rainfall     whole-genome data  →  TCGA
Flow cytometry / UMAP / t-SNE / CyTOF   single-cell data   →  GEO
CAR-T expansion / CRS timeline           CAR-T trial data   →  specialized sources
Radiomics overlays                       imaging data       →  TCIA
ASR / incidence trend lines              registry data      →  SEER / GLOBOCAN
Funnel plot / NMA                        multi-study data   →  meta-analysis datasets
```

---

## Usage

**Generate data**

```r
source("Data/generate_adam_oncviz001.R")

datasets <- generate_adam_oncviz001(
  output_dir    = "./Data/V1",
  seed          = 42L,
  n_phase1      = 20L,
  n_phase2_trt  = 40L,
  n_phase2_ctl  = 20L,
  verbose       = TRUE
)
# Output: 13 CSV files in Data/V1/  +  Data/V1.zip
```

**Validate data**

```python
# Python / Google Colab
# Set DATA_DIR = "./Data/V1" at top of script
python validate_oncology_data.py
```

**Load and subset**

```r
library(dplyr)

adsl  <- read.csv("Data/V1/ADSL.csv")
adtte <- read.csv("Data/V1/ADTTE.csv")
adrs  <- read.csv("Data/V1/ADRS.csv")
adtr  <- read.csv("Data/V1/ADTR.csv")

# Phase II Treatment arm · NSCLC only
nsclc_trt <- adsl |>
  filter(TUMORTYPE == "NSCLC", ARM == "TREATMENT", PHASE == "II")

# OS with 15 subgroup variables embedded — no join required
os <- adtte |> filter(PARAMCD == "OS")
# Subgroups: TUMORTYPE, AGEGR1, SEX, ECOG, PDL1GRP, MSISTS, TMBHIGH,
#            LIVERMETS, PRIORLINES, EGFRMUT, KRASMUT, TP53MUT, STK11MUT, STAGE

# Waterfall: best % change per patient from ADTR
best_pchg <- adtr |>
  filter(AVISITN > 0) |>
  group_by(USUBJID) |>
  slice_min(PCHG, n = 1) |>
  ungroup()

# RDI: DOSEINT already on 0–1 scale
adex <- read.csv("Data/V1/ADEX.csv")
rdi <- adex |>
  filter(ARM == "TREATMENT") |>
  group_by(USUBJID) |>
  summarise(mean_rdi = mean(DOSEINT), .groups = "drop")

# ctDNA dynamics
adbm <- read.csv("Data/V1/ADBM.csv")
ctdna <- adbm |> filter(PARAMCD == "CTDNA") |>
  left_join(adsl |> select(USUBJID, BESTRSPC), by = "USUBJID")
```

---

## Repository Structure

```
Data/
├── generate_adam_oncviz001.R    ← data generation script
└── V1/
    ├── ADSL.csv
    ├── ADRS.csv
    ├── ADTR.csv
    ├── ADAE.csv
    ├── ADLB.csv
    ├── ADTTE.csv
    ├── ADPK.csv
    ├── ADEX.csv
    ├── ADBM.csv
    ├── ADPR.csv
    ├── ADMUT.csv
    ├── ADRAND.csv
    └── ADSIG.csv

validate_oncology_data.py        ← data validation script
```

---

## Reproducibility

Fixed seed `set.seed(42)` applied before all sampling.
Produces bitwise-identical output on:

```
R      >= 4.1
dplyr  >= 1.0
tidyr  >= 1.1
purrr  >= 0.3
```

License: **CC BY 4.0** — unrestricted reuse with attribution.

---

## References

André T, Shiu K-K, Kim TW, et al. Pembrolizumab in microsatellite-instability–high advanced
colorectal cancer. *N Engl J Med.* 2020;383(23):2207–2218.

Cancer Genome Atlas Research Network. Comprehensive molecular profiling of lung adenocarcinoma.
*Nature.* 2014;511(7511):543–550.

Cerami E, Gao J, Dogrusoz U, et al. The cBio cancer genomics portal. *Cancer Discov.*
2012;2(5):401–404.

Eisenhauer EA, Therasse P, Bogaerts J, et al. New response evaluation criteria in solid tumours:
RECIST 1.1. *Eur J Cancer.* 2009;45(2):228–247.

Fayers PM, Aaronson NK, Bjordal K, et al. *The EORTC QLQ-C30 Scoring Manual.* 3rd ed.
Brussels: EORTC; 2001.

Finn RS, Qin S, Ikeda M, et al. Atezolizumab plus bevacizumab in unresectable hepatocellular
carcinoma. *N Engl J Med.* 2020;382(20):1894–1905.

Gandhi L, Rodríguez-Abreu D, Gadgeel S, et al. Pembrolizumab plus chemotherapy in metastatic
non-small-cell lung cancer. *N Engl J Med.* 2018;378(22):2078–2092.

Gao J, Aksoy BA, Dogrusoz U, et al. Integrative analysis of complex cancer genomics using
cBioPortal. *Sci Signal.* 2013;6(269):pl1.

Ling J, Johnson KA, Miao Z, et al. Metabolism and excretion of erlotinib in healthy male
volunteers. *Drug Metab Dispos.* 2006;34(3):420–426.

National Cancer Institute. *Common Terminology Criteria for Adverse Events (CTCAE) v5.0.*
Bethesda: NIH; 2017.

Robson M, Im S-A, Senkus E, et al. Olaparib for metastatic breast cancer in patients with a
germline BRCA mutation. *N Engl J Med.* 2017;377(6):523–533.

Tate JG, Bamford S, Jubb HC, et al. COSMIC: the catalogue of somatic mutations in cancer.
*Nucleic Acids Res.* 2019;47(D1):D941–D947.

TCGA Research Network. Comprehensive genomic characterization of squamous cell lung cancers.
*Nature.* 2012;489(7417):519–525.
