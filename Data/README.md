# ONCVIZ-001 · ADaM v2
### A Tumor-Stratified Synthetic Clinical Trial Dataset for Oncology Visualization

```
Study     ONCVIZ-001  ·  Vizatinib 300 mg QD vs Placebo
Design    Phase I/II Open-Label Dose-Escalation + Randomized Basket Trial
Patients  80  (Phase I = 20 TRT escalation  ·  Phase II = 40 TRT + 20 CTL)
Records   26,723 across 13 ADaM domains
Seed      42  ·  fully reproducible
Cutoff    March 5, 2026
```

---

## Study Design

This dataset models a two-part oncology basket trial:

**Phase I — Dose Escalation (n=20, Treatment only)**
Four dose levels (100 → 200 → 300 → 400 mg QD), five patients per cohort.
Generates dose-level data for DLT plots, Dose Escalation (3+3 / BOIN), Exposure-Response,
and Trough Level analyses.

**Phase II — Randomized Expansion (n=60)**
Selected dose (RP2D = 300 mg QD) versus Placebo, 2:1 ratio (40 TRT : 20 CTL).
Generates all response assessment, survival, subgroup, and biomarker charts.

This design is intentional. A single Phase I/II basket trial covering five solid tumor histologies
is one of the most common structures in modern oncology publications, and having both dose
escalation and randomized expansion data in a single consistent dataset enables charts that would
otherwise require two separate trials.

---

## Why this dataset

There is a real gap in the oncology data visualization ecosystem. Most published visualization
catalogs either use toy datasets that are too simple to reveal anything interesting, or they stitch
together figures from different studies with incompatible assumptions about trial design, patient
populations, and data structure. Reproducibility suffers. Readers cannot verify whether a spider
plot or a competing risks curve was generated from the same patients as the forest plot sitting
next to it.

The obvious alternative is real clinical data. The problem is that patient-level data from
randomized controlled trials is almost never publicly available at the resolution needed for this
kind of work. Regulatory submissions contain everything, but are not released at the individual
level. The pharmaverse `admiral` test datasets are designed for software testing, not visual
demonstration, and lack PK profiles, dose modification histories, mutation records, and immune
biomarkers. cBioPortal and TCGA provide rich genomic data but no longitudinal tumor measurements,
no QoL, no pharmacokinetics, and no cross-domain consistency.

Synthetic data is the right answer here, provided it is done carefully. The goal is not to trick
anyone into thinking this is real. The goal is a dataset that behaves statistically like real
oncology data, respects the biological relationships between domains, and is anchored to published
benchmarks so that every number can be traced to a source.

A core design decision was tumor stratification: each of the five histologies has its own
calibrated ORR, survival parameters, mutation prevalence table, and mutational signature profile.
A KRAS mutation in PDAC (93% in this dataset) means something completely different from one
in BRCA (2%). Most synthetic datasets use a single parameter set across all histologies, which
makes per-tumor subgroup analysis meaningless. This one does not.

---

## Dataset Inventory

| Domain   | Description                                   |  Rows | Cols | Notes |
|----------|-----------------------------------------------|------:|-----:|-------|
| ADSL     | Subject-level analysis dataset                |    80 |   65 | Phase I/II flag, DOSELEVEL, COMPTYPE |
| ADRS     | Tumor response per RECIST 1.1                 |   769 |   14 | Longitudinal — one row per visit per patient |
| ADTR     | Sum of longest diameters (mm)                 |   769 |   15 | Aligned visit structure with ADRS |
| ADAE     | Adverse events · MedDRA / CTCAE v5            |   752 |   19 | DOSELEVEL from ADEX join |
| ADLB     | Laboratory + immune cell parameters · 20 tests| 11,880|   21 | CD4/CD8/NK/TREG included |
| ADTTE    | Time-to-event · 7 endpoints                   |   377 |   33 | OS/PFS/EFS/TTP/DOR/TTR/DFS |
| ADPK     | Pharmacokinetics · 1-compartment model        | 1,984 |   17 | Treatment arm only |
| ADEX     | Dose exposure and modifications               | 2,361 |   20 | DOSEINT on 0–1 scale (RDI) |
| ADBM     | Biomarkers, ctDNA, immune panel               | 3,321 |   17 | Serial ctDNA with decay by BOR |
| ADPR     | Patient-reported outcomes · EORTC QLQ-C30     | 3,888 |   19 | — |
| ADMUT    | Somatic mutations · 15 genes                  |   187 |   23 | CNV records included |
| ADRAND   | Screening and randomization                   |    92 |   12 | 92 screened → 80 enrolled |
| ADSIG    | Mutational signatures · SBS                   |   263 |   13 | 12 unique SBS signatures |
| **Total**|                                               |**26,723**| | |

---

## Key Design Decisions versus v1

| Feature | v1 (400 patients) | v2 (80 patients) |
|---------|-------------------|------------------|
| N | 400 (Phase II/III) | 80 (Phase I/II) |
| Study design | Two-arm only | Phase I escalation + Phase II randomized |
| ADRS structure | One row per patient (BOR only) | **Longitudinal** — one row per visit per patient |
| BOR consistency | ADSL and ADRS generated independently | **BOR in ADSL reconciled from ADRS** after generation |
| DOSEINT scale | 0–100 (percent) | **0–1 (proportion)** — correct RDI scale |
| EFS / DFS / TTP | Absent | **Added** — all seven TTE endpoints present |
| Immune markers | Absent from ADLB | **CD4 / CD8 / NK / TREG** in ADLB |
| ctDNA | Static single visit | **Serial measurements** with response-linked decay |
| Dose levels | Single dose | **Four dose levels** (100/200/300/400 mg) |
| ADTTE subgroup cols | 15 embedded | 15 embedded (unchanged) |

---

## Tumor-Stratified Parameters

### Response and survival (Treatment arm)

| Histology | n (TRT) | n (CTL) | ORR TRT | ORR CTL | OS median TRT | ADTTE n |
|-----------|--------:|--------:|--------:|--------:|--------------:|--------:|
| NSCLC | 23 | 2 | 52% | 0% | 18.3 m | — |
| BRCA | 8 | 5 | 38% | 20% | 16.7 m | — |
| HCC | 10 | 2 | 10% | — | 4.0 m | — |
| CRC | 9 | 7 | 44% | 14% | 13.5 m | — |
| PDAC | 12 | 2 | 17% | 0% | 3.5 m | — |
| **All** | **62** | **18** | **33%** | — | — | **80** |

*Note: Small per-tumor control arm sizes reflect the 2:1 randomisation with 80 total patients.
Control arm subgroup analyses should be interpreted with caution.*

### BOR distribution (all patients)

```
CR   16 patients  (20%)
PR   10 patients  (12%)
SD   16 patients  (20%)
PD   38 patients  (48%)

Responders (CR+PR): 26 patients (33%)
```

### ADTTE endpoint coverage

```
OS      80 records   all patients
PFS     80 records   all patients
EFS     80 records   all patients (≈ PFS for solid tumors)
TTP     80 records   all patients (censors non-progression deaths)
DOR     26 records   responders only (CR+PR)
TTR     26 records   responders only
DFS      5 records   CR patients only
```

### Somatic mutation prevalence (this dataset vs TCGA benchmark)

| Gene | NSCLC | CRC | HCC | PDAC | BRCA | TCGA benchmark |
|------|------:|----:|----:|-----:|-----:|----------------|
| TP53 | ~46% | **62%** | ~28% | ~72% | ~38% | 49/68/32/72/38% |
| KRAS | ~30% | ~45% | ~3% | **93%** | ~3% | 30/42/1/91/2% |
| EGFR | ~15% | ~2% | ~1% | ~2% | ~2% | 10/1/0/1/4% |
| PIK3CA | ~7% | ~20% | ~5% | ~3% | ~35% | 7/20/5/3/35% |
| SMAD4 | ~2% | ~32% | ~3% | ~22% | ~2% | 2/32/3/22/2% |

*Small cohort (n=80) introduces sampling variance. Gene-level prevalence may deviate ±10–15%
from TCGA benchmarks within individual tumor type subgroups.*

### Mutational signatures (COSMIC v3.3 SBS)

| Histology | Dominant signatures | Biology |
|-----------|--------------------|---------| 
| NSCLC | SBS4 (35%), SBS2+SBS13 (35%) | Tobacco / APOBEC |
| CRC | SBS1+SBS5 (55%), SBS15+SBS6+SBS44 (45%) | Clock-like / MMR deficiency |
| HCC | SBS4+SBS22+SBS24 (65%) | Tobacco / aristolochic acid / aflatoxin |
| PDAC | SBS1+SBS5 (65%) | Age-related clock-like mutagenesis |
| BRCA | SBS3 (30%), SBS2+SBS13 (45%) | HR deficiency / APOBEC |

12 unique SBS signatures across all patients. Weights normalized to sum to 1.000 per patient (validated).

---

## Domain Design Notes

### ADRS — Longitudinal response (key change from v1)

ADRS now has the same visit structure as ADTR: one row per patient per assessment cycle,
generated by a Markov chain starting from the patient's initial response category.
This enables Spider Plots (which require ≥3 on-study response assessments per patient),
as well as waterfall plot construction and BOR derivation from ADRS without relying on ADSL.

BOR in ADSL (`BESTRSPC`) is derived from ADRS post-generation using the best observed response
across all visits, ensuring cross-domain consistency. The v1 mismatch (269 patients with
inconsistent BOR) is fully resolved.

### ADTTE — Seven time-to-event endpoints

```
OS   Overall Survival            all 80 patients
PFS  Progression-Free Survival   all 80 patients
EFS  Event-Free Survival         all 80 patients  (= PFS for solid tumors)
TTP  Time to Progression         all 80 patients  (censors non-progression deaths)
DOR  Duration of Response        26 responders
TTR  Time to Response            26 responders
DFS  Disease-Free Survival        5 CR patients
```

Each record embeds 15 subgroup variables from ADSL so forest plots can be built without
additional joins. Three landmark flags pre-computed: `LM6MFL`, `LM12MFL`, `LM24MFL`.

Competing events encoded in `COMPTYPE`:

```
PROGRESSION                   76 patients
DEATH_WITHOUT_PROGRESSION      1 patient   non-cancer death
CENSORED                       3 patients  administrative
```

### ADEX — Dose intensity (key change from v1)

`DOSEINT` is stored as a proportion (0.0 to 1.0), not as a percentage (0–100).
This is the correct scale for Relative Dose Intensity (RDI) plots and matches
the standard publication convention. Mean RDI across treatment cycles: **96.5%**,
reflecting the Phase I/II early-study population with fewer long-term dose modifications
compared to a Phase II/III dataset.

Dose modification logic:
```
Dose reduction      8% per-cycle probability after cycle 2
                    one level at a time: 300 → 200 → 100 mg
Dose interruption   3% per-cycle probability (one cycle, then resumption)
```

### ADLB — Immune cell markers (added in v2)

Four immune cell parameters added to ADLB alongside standard chemistry / haematology:

```
CD4    CD4+ T Cells          %
CD8    CD8+ T Cells          %
NK     NK Cells              %
TREG   Regulatory T Cells    %
```

These enable Immune Cell Infiltration Heatmaps and longitudinal immune panel plots
without requiring a separate immunophenotyping dataset.

### ADBM — Serial ctDNA dynamics (enhanced in v2)

ctDNA (`PARAMCD = "CTDNA"`) is measured at seven timepoints per patient with
response-linked kinetics: exponential decay for CR/PR patients, linear increase
for PD patients. This enables ctDNA Dynamics plots and response correlation analyses.

### ADPK — Pharmacokinetics

Treatment arm only (62 patients), two PK cycles (Cycle 1 Day 1 and Cycle 3 Day 1),
10 nominal timepoints: 0, 0.5, 1, 2, 3, 4, 6, 8, 12, 24 hours.

Summary PK parameters available:

```
CMAX    Cmax                62 patients
AUC     AUClast             62 patients
AUCINF  AUCinf              62 patients
TMAX    Time to Cmax        62 patients
THALF   Terminal half-life  62 patients
TROUGH  Trough concentration 62 patients
```

### ADRAND — Screening

```
Screened       92
Randomized     80
Screen-fail    12  (6 reason categories)
```
Enables a complete CONSORT flow diagram with reason-level breakdown.

---

## Calibration Sources

### Literature-based calibration

| Histology / Domain | Primary source | Parameters borrowed |
|--------------------|---------------|---------------------|
| NSCLC | KEYNOTE-189 (Gandhi et al., NEJM 2018) | ORR, OS/PFS HR, AE incidence |
| HCC | IMbrave150 (Finn et al., NEJM 2020) | ORR, OS HR, PFS HR |
| CRC | KEYNOTE-177 (André et al., NEJM 2020) | ORR, survival |
| BRCA | OlympiAD (Robson et al., NEJM 2017) | ORR, survival |
| PDAC | NAPOLI-1 / PRODIGE-4 | ORR, survival |
| All histologies | TCGA PanCancer Atlas (cBioPortal) | Mutation prevalences |
| All histologies | COSMIC v3.3 SBS catalogue | Signature profiles |
| PK | 1-compartment oral model (Ling et al., 2006) | CL/F=15, Vd/F=100, Ka=1.5 |
| Safety grading | NCI CTCAE v5.0 | All lab threshold values |
| QoL | EORTC QLQ-C30 manual (3rd ed.) | Baseline norms, MID=10 |

---

## Validation Summary

| Check | Result | Requirement |
|-------|--------|-------------|
| OS ≥ PFS, all patients | **0 violations** | Required |
| ADSL ↔ ADRAND exact ID match | **80 / 80** | Required |
| BOR: ADSL consistent with ADRS | **0 mismatches** | Required (fixed from v1) |
| ADLB null PARAMCD | **0** | Required |
| ADMUT gene count | **15 genes** | Required |
| ADSIG weight sum per patient | **1.000 (all 80)** | Required |
| Cross-domain extra IDs | **0 in any domain** | Required |
| ADTTE parameters present | **OS/PFS/EFS/TTP/DOR/TTR/DFS** | Required |
| DOSEINT range | **0.000 – 1.000** | Required (fixed from v1) |
| Grade ≥ 3 AE rate (treatment) | **19.5%** | 15–20% · KEYNOTE-189 |
| KRAS prevalence in PDAC | **93%** | ~90% · TCGA |
| TP53 prevalence in CRC | **62%** | ~60–65% · TCGA |
| ctDNA serial measurements | **7 timepoints/patient** | Required for dynamics plot |
| Immune markers in ADLB | **CD4/CD8/NK/TREG** | Required (added in v2) |
| Competing event categories | **3 present** | Required for CIF plots |

---

## Visualization Coverage

### Fully supported

| Category | Chart types | Primary domain |
|----------|-------------|----------------|
| Response Assessment | Waterfall, Spider, Swimmer, BOR, SLD over time | ADTR · ADRS |
| Survival | KM (OS/PFS/EFS/TTP/DOR/TTR/DFS), Landmark, CIF, RMST | ADTTE |
| Genomics | OncoPrint, VAF, TMB, MSI, ctDNA dynamics, Mutational Signature | ADMUT · ADSIG · ADBM |
| Safety | AE bar chart, Toxicity heatmap, DLT plot, Dose Escalation | ADAE · ADLB |
| PK/PD | Concentration-time, Trough, Exposure-efficacy, Waterfall+PK overlay | ADPK · ADEX |
| Biomarker | Forest plot (15 subgroups), Immune cell panel, PD-L1 by response | ADTTE · ADBM · ADLB |
| QoL / PRO | PRO trajectories, MID responder rate, Deterioration-free | ADPR |
| Trial Design | CONSORT, Enrollment curve, Dose intensity, RDI | ADRAND · ADEX |

### Not supported (structural incompatibility)

```
CNV / Circos / Manhattan / Rainfall    whole-genome data  →  TCGA
Flow cytometry / UMAP / t-SNE / CyTOF  single-cell data   →  GEO
CAR-T expansion / CRS timeline          CAR-T trial data   →  specialized sources
Radiomics overlays                      imaging data       →  TCIA
ASR / incidence trend lines             registry data      →  SEER / GLOBOCAN
Funnel plot / NMA                       multi-study data   →  meta-analysis datasets
```

---

## Usage

**Generate data**

```r
source("generate_adam_oncviz001.R")

datasets <- generate_adam_oncviz001(
  output_dir    = "./data_v2",
  seed          = 42L,
  n_phase1      = 20L,    # Phase I dose-escalation patients (Treatment only)
  n_phase2_trt  = 40L,    # Phase II Treatment arm
  n_phase2_ctl  = 20L,    # Phase II Control arm
  verbose       = TRUE
)
# Output: 13 CSV files + data_v2.zip
```

**Load and subset**

```r
library(dplyr)

adsl  <- read.csv("data_v2/ADSL.csv")
adtte <- read.csv("data_v2/ADTTE.csv")
adrs  <- read.csv("data_v2/ADRS.csv")
adtr  <- read.csv("data_v2/ADTR.csv")

# Phase II Treatment arm NSCLC only
nsclc_trt <- adsl |> filter(TUMORTYPE == "NSCLC", ARM == "TREATMENT", PHASE == "II")

# OS with all 15 subgroup variables embedded — no join required
os <- adtte |> filter(PARAMCD == "OS")
# Available subgroups: TUMORTYPE, AGEGR1, SEX, ECOG, PDL1GRP, MSISTS,
# TMBHIGH, LIVERMETS, PRIORLINES, EGFRMUT, KRASMUT,
# TP53MUT, STK11MUT, STAGE

# Longitudinal response — aligned with ADTR visits
adrs_trt <- adrs |>
  filter(ARM == "TREATMENT", AVISITN > 0, AVALC != "NE") |>
  left_join(adsl |> select(USUBJID, TUMORTYPE, BESTRSPC), by = "USUBJID")

# Waterfall: best % change from ADTR
best_pchg <- adtr |>
  filter(AVISITN > 0) |>
  group_by(USUBJID) |>
  slice_min(PCHG, n = 1) |>
  ungroup()

# RDI: DOSEINT is already on 0–1 scale
adex <- read.csv("data_v2/ADEX.csv")
rdi <- adex |>
  filter(ARM == "TREATMENT") |>
  group_by(USUBJID) |>
  summarise(mean_rdi = mean(DOSEINT), .groups = "drop")
```

---

## Reproducibility

The fixed seed `set.seed(42)` is applied before any sampling.
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
