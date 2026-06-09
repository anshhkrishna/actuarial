# NHANES Epigenetic Aging Clock Validation

Survival analysis validation of epigenetic aging clocks (GrimAge, GrimAge2, Horvath, Hannum, PhenoAge, DunedinPoAm) against mortality outcomes in NHANES.

## Project Status

**Phase: Analytic Cohort Construction Complete ✓**

- [x] Data acquisition (DNAm, demographics, mortality linkage)
- [x] Data merging and integration
- [x] Cohort construction with explicit filtering
- [x] Survival variables construction
- [x] Predictor standardization (z-scoring)
- [x] Cohort descriptive statistics
- [x] Train/test split (stratified, 70/30)
- [ ] Cox proportional hazards modeling
- [ ] Model validation and performance metrics

## Project Structure

```
actuarial-validation/
├── data/
│   ├── raw/           # NHANES downloaded files
│   │   ├── dnmepi.sas7bdat (1.19 MB, 4,449 rows)
│   │   ├── DEMO_1999.xpt (10.97 MB, 9,965 rows)
│   │   ├── DEMO_2001.xpt (3.12 MB, 11,039 rows)
│   │   ├── NHANES_1999_2000_MORT_2019_PUBLIC.dat (476 KB)
│   │   └── NHANES_2001_2002_MORT_2019_PUBLIC.dat (528 KB)
│   └── processed/     # Merged analytic datasets
│       ├── analytic_cohort.csv (2,532 rows × 29 columns)
│       └── train_test_split.json
├── src/
│   ├── __init__.py
│   ├── download.py    # Data acquisition & documentation
│   ├── verify_variables.py (verification script)
│   └── build_cohort.py (cohort construction)
├── notebooks/         # Analysis notebooks (coming)
├── figures/           # Output plots (coming)
├── requirements.txt   # Python dependencies
└── README.md          # This file
```

## Data Sources (Verified URLs & Specifications)

### NHANES DNA Methylation (DNAm) Epigenetic Biomarkers
- **File:** `dnmepi.sas7bdat` (combined 1999-2002 cycles)
- **URL:** https://wwwn.cdc.gov/nchs/data/nhanes/dnam/dnmepi.sas7bdat
- **Size:** 1.19 MB
- **Rows:** 4,449 (ages 50+)
- **Array:** Illumina MethylationEPIC v1.0 (~850K CpG sites)
- **Format:** SAS7BDAT binary

### NHANES Demographics
- **1999-2000:** `DEMO.xpt` (9,965 rows)
  - URL: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/1999/DataFiles/DEMO.xpt
- **2001-2002:** `DEMO_B.xpt` (11,039 rows)
  - URL: https://wwwn.cdc.gov/Nchs/Data/Nhanes/Public/2001/DataFiles/DEMO_B.xpt
- **Format:** SAS Transport (.xpt)
- **Key variables:** SEQN, RIDAGEYR, RIAGENDR, SDMVPSU, SDMVSTRA

### Linked Mortality (National Death Index, NDI)
- **1999-2000:** https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/NHANES_1999_2000_MORT_2019_PUBLIC.dat
- **2001-2002:** https://ftp.cdc.gov/pub/Health_Statistics/NCHS/datalinkage/linked_mortality/NHANES_2001_2002_MORT_2019_PUBLIC.dat
- **Format:** Fixed-width ASCII
- **Linkage:** Through 2019
- **Key variables:** SEQN, ELIGSTAT, MORTSTAT, PERMTH_INT, PERMTH_EXM (person-months)

## Epigenetic Clock Variables (Empirically Verified)

All variable names verified by loading actual SAS7BDAT file and inspecting column names + labels.

| Clock | **Column Name** | Description | r with Age |
|-------|---------|-------------|-----------|
| **Horvath** | `HorvathAge` | Original epigenetic age (51 tissues) | 0.793 |
| **Hannum** | `HannumAge` | Epigenetic age (whole blood) | 0.808 |
| **Skin/Blood** | `SkinBloodAge` | Horvath age (skin & blood tissues) | 0.850 |
| **PhenoAge** | `PhenoAge` | Phenotypic age (Levine) | 0.762 |
| **GrimAge** | `GrimAgeMort` | Mortality predictor (Horvath 2018) | 0.840 |
| **GrimAge2** | `GrimAge2Mort` | Updated mortality predictor (2024) | 0.804 |
| **DunedinPoAm** | `DunedinPoAm` | Pace-of-aging (biological tempo) | **0.036** ✓ |
| **DNAmTL** | `HorvathTelo` | Telomere length (Horvath) | **-0.582** ✓ |

**Key Identifiers & Weights:**
- Participant ID: `SEQN`
- Survey weight (DNAm subsample): `WTDN4YR` (4-year combined weight, not WTINT2YR)

## Analytic Cohort Summary

**Final sample: n=2,532** (eligible NHANES participants aged 50+ with DNAm and mortality linkage)

### Filtering steps:
- DNAm subsample: 4,449 → 2,532 (57% have GrimAgeMort phenotype data)
- Merged with demographics: 4,449 → 4,449 (100% match)
- Eligible (ELIGSTAT=1): 4,446 → 4,446 (dropped 3 ineligible)
- GrimAgeMort non-null: 4,446 → 2,532 (dropped 1,914 without biomarker)
- Positive follow-up: 2,532 → 2,532 (no drops)

### Demographics
- **Chronological age:** 66.1 ± 10.1 years (range 50–85, top-coded at 85)
- **Female:** 49.2%
- **Follow-up time:** Median 17.1 years, max 20.8 years

### Mortality Outcomes
- **Events (deaths):** 1,361 (53.8%)
- **Censored (alive):** 1,171 (46.2%)

### Mortality Validation ✓
Deceased participants older on both measures:
- **Chronological age:** Dead 71.1 vs. Alive 60.3 years
- **GrimAgeMort:** Dead 71.1 vs. Alive 61.0

### Clock Characteristics
- **Age-trained clocks:** r(RIDAGEYR) = 0.76–0.85 (Horvath, Hannum, etc.)
- **Pace-of-aging (DunedinPoAm):** r(RIDAGEYR) = 0.036 ✓ (appropriately uncorrelated)
- **Telomere length (HorvathTelo):** r(RIDAGEYR) = −0.58 ✓ (inverse, as expected)

### Standardized Predictors
- Z-scored versions (mean 0, SD 1) created for all clocks + chronological age
- Variable names: `z_HorvathAge`, `z_GrimAgeMort`, etc.

### Train/Test Split (Stratified on Event)
- **Train:** 1,772 rows, 952 events (53.7%)
- **Test:** 760 rows, 409 events (53.8%)
- Random state: 42

## Setup & Usage

### Installation
```bash
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

### Download Data
```bash
python src/download.py  # Downloads + validates NHANES data
```

### Build Analytic Cohort
```bash
python src/build_cohort.py  # Merges data, constructs cohort, saves to CSV
```

### Output Files
- **`data/processed/analytic_cohort.csv`** — Full cohort (2,532 rows × 29 columns)
  - Raw predictors: RIDAGEYR, 8 epigenetic clocks, survey design variables
  - Z-scored predictors: z_RIDAGEYR, z_HorvathAge, z_GrimAgeMort, etc.
  - Survival variables: time_years, event (0=censored, 1=dead)
  - Metadata: SEQN, WTDN4YR, MORTSTAT, PERMTH_EXM, etc.

- **`data/processed/train_test_split.json`** — Train/test indices
  - Stratified on event status
  - 70/30 split, random_state=42

## Analysis Pipeline (Next Steps)

1. **Cox Proportional Hazards Models** (on train set)
   - Univariate: each clock vs. mortality
   - Adjusted: each clock + chronological age + sex
   - Full: all clocks + demographics + survey adjustment

2. **Model Validation** (on test set)
   - Concordance index (C-index)
   - Calibration (observed vs. predicted)
   - Comparative performance (clock vs. age alone)

3. **Visualization**
   - Kaplan-Meier curves (stratified by clock quartiles)
   - Forest plots (hazard ratios + 95% CI)
   - ROC curves, performance comparisons

## Survey Design Considerations

- **Stratification:** SDMVSTRA (strata)
- **Clustering:** SDMVPSU (primary sampling units)
- **Sample weight:** WTDN4YR (4-year combined DNAm weight)
- Recommend using `lifelines` or `statsmodels` with survey weights for weighted Cox models

## Known Constraints & Data Characteristics

1. **Survey Design**
   - Complex survey: stratified, multistage cluster sampling
   - Use WTDN4YR (not WTINT2YR) for DNAm subsample analyses
   - Must account for SDMVPSU (clustering) and SDMVSTRA (stratification)

2. **Mortality Linkage**
   - NDI linkage through 2019 (20-year follow-up)
   - Public-use file only (restricted-access data not available)
   - ELIGSTAT coding: 1=eligible, 2=under 18, 3=ineligible
   - Only ELIGSTAT=1 participants in analytic cohort

3. **Sample Composition**
   - Older adult cohort (aged 50–85)
   - High event rate (53.8%, appropriate for ~20-year follow-up)
   - DNAm subsample: 2,532 / 4,449 (57% have GrimAgeMort data)
   - ~46% female, ~54% male

4. **DNAm Array & Phenotypes**
   - Illumina MethylationEPIC v1.0 (~850K sites)
   - Includes 8 epigenetic clocks + immune cell proportions + biomarkers
   - No CpG-level matrix in public-use file (restricted access)

## References

- Levine, M. E., et al. (2018). An epigenetic biomarker of aging for lifespan and healthspan. *Aging Cell*, 17(4), e12759.
- Horvath, S., & Raj, K. (2018). DNA methylation-based biomarkers and the epigenetic clock theory of aging. *Trends in Genetics*, 34(11), 856–862.
- NHANES Survey Methods and Analytic Guidelines: https://www.cdc.gov/nchs/nhanes/analyticguidelines.aspx
- NCHS Data Linkage – Linked Mortality Files: https://www.cdc.gov/nchs/data-linkage/mortality-public.htm
- CDC DNAm Epigenetic Biomarkers: https://wwwn.cdc.gov/nchs/nhanes/dnam/
