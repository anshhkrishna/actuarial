"""
Build the merged analytic cohort for survival analysis validation.

Steps:
1. Load three data sources (DNAm, Demographics, Mortality)
2. Merge with explicit filtering
3. Construct survival variables
4. Standardize predictors
5. Print descriptive statistics
6. Create train/test split
7. Save cohort
"""

import pandas as pd
import numpy as np
from pathlib import Path
from sklearn.model_selection import train_test_split
from sklearn.preprocessing import StandardScaler

# Paths
RAW_DATA_DIR = Path("data/raw")
PROC_DATA_DIR = Path("data/processed")
PROC_DATA_DIR.mkdir(exist_ok=True)

print("\n" + "="*80)
print("STEP 1: LOAD THREE DATA SOURCES")
print("="*80)

# Load DNAm
print("\n[1/3] DNAm (dnmepi.sas7bdat)")
import pyreadstat
dnam_df, _ = pyreadstat.read_sas7bdat(str(RAW_DATA_DIR / "dnmepi.sas7bdat"))

# Select required columns
dnam_cols = ["SEQN", "WTDN4YR", "HorvathAge", "HannumAge", "SkinBloodAge",
             "PhenoAge", "GrimAgeMort", "GrimAge2Mort", "DunedinPoAm", "HorvathTelo"]
dnam_df = dnam_df[dnam_cols].copy()
print(f"  ✓ Loaded {dnam_df.shape[0]:,} rows × {dnam_df.shape[1]} columns")

# Load Demographics
print("\n[2/3] Demographics (1999-2000 and 2001-2002)")
demo_1999, _ = pyreadstat.read_xport(str(RAW_DATA_DIR / "DEMO_1999.xpt"))
demo_2001, _ = pyreadstat.read_xport(str(RAW_DATA_DIR / "DEMO_2001.xpt"))

demo_cols = ["SEQN", "RIDAGEYR", "RIAGENDR", "SDMVPSU", "SDMVSTRA"]
demo_1999 = demo_1999[demo_cols].copy()
demo_2001 = demo_2001[demo_cols].copy()

demo_df = pd.concat([demo_1999, demo_2001], ignore_index=True)
print(f"  ✓ Loaded {demo_df.shape[0]:,} rows (1999: {demo_1999.shape[0]:,}, 2001: {demo_2001.shape[0]:,})")

# Load Mortality (parse fixed-width)
print("\n[3/3] Linked Mortality (1999-2000 and 2001-2002)")

colspecs = [(0, 6), (14, 15), (15, 16), (16, 19), (19, 20), (20, 21), (42, 45), (45, 48)]
colnames = ["SEQN", "ELIGSTAT", "MORTSTAT", "UCOD_LEADING", "DIABETES", "HYPERTEN", "PERMTH_INT", "PERMTH_EXM"]

mort_1999 = pd.read_fwf(
    RAW_DATA_DIR / "NHANES_1999_2000_MORT_2019_PUBLIC.dat",
    colspecs=colspecs, names=colnames, na_values=['.', ' '],
    dtype={'SEQN': str, 'MORTSTAT': 'Int64', 'ELIGSTAT': 'Int64'}, skipinitialspace=True
)

mort_2001 = pd.read_fwf(
    RAW_DATA_DIR / "NHANES_2001_2002_MORT_2019_PUBLIC.dat",
    colspecs=colspecs, names=colnames, na_values=['.', ' '],
    dtype={'SEQN': str, 'MORTSTAT': 'Int64', 'ELIGSTAT': 'Int64'}, skipinitialspace=True
)

mort_df = pd.concat([mort_1999, mort_2001], ignore_index=True)
mort_cols = ["SEQN", "ELIGSTAT", "MORTSTAT", "PERMTH_EXM", "PERMTH_INT"]
mort_df = mort_df[mort_cols].copy()

print(f"  ✓ Loaded {mort_df.shape[0]:,} rows (1999: {mort_1999.shape[0]:,}, 2001: {mort_2001.shape[0]:,})")

# Ensure SEQN is consistent type for merging
# DNAm SEQN is float; convert to zero-padded 5-digit string
dnam_df["SEQN"] = pd.to_numeric(dnam_df["SEQN"], errors='coerce').astype('Int64').astype(str).str.zfill(5)
# Demo SEQN is float; convert to zero-padded 5-digit string
demo_df["SEQN"] = pd.to_numeric(demo_df["SEQN"], errors='coerce').astype('Int64').astype(str).str.zfill(5)
# Mort SEQN is already zero-padded string, keep as-is
mort_df["SEQN"] = mort_df["SEQN"].astype(str).str.strip()

print("\n" + "="*80)
print("STEP 2: MERGE DATASETS")
print("="*80)

print(f"\nStarting: DNAm n={dnam_df.shape[0]:,}")

# Inner join DNAm + DEMO
print(f"\n[2a] Inner join DNAm to DEMO on SEQN")
cohort = dnam_df.merge(demo_df, on="SEQN", how="inner")
print(f"  After join: n={cohort.shape[0]:,}")
print(f"  Dropped: {dnam_df.shape[0] - cohort.shape[0]:,} (DNAm samples not in DEMO)")

# Left join mortality
print(f"\n[2b] Left join mortality on SEQN")
cohort_before_mort = cohort.shape[0]
cohort = cohort.merge(mort_df, on="SEQN", how="left")
print(f"  After join: n={cohort.shape[0]:,}")
print(f"  All rows retained (left join)")

print("\n" + "="*80)
print("STEP 3: BUILD ANALYTIC COHORT (EXPLICIT FILTERING)")
print("="*80)

print(f"\nStarting: n={cohort.shape[0]:,}")

# Filter 1: ELIGSTAT == 1
n_before = cohort.shape[0]
cohort = cohort[cohort["ELIGSTAT"] == 1].copy()
print(f"\n[3a] Keep ELIGSTAT == 1 (eligible)")
print(f"  Dropped: {n_before - cohort.shape[0]:,}")
print(f"  Remaining: {cohort.shape[0]:,}")

# Filter 2: Non-null GrimAgeMort (defines DNAm subsample)
n_before = cohort.shape[0]
cohort = cohort[cohort["GrimAgeMort"].notna()].copy()
print(f"\n[3b] Keep GrimAgeMort non-null (DNAm subsample)")
print(f"  Dropped: {n_before - cohort.shape[0]:,}")
print(f"  Remaining: {cohort.shape[0]:,} ← Expected ~2,500 (DNAm subsample size)")

# Filter 3: Non-null and non-zero PERMTH_EXM
n_before = cohort.shape[0]
cohort = cohort[(cohort["PERMTH_EXM"].notna()) & (cohort["PERMTH_EXM"] > 0)].copy()
print(f"\n[3c] Keep PERMTH_EXM > 0 (have follow-up from exam)")
print(f"  Dropped: {n_before - cohort.shape[0]:,}")
print(f"  Remaining: {cohort.shape[0]:,}")

# Check RIDAGEYR range
print(f"\n[3d] RIDAGEYR (chronological age) range check:")
age_min, age_max = cohort["RIDAGEYR"].min(), cohort["RIDAGEYR"].max()
print(f"  Min: {age_min:.0f}, Max: {age_max:.0f}")
print(f"  Expected: ~50-85 (top-coded at 85)")
if age_min >= 40 and age_max <= 90:
    print(f"  ✓ Range reasonable for older adult cohort")

print(f"\n{'='*80}")
print(f"FINAL ANALYTIC COHORT: n={cohort.shape[0]:,}")
print(f"{'='*80}")

print("\n" + "="*80)
print("STEP 4: CONSTRUCT SURVIVAL VARIABLES")
print("="*80)

# Time in years
cohort["time_years"] = cohort["PERMTH_EXM"] / 12.0

# Event (1=dead, 0=censored)
cohort["event"] = cohort["MORTSTAT"].fillna(0).astype(int)

# Summary
n_total = cohort.shape[0]
n_events = cohort["event"].sum()
pct_events = 100 * n_events / n_total

time_median = cohort["time_years"].median()
time_max = cohort["time_years"].max()

print(f"\nSurvival variables constructed:")
print(f"  Total n: {n_total:,}")
print(f"  Events (deaths): {n_events:,} ({pct_events:.1f}%)")
print(f"  Censored: {n_total - n_events:,} ({100-pct_events:.1f}%)")
print(f"  Median follow-up: {time_median:.1f} years")
print(f"  Max follow-up: {time_max:.1f} years")

print("\n" + "="*80)
print("STEP 5: STANDARDIZE PREDICTORS")
print("="*80)

# Clock variables to standardize
clock_vars = ["HorvathAge", "HannumAge", "SkinBloodAge", "PhenoAge",
              "GrimAgeMort", "GrimAge2Mort", "DunedinPoAm", "HorvathTelo"]

# Standardize age and clocks
scaler = StandardScaler()
vars_to_scale = ["RIDAGEYR"] + clock_vars

scaled_data = scaler.fit_transform(cohort[vars_to_scale])

# Create z-scored columns
for i, var in enumerate(vars_to_scale):
    cohort[f"z_{var}"] = scaled_data[:, i]

print(f"Created standardized (z-scored) versions of:")
print(f"  - Chronological age (RIDAGEYR)")
print(f"  - {len(clock_vars)} epigenetic clocks")

print(f"\nExample (GrimAgeMort):")
print(f"  Raw: mean={cohort['GrimAgeMort'].mean():.2f}, sd={cohort['GrimAgeMort'].std():.2f}")
print(f"  Standardized: mean={cohort['z_GrimAgeMort'].mean():.6f}, sd={cohort['z_GrimAgeMort'].std():.6f}")

print("\n" + "="*80)
print("STEP 6: COHORT DESCRIPTIVE TABLE")
print("="*80)

# Chronological age
age_mean = cohort["RIDAGEYR"].mean()
age_sd = cohort["RIDAGEYR"].std()
female_pct = 100 * (cohort["RIAGENDR"] == 2).sum() / cohort.shape[0]

print(f"\nDemographics:")
print(f"  Chronological age: {age_mean:.1f} ± {age_sd:.1f} years")
print(f"  Female: {female_pct:.1f}%")

# Epigenetic clocks
print(f"\nEpigenetic clocks (raw scale):")
print(f"  {'Clock':<20} {'Mean':>8} {'SD':>8}")
print(f"  {'-'*37}")
for var in clock_vars:
    mean_val = cohort[var].mean()
    sd_val = cohort[var].std()
    print(f"  {var:<20} {mean_val:8.1f} {sd_val:8.1f}")

# Correlations with chronological age
print(f"\nPearson correlation with chronological age (RIDAGEYR):")
print(f"  {'Clock':<20} {'Correlation':>12}")
print(f"  {'-'*33}")
for var in clock_vars:
    corr = cohort["RIDAGEYR"].corr(cohort[var])
    print(f"  {var:<20} {corr:12.3f}")

# Events vs non-events
print(f"\nComparison by mortality status (among eligible with follow-up):")
print(f"  Variable            Alive (n={n_total - n_events:,})        Dead (n={n_events:,})")
print(f"  {'-'*70}")

event_groups = cohort.groupby("event")[["RIDAGEYR", "GrimAgeMort"]].mean()
for idx, row in event_groups.iterrows():
    event_label = "Dead" if idx == 1 else "Alive"
    print(f"  RIDAGEYR ({event_label}): {row['RIDAGEYR']:>25.1f}")
print()
for idx, row in event_groups.iterrows():
    event_label = "Dead" if idx == 1 else "Alive"
    print(f"  GrimAgeMort ({event_label}): {row['GrimAgeMort']:>20.1f}")

print(f"\nExpected: Dead group should be older on both measures")
if event_groups.loc[1, "RIDAGEYR"] > event_groups.loc[0, "RIDAGEYR"]:
    print(f"  ✓ Deceased participants are older (chronological age)")
if event_groups.loc[1, "GrimAgeMort"] > event_groups.loc[0, "GrimAgeMort"]:
    print(f"  ✓ Deceased participants have higher GrimAgeMort")

print("\n" + "="*80)
print("STEP 7: TRAIN/TEST SPLIT")
print("="*80)

# Stratified split (70/30) on event
X = cohort.drop(["event", "time_years"], axis=1)
y = cohort["event"]

train_idx, test_idx = train_test_split(
    np.arange(len(cohort)),
    test_size=0.3,
    stratify=y,
    random_state=42
)

cohort_train = cohort.iloc[train_idx].copy()
cohort_test = cohort.iloc[test_idx].copy()

print(f"\nTrain/test split (stratified on event, random_state=42):")
print(f"\n  Train set:")
print(f"    n: {len(cohort_train):,}")
print(f"    Events: {cohort_train['event'].sum():,} ({100*cohort_train['event'].sum()/len(cohort_train):.1f}%)")
print(f"    Censored: {(1-cohort_train['event']).sum():,} ({100*(1-cohort_train['event']).sum()/len(cohort_train):.1f}%)")

print(f"\n  Test set:")
print(f"    n: {len(cohort_test):,}")
print(f"    Events: {cohort_test['event'].sum():,} ({100*cohort_test['event'].sum()/len(cohort_test):.1f}%)")
print(f"    Censored: {(1-cohort_test['event']).sum():,} ({100*(1-cohort_test['event']).sum()/len(cohort_test):.1f}%)")

print(f"\n  Combined: n={len(cohort_train) + len(cohort_test):,} (=total)")

# Save cohort
print("\n" + "="*80)
print("STEP 8: SAVE ANALYTIC COHORT")
print("="*80)

cohort_path = PROC_DATA_DIR / "analytic_cohort.csv"
cohort.to_csv(cohort_path, index=False)
print(f"\n✓ Saved full cohort: {cohort_path}")
print(f"  {cohort.shape[0]:,} rows × {cohort.shape[1]} columns")

# Save split indices
split_data = {
    "train_indices": train_idx.tolist(),
    "test_indices": test_idx.tolist(),
}

import json
split_path = PROC_DATA_DIR / "train_test_split.json"
with open(split_path, 'w') as f:
    json.dump(split_data, f)
print(f"✓ Saved split indices: {split_path}")

print("\n" + "="*80)
print("ANALYTIC COHORT CONSTRUCTION COMPLETE")
print("="*80)
print(f"\nFiles saved:")
print(f"  - {cohort_path}")
print(f"  - {split_path}")
print(f"\nReady for Cox proportional hazards modeling.")
print("="*80 + "\n")

