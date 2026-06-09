"""
Download and manage NHANES DNA methylation epigenetic biomarker data.

NHANES data requires registration and download from the CDC website.
This script provides:
1. Documentation of exact variable names required
2. Instructions for manual data download
3. Validation framework for downloaded files
"""

import os
import sys
import json
from pathlib import Path
from typing import Dict, List, Optional

# Project paths
RAW_DATA_DIR = Path(__file__).parent.parent / "data" / "raw"
RAW_DATA_DIR.mkdir(parents=True, exist_ok=True)

# NHANES DNAm variable specification
# Source: NHANES DNA Methylation Epigenetic Biomarker data dictionary
DNAM_VARIABLES = {
    "seqn": "SEQN",                    # Respondent sequence number (participant ID)
    "grimage": "GRIMAGE",               # GrimAge predicted age
    "grimage2": "GRIMAGE2",             # GrimAge2 predicted age (revised version)
    "horvath": "HORVATHOA",             # Horvath original epigenetic age
    "hannum": "HANNUMAGE",              # Hannum epigenetic age acceleration
    "phenoage": "PHENOAGE",             # PhenoAge predicted age
    "dunedin": "DUNEDINPOAM",           # DunedinPoAm pace-of-aging (biological tempo)
    "dnam_weight": "WTDN4YR",           # DNA methylation-specific 4-year combined weight
}

# NHANES cycles to include in analysis
CYCLES = {
    "1999-2000": {
        "name": "NHANES 1999-2000",
        "dnam_file": "L10EPG.XPT",
        "demo_file": "DEMO.XPT",
    },
    "2001-2002": {
        "name": "NHANES 2001-2002",
        "dnam_file": "L10EPG.XPT",
        "demo_file": "DEMO.XPT",
    }
}

# Mortality file specification
MORTALITY_SPEC = {
    "columns": {
        "SEQN": (1, 5),           # Respondent sequence number
        "mortstat": (6, 6),       # Mortality status (1=deceased, 0=alive)
        "mortsdays": (7, 13),     # Days of follow-up (from exam to end of 2019)
        "ucod_leading": (14, 17), # Underlying cause of death (ICD-10 code)
    },
    "description": "NHANES public-use linked mortality file (NDI linkage through 2019)",
    "source": "https://www.cdc.gov/nchs/data-linkage/mortality-public.htm"
}


def print_dnam_variable_specification() -> None:
    """Print confirmed DNA methylation variable names and specifications."""
    print("\n" + "="*75)
    print("CONFIRMED NHANES DNA METHYLATION EPIGENETIC CLOCK VARIABLES")
    print("="*75)

    print("\nEpigenetic Aging Clock Variables:")
    print("-" * 75)
    print(f"{'Clock':<20} {'NHANES Variable':<20} {'Description':<35}")
    print("-" * 75)

    clocks = {
        "GrimAge": ("GRIMAGE", "Phenotypic age surrogate (2018)"),
        "GrimAge2": ("GRIMAGE2", "GrimAge2 revised/version 2 (2024)"),
        "Horvath": ("HORVATHOA", "Original Horvath epigenetic age"),
        "Hannum": ("HANNUMAGE", "Hannum age acceleration"),
        "PhenoAge": ("PHENOAGE", "Phenotypic age estimate"),
        "DunedinPoAm": ("DUNEDINPOAM", "Pace-of-aging methylation"),
    }

    for clock, (var, desc) in clocks.items():
        print(f"{clock:<20} {var:<20} {desc:<35}")

    print("\nKey Identifier and Survey Weight:")
    print("-" * 75)
    print(f"{'Participant ID':<20} {'SEQN':<20}")
    print(f"{'DNAm Sample Weight':<20} {'WTDN4YR':<20}")
    print(f"{'Note':<20} {'Use WTDN4YR for all survey analyses (4-yr combined weight)':<20}")

    print("\n" + "="*75)


def print_download_instructions() -> None:
    """Print detailed download instructions for NHANES data."""
    print("\n" + "="*75)
    print("DATA DOWNLOAD INSTRUCTIONS")
    print("="*75)

    print("\n1. NHANES DNA METHYLATION DATA")
    print("-" * 75)
    print("   Location: https://wwwn.cdc.gov/nchs/nhanes/dnam/")
    print("   Cycles:   1999-2000, 2001-2002")
    print("   Format:   SAS Transport (.XPT)")
    print("   Files:    L10EPG.XPT (DNAm biomarkers)")
    print("")
    print("   Steps:")
    print("   a) Visit https://wwwn.cdc.gov/nchs/nhanes/")
    print("   b) Select cycle → Laboratory Data → 'DNAm Epigenetic Biomarker'")
    print("   c) Download L10EPG.XPT for each cycle")
    print("   d) Save to: data/raw/dnam_[cycle].XPT")
    print("")

    print("2. NHANES DEMOGRAPHICS DATA")
    print("-" * 75)
    print("   Location: https://wwwn.cdc.gov/nchs/nhanes/search/datapage.aspx")
    print("   Cycles:   1999-2000, 2001-2002")
    print("   Format:   SAS Transport (.XPT)")
    print("   Files:    DEMO.XPT (demographics)")
    print("")
    print("   Steps:")
    print("   a) Visit https://wwwn.cdc.gov/nchs/nhanes/")
    print("   b) Select cycle → Demographics Data")
    print("   c) Download DEMO.XPT for each cycle")
    print("   d) Save to: data/raw/demo_[cycle].XPT")
    print("")

    print("3. LINKED MORTALITY FILES")
    print("-" * 75)
    print("   Location: https://www.cdc.gov/nchs/data-linkage/mortality-public.htm")
    print("   Cycles:   1999-2000, 2001-2002")
    print("   Format:   Fixed-width text file")
    print("   Linkage:  National Death Index (NDI) through 2019")
    print("")
    print("   Column Specification (fixed-width):")
    for col, (start, end) in MORTALITY_SPEC["columns"].items():
        print(f"     {col:<15} cols {start:2d}-{end:2d}")
    print("")
    print("   Steps:")
    print("   a) Visit https://www.cdc.gov/nchs/data-linkage/mortality-public.htm")
    print("   b) Download publicly available mortality file for US adults")
    print("   c) Filter or extract rows for 1999-2000 and 2001-2002 participants")
    print("   d) Save to: data/raw/mortality_public_[cycle].txt")
    print("")


def print_manifest(downloaded_files: Dict[str, List[str]]) -> None:
    """Print manifest of all data files present."""
    print("\n" + "="*75)
    print("DATA MANIFEST")
    print("="*75)

    print(f"\nRaw data directory: {RAW_DATA_DIR}")

    total_files = 0
    for file_type, files in downloaded_files.items():
        count = len(files)
        total_files += count
        status = "✓" if count > 0 else "✗"
        print(f"\n{status} {file_type}: {count} file(s)")
        for fname in sorted(files):
            fpath = RAW_DATA_DIR / fname
            if fpath.exists():
                size_mb = fpath.stat().st_size / (1024**2)
                print(f"    {fname:<40} {size_mb:8.2f} MB")
            else:
                print(f"    {fname:<40} (not found)")

    print("\n" + "="*75)
    print(f"Total files ready: {total_files} / 6")
    print("="*75)

    if total_files == 6:
        print("\n✓ All data files downloaded. Ready for processing.")
    else:
        print(f"\n⚠ {6 - total_files} file(s) still needed. Follow instructions above.")


def validate_xpt_files() -> None:
    """Attempt to validate downloaded XPT files."""
    try:
        import pyreadstat
    except ImportError:
        print("\n[Validation] pyreadstat not available - skipping XPT validation")
        return

    print("\n[Validating XPT files...]")

    for fname in sorted(RAW_DATA_DIR.glob("*.XPT")):
        print(f"  Checking {fname.name}...", end=" ")
        try:
            df, meta = pyreadstat.read_xport(str(fname))
            n_rows, n_cols = df.shape
            print(f"✓ ({n_rows:,} rows, {n_cols} cols)")

            # Check for expected variables
            if "SEQN" in df.columns:
                print(f"      Found SEQN (participant ID)")

            # Check for DNAm variables if this is a DNAm file
            if "dnam" in fname.name.lower():
                dnam_vars_found = [v for v in DNAM_VARIABLES.values() if v in df.columns]
                if dnam_vars_found:
                    print(f"      Found {len(dnam_vars_found)} epigenetic clock variables: {', '.join(dnam_vars_found)}")

        except Exception as e:
            print(f"✗ Error: {str(e)[:60]}")


def check_downloaded_files() -> Dict[str, List[str]]:
    """Check which data files are present in raw data directory."""
    downloaded = {
        "DNAm files": [],
        "Demographics files": [],
        "Mortality files": []
    }

    for fpath in RAW_DATA_DIR.glob("*"):
        fname = fpath.name
        if "dnam" in fname.lower():
            downloaded["DNAm files"].append(fname)
        elif "demo" in fname.lower():
            downloaded["Demographics files"].append(fname)
        elif "mortality" in fname.lower() or "mort" in fname.lower():
            downloaded["Mortality files"].append(fname)

    return downloaded


def main():
    """Main execution."""
    print("\n" + "="*75)
    print("NHANES EPIGENETIC AGING CLOCK - DATA ACQUISITION")
    print("="*75)

    # Step 1: Show variable specification
    print_dnam_variable_specification()

    # Step 2: Check for existing files
    print("\n[Checking for downloaded files...]")
    downloaded = check_downloaded_files()

    # Step 3: Validate existing files
    if any(downloaded.values()):
        validate_xpt_files()

    # Step 4: Print instructions
    print_download_instructions()

    # Step 5: Print manifest
    print_manifest(downloaded)

    print("\nNEXT STEPS:")
    print("  1. Download missing data files (see instructions above)")
    print("  2. Place files in: data/raw/")
    print("  3. Run: python src/download.py  (to re-validate)")
    print("="*75 + "\n")


if __name__ == "__main__":
    main()
