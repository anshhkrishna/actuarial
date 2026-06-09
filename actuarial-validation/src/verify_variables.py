"""
Fetch NHANES DNAm .XPT files and empirically verify column names and labels.

This script:
1. Downloads the actual DNAm .XPT files from NHANES
2. Loads them with pyreadstat to access column metadata
3. Prints all columns and labels
4. Maps target epigenetic clocks to actual column names
5. Verifies WTDN4YR and SEQN are present
"""

import sys
from pathlib import Path
import requests
import pyreadstat
import pandas as pd

RAW_DATA_DIR = Path(__file__).parent.parent / "data" / "raw"
RAW_DATA_DIR.mkdir(parents=True, exist_ok=True)

# Target clocks we're looking for
TARGET_CLOCKS = [
    "Horvath", "HorvathAge", "Hannum", "SkinBlood",
    "PhenoAge", "GrimAge", "GrimAge2", "DunedinPoAm", "DNAmTL"
]

# NHANES cycles
CYCLES = ["1999-2000", "2001-2002"]
DNAM_FILENAME = "L10EPG.XPT"
NHANES_BASE_URL = "https://wwwn.cdc.gov/Nchs/Nhanes"


def download_xpt_file(cycle: str) -> Path:
    """Download DNAm .XPT file for a given cycle."""
    url = f"{NHANES_BASE_URL}/{cycle}/{DNAM_FILENAME}"
    filepath = RAW_DATA_DIR / f"dnam_{cycle}.XPT"

    print(f"\nDownloading: {url}")
    try:
        response = requests.get(url, timeout=30, stream=True)
        response.raise_for_status()

        with open(filepath, 'wb') as f:
            for chunk in response.iter_content(chunk_size=8192):
                if chunk:
                    f.write(chunk)

        size_mb = filepath.stat().st_size / (1024**2)
        print(f"✓ Downloaded {filepath.name} ({size_mb:.2f} MB)")
        return filepath

    except requests.exceptions.HTTPError as e:
        print(f"✗ HTTP {e.response.status_code}: {url}")
        return None
    except Exception as e:
        print(f"✗ Error: {e}")
        return None


def inspect_xpt_file(filepath: Path) -> None:
    """Load XPT file and print column names + labels."""
    print(f"\n{'='*80}")
    print(f"INSPECTING: {filepath.name}")
    print(f"{'='*80}")

    try:
        df, meta = pyreadstat.read_xport(str(filepath))

        print(f"\nFile shape: {df.shape[0]:,} rows × {df.shape[1]} columns")

        # Print all column names and labels
        print(f"\n{'Column Name':<20} {'Column Label':<60}")
        print("-" * 80)

        for col in df.columns:
            label = meta.column_labels.get(col, "[no label]") if meta.column_labels else "[no label]"
            print(f"{col:<20} {label:<60}")

        # Check for required variables
        print(f"\n{'='*80}")
        print("VERIFICATION OF REQUIRED VARIABLES:")
        print(f"{'='*80}")

        required_vars = ["SEQN", "WTDN4YR"]
        for var in required_vars:
            status = "✓" if var in df.columns else "✗"
            print(f"{status} {var:<20} {'FOUND' if var in df.columns else 'MISSING'}")

        # Look for target clocks (case-insensitive search)
        print(f"\n{'='*80}")
        print("TARGET EPIGENETIC CLOCKS:")
        print(f"{'='*80}")

        found_clocks = {}
        df_columns_upper = {col.upper(): col for col in df.columns}

        for clock in TARGET_CLOCKS:
            clock_upper = clock.upper()
            if clock_upper in df_columns_upper:
                actual_col = df_columns_upper[clock_upper]
                label = meta.column_labels.get(actual_col, "[no label]") if meta.column_labels else "[no label]"
                found_clocks[clock] = (actual_col, label)
                print(f"✓ {clock:<15} → {actual_col:<20} | {label}")
            else:
                # Try substring search
                matches = [col for col in df.columns if clock.upper() in col.upper()]
                if matches:
                    for match in matches:
                        label = meta.column_labels.get(match, "[no label]") if meta.column_labels else "[no label]"
                        found_clocks[f"{clock}*"] = (match, label)
                        print(f"◐ {clock:<15} → {match:<20} | {label} (substring match)")
                else:
                    print(f"✗ {clock:<15} NOT FOUND")

        # Print summary mapping table
        if found_clocks:
            print(f"\n{'='*80}")
            print("VERIFIED MAPPING (Target Clock → Actual Column → Label):")
            print(f"{'='*80}")
            print(f"{'Target Clock':<20} {'Actual Column':<20} {'Column Label':<40}")
            print("-" * 80)
            for target, (actual, label) in sorted(found_clocks.items()):
                print(f"{target:<20} {actual:<20} {label[:40]:<40}")

        return found_clocks

    except Exception as e:
        print(f"✗ Error reading file: {e}")
        return {}


def main():
    print("\n" + "="*80)
    print("NHANES DNAm VARIABLE VERIFICATION")
    print("="*80)

    all_findings = {}

    for cycle in CYCLES:
        print(f"\n{'='*80}")
        print(f"CYCLE: {cycle}")
        print(f"{'='*80}")

        # Download file
        filepath = download_xpt_file(cycle)

        if filepath and filepath.exists():
            # Inspect file
            findings = inspect_xpt_file(filepath)
            all_findings[cycle] = findings
        else:
            print(f"\n✗ Could not download {cycle}. Exact URL attempted:")
            url = f"{NHANES_BASE_URL}/{cycle}/{DNAM_FILENAME}"
            print(f"  {url}")

    # Final summary
    print(f"\n\n{'='*80}")
    print("FINAL SUMMARY")
    print(f"{'='*80}")

    for cycle, clocks in all_findings.items():
        if clocks:
            print(f"\n✓ {cycle}: Found {len(clocks)} target clocks")
        else:
            print(f"\n✗ {cycle}: No target clocks found")

    print(f"\n{'='*80}\n")


if __name__ == "__main__":
    main()
