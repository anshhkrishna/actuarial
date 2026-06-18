#!/usr/bin/env python3
"""
End-to-end pipeline runner for the NHANES epigenetic clock mortality
validation. Runs all four stages in order:

    1. download  — document data sources + validate raw files present
    2. verify    — fetch/inspect DNAm files and confirm variable names
    3. build     — merge sources, build analytic cohort + train/test split
    4. analysis  — fit Cox models, evaluate test C-index, save results

Raw NHANES files cannot be redistributed and (for the mortality files)
require agreeing to NCHS data-use terms, so they are NOT in the repo.
Place the files listed below in data/raw/ before running. See README.md
for the exact download URLs.

Usage:
    python run_all.py
"""

import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SRC = ROOT / "src"

# Raw files required by build_cohort.py
REQUIRED_RAW = [
    "dnmepi.sas7bdat",
    "DEMO_1999.xpt",
    "DEMO_2001.xpt",
    "NHANES_1999_2000_MORT_2019_PUBLIC.dat",
    "NHANES_2001_2002_MORT_2019_PUBLIC.dat",
]


def preflight():
    raw = ROOT / "data" / "raw"
    missing = [f for f in REQUIRED_RAW if not (raw / f).exists()]
    if missing:
        print("=" * 78)
        print("MISSING RAW DATA — cannot build the cohort.")
        print("=" * 78)
        print("Place these files in data/raw/ (download URLs in README.md):")
        for f in missing:
            print(f"  - {f}")
        print("\nNHANES raw files are not redistributed in this repo.")
        sys.exit(1)


def run(label, *cmd, fatal=True):
    print("\n" + "#" * 78)
    print(f"# STAGE: {label}")
    print("#" * 78, flush=True)
    result = subprocess.run([sys.executable, *cmd], cwd=ROOT)
    if result.returncode != 0:
        msg = f"Stage '{label}' exited with code {result.returncode}"
        if fatal:
            print(msg, file=sys.stderr)
            sys.exit(result.returncode)
        print(f"[warning] {msg} (non-fatal, continuing)", file=sys.stderr)


def main():
    preflight()
    # 1. document sources + validate raw files (non-fatal: informational)
    run("download", str(SRC / "download.py"), fatal=False)
    # 2. verify DNAm variable names (needs network; non-fatal if offline)
    run("verify", str(SRC / "verify_variables.py"), fatal=False)
    # 3. build the analytic cohort + train/test split
    run("build_cohort", str(SRC / "build_cohort.py"))
    # 4. survival analysis + figures + results table
    run("analysis", str(SRC / "analysis.py"))

    print("\n" + "=" * 78)
    print("PIPELINE COMPLETE")
    print("=" * 78)
    print("Outputs:")
    print("  results/cindex_comparison.csv")
    print("  results/ph_assumption_check.txt")
    print("  figures/incremental_cindex.png")
    print("  figures/km_by_grimage_tertile.png")


if __name__ == "__main__":
    main()
