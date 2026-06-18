"""
Survival analysis: do DNA-methylation aging clocks add mortality
discrimination beyond chronological age and sex?

Fits Cox proportional hazards models with lifelines on a fixed TRAIN
split and evaluates Harrell's concordance (C-index) on the held-out
TEST split.

Models per clock:
  - base:        z_RIDAGEYR + sex
  - age+sex+clock: z_RIDAGEYR + sex + z_<clock>
  - clock alone:   z_<clock>

Outputs:
  - results/cindex_comparison.csv
  - figures/incremental_cindex.png
  - figures/km_by_grimage_tertile.png
"""

import json
from pathlib import Path

import numpy as np
import pandas as pd
import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from lifelines import CoxPHFitter, KaplanMeierFitter
from lifelines.plotting import add_at_risk_counts

ROOT = Path(__file__).resolve().parent.parent
PROC_DIR = ROOT / "data" / "processed"
RESULTS_DIR = ROOT / "results"
FIG_DIR = ROOT / "figures"
RESULTS_DIR.mkdir(exist_ok=True)
FIG_DIR.mkdir(exist_ok=True)

COHORT_PATH = PROC_DIR / "analytic_cohort.csv"
SPLIT_PATH = PROC_DIR / "train_test_split.json"

CLOCKS = [
    "HorvathAge",
    "HannumAge",
    "SkinBloodAge",
    "PhenoAge",
    "GrimAgeMort",
    "GrimAge2Mort",
    "DunedinPoAm",
    "HorvathTelo",
]

DURATION = "time_years"
EVENT = "event"


def load_data():
    """Load cohort and the fixed train/test positional split."""
    cohort = pd.read_csv(COHORT_PATH)
    # sex as a clean binary covariate: 1 = female (RIAGENDR==2), 0 = male
    cohort["sex_female"] = (cohort["RIAGENDR"] == 2).astype(int)

    with open(SPLIT_PATH) as f:
        split = json.load(f)
    train = cohort.iloc[split["train_indices"]].copy()
    test = cohort.iloc[split["test_indices"]].copy()
    return cohort, train, test


def fit_cox(train, covariates):
    """Fit a Cox model on the given covariates; return the fitted model."""
    cols = covariates + [DURATION, EVENT]
    cph = CoxPHFitter()
    cph.fit(train[cols], duration_col=DURATION, event_col=EVENT)
    return cph


def test_cindex(cph, test, covariates):
    """Harrell's C-index of a fitted model on the held-out test set."""
    cols = covariates + [DURATION, EVENT]
    return cph.score(test[cols], scoring_method="concordance_index")


def run_models(train, test):
    """Fit base, age+sex+clock, and clock-alone models; collect results."""
    base_cov = ["z_RIDAGEYR", "sex_female"]
    base = fit_cox(train, base_cov)
    c_base = test_cindex(base, test, base_cov)
    print(f"\nBase model (age + sex): test C-index = {c_base:.4f}")

    rows = []
    for clock in CLOCKS:
        zc = f"z_{clock}"

        # age + sex + clock
        cov = base_cov + [zc]
        m = fit_cox(train, cov)
        c_with = test_cindex(m, test, cov)
        summ = m.summary.loc[zc]
        hr = float(np.exp(summ["coef"]))
        ci_low = float(np.exp(summ["coef lower 95%"]))
        ci_high = float(np.exp(summ["coef upper 95%"]))
        p = float(summ["p"])

        # clock alone
        m_alone = fit_cox(train, [zc])
        c_alone = test_cindex(m_alone, test, [zc])

        rows.append(
            {
                "clock": clock,
                "C_base": c_base,
                "C_with_clock": c_with,
                "delta_C": c_with - c_base,
                "HR": hr,
                "CI_low": ci_low,
                "CI_high": ci_high,
                "p": p,
                "C_clock_alone": c_alone,
            }
        )

    results = pd.DataFrame(rows).sort_values("delta_C", ascending=False).reset_index(drop=True)
    return base, base_cov, results


def check_ph(base, base_cov, train, best_clock):
    """Run PH assumption checks for the base model and the best clock model.

    Returns a text report; violations are reported, not hidden.
    """
    import io
    import contextlib

    report = io.StringIO()

    def section(title, cph, cols):
        report.write("=" * 78 + "\n")
        report.write(title + "\n")
        report.write("=" * 78 + "\n")
        buf = io.StringIO()
        with contextlib.redirect_stdout(buf):
            try:
                cph.check_assumptions(train[cols + [DURATION, EVENT]], p_value_threshold=0.05)
            except Exception as e:  # pragma: no cover - defensive
                buf.write(f"[check_assumptions raised: {e}]\n")
        out = buf.getvalue().strip()
        report.write((out if out else "No PH violations flagged at p < 0.05.") + "\n\n")

    section("PROPORTIONAL HAZARDS CHECK — BASE MODEL (age + sex)", base, base_cov)

    zc = f"z_{best_clock}"
    best_cov = base_cov + [zc]
    best_model = fit_cox(train, best_cov)
    section(
        f"PROPORTIONAL HAZARDS CHECK — BEST CLOCK MODEL (age + sex + {zc})",
        best_model,
        best_cov,
    )

    return report.getvalue()


def plot_incremental_cindex(results):
    """Bar chart of incremental test C-index (delta_C) by clock."""
    df = results.sort_values("delta_C")
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = ["#2a7fb8" if d >= 0 else "#c0504d" for d in df["delta_C"]]
    ax.barh(df["clock"], df["delta_C"], color=colors)
    ax.axvline(0, color="black", linewidth=1)
    ax.set_xlabel("Δ test C-index vs. age + sex base model")
    ax.set_title("Incremental discrimination added by each clock\n(age + sex + clock vs. age + sex)")
    for y, (d,) in enumerate(zip(df["delta_C"])):
        ax.text(d + (0.0005 if d >= 0 else -0.0005), y, f"{d:+.4f}",
                va="center", ha="left" if d >= 0 else "right", fontsize=8)
    ax.margins(x=0.15)
    fig.tight_layout()
    out = FIG_DIR / "incremental_cindex.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved {out}")


def plot_km_grimage_tertile(cohort):
    """Kaplan-Meier survival curves by GrimAgeMort tertile, with at-risk counts."""
    df = cohort.copy()
    df["tertile"] = pd.qcut(df["GrimAgeMort"], 3, labels=["Low", "Middle", "High"])

    fig, ax = plt.subplots(figsize=(8, 6))
    fitters = []
    colors = {"Low": "#2a7fb8", "Middle": "#f0a202", "High": "#c0504d"}
    for label in ["Low", "Middle", "High"]:
        m = df[df["tertile"] == label]
        kmf = KaplanMeierFitter()
        kmf.fit(m[DURATION], m[EVENT], label=f"{label} GrimAgeMort (n={len(m)})")
        kmf.plot_survival_function(ax=ax, ci_show=True, color=colors[label])
        fitters.append(kmf)

    ax.set_xlabel("Years of follow-up")
    ax.set_ylabel("Survival probability")
    ax.set_title("Survival by GrimAgeMort tertile — NHANES 1999–2002, age 50+")
    ax.set_ylim(0, 1)
    add_at_risk_counts(*fitters, ax=ax)
    fig.tight_layout()
    out = FIG_DIR / "km_by_grimage_tertile.png"
    fig.savefig(out, dpi=150)
    plt.close(fig)
    print(f"Saved {out}")


def main():
    cohort, train, test = load_data()
    print(f"Loaded cohort n={len(cohort)} (train={len(train)}, test={len(test)})")
    print(f"Events: train={int(train[EVENT].sum())}, test={int(test[EVENT].sum())}")

    base, base_cov, results = run_models(train, test)

    # Save results table
    out_csv = RESULTS_DIR / "cindex_comparison.csv"
    results.to_csv(out_csv, index=False)

    # Print final table
    print("\n" + "=" * 78)
    print("RESULTS — incremental test C-index over age + sex (sorted by delta_C)")
    print("=" * 78)
    disp = results.copy()
    with pd.option_context("display.width", 200, "display.max_columns", None):
        print(
            disp.to_string(
                index=False,
                formatters={
                    "C_base": "{:.4f}".format,
                    "C_with_clock": "{:.4f}".format,
                    "delta_C": "{:+.4f}".format,
                    "HR": "{:.3f}".format,
                    "CI_low": "{:.3f}".format,
                    "CI_high": "{:.3f}".format,
                    "p": "{:.2e}".format,
                    "C_clock_alone": "{:.4f}".format,
                },
            )
        )
    print(f"\nSaved {out_csv}")

    # PH checks for base + best clock
    best_clock = results.iloc[0]["clock"]
    ph_report = check_ph(base, base_cov, train, best_clock)
    print("\n" + ph_report)
    (RESULTS_DIR / "ph_assumption_check.txt").write_text(ph_report)
    print(f"Saved {RESULTS_DIR / 'ph_assumption_check.txt'}")

    # Figures
    plot_incremental_cindex(results)
    plot_km_grimage_tertile(cohort)


if __name__ == "__main__":
    main()
