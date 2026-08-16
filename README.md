# actuarial

cox survival validation of 8 dna methylation aging clocks against 20 year nhanes mortality
follow up. n = 2,532 participants aged 50+, 1,361 deaths, median follow up 17.1 years.

**result:** mortality trained clocks (grimage2, grimage) carry real signal beyond age and
sex (hr 2.06 per sd, p = 2.5e-36) but the discrimination gain is modest: +0.024 c index
over an age + sex baseline. the 5 age trained clocks add essentially nothing. full numbers
in [`RESULTS.md`](RESULTS.md).

## stack

| layer | tools |
|---|---|
| survival models | lifelines, cox proportional hazards |
| data | pandas, numpy |
| stats | statsmodels, scipy |
| figures | matplotlib |
| source data | nhanes 1999 to 2002 dnam, nchs linked mortality files |

## pipeline

`download` → `verify` → `build_cohort` → `analysis`

merges dnam + demographics + mortality on `SEQN`, applies eligibility filters, z scores
predictors, stratified 70/30 split, fits cox on train, evaluates harrell c index on the
held out test split, checks the ph assumption.

## run

```bash
pip install -r requirements.txt
# place raw nhanes files in data/raw/ first
python run_all.py        # or: make all
```

writes `results/cindex_comparison.csv`, `results/ph_assumption_check.txt`, and figures.

## license

mit, and it covers the code only. the underlying nhanes data is produced by the u.s.
national center for health statistics and is subject to nchs terms of use — see
[`CITATION.md`](CITATION.md).
