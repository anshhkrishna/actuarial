# Results

- Base Cox model (chronological age + sex) reaches a **test C-index of 0.746** for ~17–20 year all-cause mortality in NHANES 1999–2002 (n=2,532, age 50+).
- **GrimAge2Mort** adds the most discrimination: test C-index rises to **0.769 (Δ +0.0235)**, HR per SD = **2.06 (95% CI 1.84–2.30, p≈2.5e-36)**; GrimAgeMort is close behind (0.766, Δ +0.0202).
- The pace-of-aging clock **DunedinPoAm** adds Δ +0.0116; the four age-trained clocks (Hannum, PhenoAge, Horvath, SkinBlood) add ≤0.0035, and SkinBlood is not significant beyond age+sex (p=0.25).
- **Takeaway:** mortality-trained clocks carry real signal beyond age and sex, but the absolute gain in discrimination is modest (≤0.024 C-index).
- Proportional-hazards check: the base model is fine; in the GrimAge2Mort model both age and the clock flag mild PH violations (p≈0.008 and 0.0004) — reported in `results/ph_assumption_check.txt`, not hidden.
