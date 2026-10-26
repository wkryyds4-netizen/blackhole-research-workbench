# Official O3 mixture: conditional all-69-source comparison

27 October 2026. A new local research stage after the goal resumed. Published v0.3 artifacts remain unchanged; this follow-up is not yet incorporated into that paper or the workbench.

The recommended four-population mixture in [fixed official release7890437](https://zenodo.org/records/7890437) replaces only the pure-BBH O3 selection term. All69 original-prior event integrals, both six-model families, the early O1/O2 term and time fraction remain fixed. The mixture contains1,141,491 accepted rows from86,318,296 generated trials. Every one of the284,672 pure-BBH rows has a unique exact match, with40 common physical/search fields identical.

**Stage finding: the change of official injection representation moves g posteriors much less than the tested magnitude-model contrast; the tiny native-r ordering remains search-policy dependent. The mixture does not close the target-support gap.** This compares related fixed banks, not independent sensitivity measurements or calibrated campaign errors.

## Posterior and model comparisons

pp denotes percentage points. TV is normalized posterior total variation on0–1. Different scientific functionals are never subtracted. The quantile-shift column covers only the stored5/50/95 quantiles, not a uniform bound on all quantiles.

| Metric | O3 policy | Max \|median shift\| (pp) | Max \|5/50/95 shift\| (pp) | Max posterior TV |
|---|---|---|---|---|
| g | matched_filter | 0.137027 | 0.196019 | 0.00543130 |
| g | with_cwb | 0.123030 | 0.176349 | 0.00488616 |
| r | matched_filter | 0.068321 | 0.137374 | 0.00977059 |
| r | with_cwb | 0.074574 | 0.146202 | 0.01042921 |

The g contrast subtracts the Beta(1.5,8.5) median from Beta(2,5) at fixed λ. The r contrast is σ=0.3 minus σ=1.2 at fixed magnitude family. Here g=P(χ_eff<0), while r=P(cos θ₁<0) retains mixed native references. Neither is an identified formation-channel fraction.

| Metric | O3 policy | Shape / magnitude | Pure BBH (pp) | Mixture (pp) |
|---|---|---|---|---|
| g | matched_filter | 0 | +12.097983 | +11.989168 |
| g | matched_filter | 4 | +12.892206 | +12.784112 |
| g | matched_filter | 16 | +14.342024 | +14.232472 |
| g | with_cwb | 0 | +11.907912 | +11.829218 |
| g | with_cwb | 4 | +12.682254 | +12.605962 |
| g | with_cwb | 16 | +14.092920 | +14.018385 |
| r | matched_filter | Beta(1.5,8.5) | -0.060242 | -0.021579 |
| r | matched_filter | Beta(2,5) | +6.353620 | +6.383751 |
| r | with_cwb | Beta(1.5,8.5) | +0.006044 | +0.052012 |
| r | with_cwb | Beta(2,5) | +6.293535 | +6.331654 |

The saved detected contribution in all newly covered target support is zero for both policies. Consequently, full-union and old-BBH-domain-restricted mixture results agree. This numerical equality does not establish zero physical detection probability. Shared rows make these related estimators; their difference has no assigned significance.

## Added support and remaining omissions

Within the fixed target m₁≤100,m₂≥3 solar masses, BNS/NSBH proposals add no support. IMBH adds only m₁≥90,m₂≥10 at1.9<z≤2.3. This mass region comprises0.306136% of the target mass distribution, adding just0.102273 percentage points of total population probability. The proposal union still misses33.305257% of the target probability mass. These are source-population masses, not detection fractions or selection-error bars.

There are178 saved rows in the newly covered target region and no saved detections under either search policy. Weighted proposal sampling and hopeless selection prevent using this as178 independent direct-target Bernoulli trials.

A subsequent source diagnostic also finds that, in260,145 overlap rows supported only by the BBH proposal, the released mixture joint PDF is not exactly the separate BBH PDF times the generated-trial fraction. Its maximum relative residual is6.4486%, with median ratio0.9999995. The numerical construction history remains unverified. Each official joint PDF is retained as provided; the diagnostic does not recover the actual generation configuration.

## Verification, prior work and reproduction

The independent validator reconstructs Cartesian target/joint-proposal ratios, early/O3 integrals and all69 multiplicative event factors without producer weighting/posterior functions. Maximum relative selection difference is8.89e-15; relative-to-peak density difference1.94e-12; quantile difference4.17e-15. Piecewise-polynomial/Gauss–Legendre support quadratures differ by4.14e-10. This is numerical implementation verification with unchanged event-integral inputs, not independent observations.

Mixture injection use follows the official release and[Essick2021](https://doi.org/10.3847/2515-5172/ac2ba7), not a new method. See[Essick2023](https://arxiv.org/html/2307.02765v1) for hopeless cuts and the distinction between real searches and semianalytic approximations. The incremental result is this specific all69/12-model comparison, row identity, posterior distances and quantified remaining support.

Reproduction needs the fixed raw data and prior all69 caches, plus population-requirements.txt and Matplotlib. Raw HDF, row maps and numerical checkpoints are excluded from Git; official MD5 and local SHA256 are recorded. Fetching this one file does not restore all preceding research inputs.

```text
python research/population-mixture/fetch_mixture.py
python research/population-mixture/inspect_mixture.py
python research/population-mixture/study_mixture.py
python research/population-mixture/validate_mixture.py
python research/population-mixture/plot_mixture.py
python research/population-mixture/summarize_mixture.py
```
[Results](results.json) · [Independent calculation](validation.json) · [Source audit](provenance.json) · [Before-comparison design](design.md) · [Vector PDF](mixture-comparison.pdf)

![Official O3 mixture comparison](mixture-comparison.svg)
