# Five-channel conditional inference: 24 complete comparisons

2026-10-27. All69 events and11,180,578 posterior rows are retained. Event/selection integrals, two posterior algorithms and the common-term decomposition are complete. Results are conditional on the restricted domain, fixed birth spin0.1/CE efficiency2 and declared spin completions; they are not cosmic formation-fraction measurements.

| Comparison | Largest median movement / pp | Channel |
|---|---:|---|
| Regularizer probability | 8.654 | CE |
| Hidden-spin completion | 3.415 | CE |
| O3 bank replacement | 0.614 | CE |
| Search policy | 0.427 | NSC |

Maxima occur at separate grid points and are neither additive errors nor credible intervals. The largest independent-algorithm quantile difference is 0.169pp; the maximum four-replicate importance span is 0.455pp. Replacing normalizers with independent draw halves shifts quantiles by at most 0.051pp.

For the official-mixture matched-filter case with the uniform base completion, the CE intrinsic median is83.85% without regularization and76.17% with a1% common term. The computed efficiency of the common term is about8.06%, versus0.139% for unregularized CE, so its median responsibility in the detected model mixture reaches12.14%. A small density regularizer therefore need not be a small selection correction. This does not measure a12.14% unknown physical channel.

The original `detected_quantiles` are **regularized channel labels**; shared U is allocated to labels according to beta. `detected-component-results.json` separately reports model terms, U, and the conditional detected model terms. Marginal medians need not sum to1.

## Limits and next step

The minimum event contribution concentration is5.24, and selection concentration218.21. These are not PE-chain/campaign ESS values and do not trigger deletion. Stable posterior computation does not calibrate bank tails. Historical injection reference/screening, missing PE support, waveform error and true channel-specific component spins remain unresolved. Both completions are explicit stress models, not verified formation physics.

The next substantive step is paired finite-bank propagation for these channel integrals, without reusing old g/r thinning intervals, then channel-specific spin completions and hyperparameter variations.

## Figures, papers and evidence

- [Posterior figure](channel-posteriors.svg)
- [Sensitivity and regularizer](channel-sensitivity.svg)
- [中文 LaTeX PDF](../../output/pdf/blackhole-channel-completion-latex-zh-2026-10-27.pdf) / [English LaTeX PDF](../../output/pdf/blackhole-channel-completion-latex-en-2026-10-27.pdf)
- [LaTeX projects and shared facts](../paper/channels-2026-10-27/)
- [Method and literature comparison](measure-method.md)
- [Complete synthesis](branching-synthesis.json)
- [Independent posterior validation](observed-mcmc-validation.json)


## Reproduction

Restore model products and the previously verified69-event/injection caches first; see the directory README. Run from the repository root. Cache protocols include source hashes; preserve this snapshot before recalibration or design changes, and never mix protocols. Do not execute root model.py. Original direct-product and short-chain failures remain preserved.

```text
python research/population-channels/validate_fast_measure.py
python research/population-channels/validate_channel_weights.py
python research/population-channels/run_channel_integrals.py
python research/population-channels/validate_channel_products.py
python research/population-channels/validate_branching_sampler.py
python research/population-channels/fit_branching.py
python research/population-channels/validate_independent_mcmc.py --refined
python research/population-channels/validate_observed_mcmc.py
python research/population-channels/derive_channel_components.py
python research/population-channels/summarize_branching.py
python research/paper/channels-2026-10-27/build_note.py --render
```
