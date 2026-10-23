# Contributing

Please separate numerical methods, dataset adapters and UI changes. Include
units, reference frames, source versions and scientific limitations with new
data. Do not label archived data as held-out validation without verifying model
training membership.

Run `npm ci`, then `npm test` for the portable frontend contracts. For browser
checks, install Playwright Chromium with `npx playwright install chromium`,
serve `simulator/public` on loopback port 8765, set
`PLAYWRIGHT_CHANNEL=chromium`, and run `npm run test:release`.

Large raw datasets, Python environments, downloaded third-party repositories,
credentials and data lacking redistribution permission must stay out of Git.
Use fixed source URLs and checksums instead. The public dataset catalog is
`simulator/public/datasets/catalog.public.json`; the ignored local `catalog.json`
may additionally refer to privately held research snapshots.

Numerical changes need tests against independent calculations or known
synthetic cases. A lower mismatch alone does not establish improved physics.

For the population studies, follow `research/REPRODUCING_POPULATION.md` and install
`research/population-requirements.txt` in an isolated environment. Stage-one
functional checks and the synthetic Monte Carlo checks use tracked inputs;
real-data stage-two runs need separately retrieved, checksum-verified tutorial
arrays. Do not silently overwrite the frozen `results/model.json` or replace a
recorded input checksum when reproducing a paper result.

Preserve the distinction between an empirical resampling diagnostic and actual
astrophysical uncertainty. Report shared-sample covariance, input versions,
selection rules and prior support when comparing population models. Additional
bootstrap repetitions do not add independent observations or recover unseen
posterior tails. Run short checks into a separate output directory; use the
published seeds and repetition counts when comparing with the full results.

Documentary snapshots should be fetched from their public source, not copied
wholesale into a contribution. If an upstream excerpt is necessary, retain its
license and attribution in `NOTICE.md`. Keep original HDF/pickle files and direct
array conversions out of Git, even when they are publicly downloadable.

## Browser release checks on hosted runners

The GitHub job sets `WORKBENCH_DISABLE_WEBGL=1` to test the supported coordinate
fallback on machines without a hardware GPU. The test asserts that optics are
unavailable, coordinates remain visible, and all public scenes and analysis
views still work. This does not validate optical image fidelity. Omit that
environment variable for the default graphics path; both paths were checked
locally for the population-paper release. The default graphics path timed out
during scene selection on the hosted Linux runner; the local default path passed.
