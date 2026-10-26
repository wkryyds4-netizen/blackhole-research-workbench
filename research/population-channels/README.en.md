# Physical formation-channel interface: paused checkpoint

Saved on 2026-10-27 at the user's request. [中文](README.md).
This checkpoint records model integration and numerical diagnostics.
No formation-channel fractions have been fitted to observed events.

## Saved work

The original adapter loads the five CE, CHE, GC, NSC and SMT flow models from
[Colloms et al.'s public release](https://zenodo.org/records/14967688).
Their coordinates are source-frame chirp mass, mass ratio, effective spin and
redshift. Raw products, weights and downloaded source snapshots stay local;
Git contains provenance, checksums, original tools and aggregate diagnostics.
`channel_flow.py` uses Glasflow and strict numeric state-dictionary loading,
without executing downloaded AMAZE code.

`channel-preview.json` summarizes 491,520 emulator draws: three birth-spin
values, five channels and 32,768 rows per case. These are model previews,
not observational posteriors. `projected_measure.py` and
`measure-validation.json` implement and check a four-to-seven-dimensional
measure bridge under explicit nuisance-spin conditional distributions.
The bridge has not yet been applied to the 69-event channel likelihood.
`sign-threshold-results.json` is a post-preview exploratory diagnostic using
fixed draws; it was not preregistered.

## Retained limitations

- Not all original runtime gates passed. The NSC sampled-latent density error
  of 3.58e-8 exceeds the preset 1e-8 gate. Checks referenced to the stored
  physical coordinates pass, while `runtime-failures.json` preserves the
  original failure and the unchanged threshold.
- Reference-sample diagnostics support CE context `[chi_birth, ln(alpha_CE)]`.
  The historical training invocation remains unverified; a current example
  configuration is not evidence of that invocation.
- The released GC mapping allows q up to 1.001. About 0.24%–0.94% of GC preview
  rows have q>1. These rows are retained and counted; a physical truncation and
  renormalization have not been implemented.
- At zero birth spin, the GC strict negative-spin fraction is about 49.51%,
  versus 6.92% below -0.01, with 86.61% inside the corresponding neutral band.
  Continuous-flow broadening affects sign statistics. Neither the neutral
  band nor the negative fraction measures a true zero-spin atom, and neither
  converts the previous g posterior into channel fractions.
- `measure-history.json` retains the original unsplit-quadrature failure.
  Splitting at zero passes the original gate. These mathematical checks do
  not constitute an observational channel analysis.

## Reproduction and resumption

The local environment uses Python 3.14, CPU Torch 2.13.0, Glasflow 0.4.1,
NumPy, SciPy, h5py and Requests. Source tools are `fetch_sources.py`,
`fetch_products.py` and `inspect_products.py`; pinned URLs, revisions and hashes
are recorded in `sources/manifest.json` and `product-inspection.json`.
Raw products and sample caches are excluded from Git. This is not a verified
one-command clean-environment reproduction package.

After restoring and checking those products, diagnostics are rebuilt with
`check_runtime.py`, `sample_contexts.py`, `check_measure.py` and
`diagnose_sign_threshold.py`. Preserve this checkpoint before rerunning:
these tools overwrite their same-name outputs. The CE raw-context comparison
in the runtime check diagnoses mappings; the physical preview follows the
log-alpha convention in `preview-design.json`.

On resumption, resolve GC physical-domain normalization, a fast validated
spin-projection density, consistent event/injection nuisance conditionals and
common support before fitting a full channel likelihood. Density regularization
must enter both event and selection terms. Channel inference remains unfinished.

Primary comparisons are [Zevin et al.](https://arxiv.org/abs/2011.10057v2),
[Colloms et al.](https://arxiv.org/abs/2503.03819v1) and the pinned
[AMAZE source revision](https://github.com/michaelzevin/AMAZE/tree/4b3a0c78eec8ca4780f6e78fd50f0146ced33878).
The inspected source snapshot has no identified explicit license and is not
redistributed here.
