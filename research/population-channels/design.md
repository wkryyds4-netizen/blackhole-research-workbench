# Physical formation-channel interface and source study

2026-10-27. Recorded before downloading physical-channel model samples or
evaluating any new population likelihood. This follows the completed invariant
g and54-model studies; it does not declare those descriptive models physical.

## Required evidence

1. Inspect primary authors' physical predictions and analysis code, starting
   with Zevin et al.2011.10057v2 / Zenodo4448170 / AMA​ZE. Pin actual record and
   code revisions. The linked repository may have moved; resolve its canonical
   name from primary metadata. Compare newer public emulation work if needed.
2. Identify what a channel file represents: intrinsic merger measure,
   cosmological event weighting, detector-frame/source-frame masses, spin
   variables/reference, finite sample weights and whether sensitivity is
   already included. Separate intrinsic and detectable mixture fractions.
3. Build an original data adapter, verify its normalization and source
   checksums, and produce physical joint/marginal diagnostics before inference.
   Do not load upstream pickle objects or execute downloaded notebooks/code.
4. Match the target measure to the event PE prior and the injection proposal.
   A four-dimensional(chirp mass,q,chi_eff,z) channel density cannot simply be
   divided by our existing seven-dimensional PE prior. Missing nuisance-spin
   conditionals or different prior support must be explicitly reconstructed
   or kept as an unresolved physical-inference limitation.
5. Before a new fit, specify which physical-library assumptions, missing
   channels, density estimation choices and selection approximations are
   varied. Match event and sensitivity definitions and validate the complete
   conditional likelihood independently. Do not report channel fractions from
   g alone, relabel phenomenological spin components, or reuse a selection
   fraction estimated for another detector configuration without validation.

## Related correlation test

If a mass-dependent-g stress model is added, keep its marginal global-g prior
matched to the constant-g control. Independent Uniform bin fractions induce a
different global prior. Define bin boundaries, mass weights, the global
functional and no-information/nested-control gates before fitting. Such a
stress test is not itself a physical formation-channel library.

## Public scope

Keep raw source files/snapshots local and ignored. Track source URLs, hashes,
license declarations, original adapters, aggregate physical diagnostics and
test results. Existing v0.3 and completed conditional studies stay unchanged.
This phase must advance the actual formation-channel question, rather than
indefinitely adding auxiliary point-estimator audits.
