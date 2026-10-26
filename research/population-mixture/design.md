# Official full-O3 mixture follow-up

Recorded 2026-10-27 before reading the mixture HDF or evaluating its population
integrals. This is a source-driven extension of previously observed results,
not a claim of preregistration or a new statistical method.

## Question and fixed inputs

Does using the officially recommended four-population O3 mixture, in place of
the pure BBH O3 bank, materially change the existing conditional all-69-source
results? Which parts of the known missing redshift/spin support does the union
actually cover? The mixture is not an independent observation of sensitivity.

Use only Zenodo record 7890437, `endo3_mixture-LIGO-T2100113-v12.hdf5`,
410957336 bytes, official MD5 `ebbae2546c847a8656022f4c067fb266`.
The fixed upstream record and distribution definition have already been
retrieved and checksummed in population-selection-audit/sources/manifest.json.
Retain all previous results and the published v0.3 artifacts byte-for-byte.

Keep the 69 event integrals, their priors, all six native-tilt models and all
six invariant-effective-spin controls unchanged. Keep the existing early O1/O2
term and O3 analysis-time fraction. Evaluate matched-filter union IFAR>1 yr
without cWB, and the separate union including cWB. Keep their meanings separate.

## Gates and comparisons

1. Verify the full file against official size/MD5 and record SHA256. Inspect
   root/group metadata, population labels, joint sampling_pdf and mixture_weight.
   Do not execute the upstream notebook or equate software defaults to an
   established historical injection waveform/reference.
2. Map every pure-BBH row into the mixture using exact physical coordinates;
   verify uniqueness, all available shared physical/search fields and source
   run times. Report any omitted/changed rows explicitly. The mapping is a
   covariance/source diagnostic, not proof of random independence.
3. Use the joint mixture proposal and official per-row mixture weights in
   the Monte Carlo integral. Do not factor the joint mixture density. Verify
   trial counts, time/trial normalization and coordinate Jacobians. Separately
   show a BBH-only-restricted target control, holding the event likelihood
   fixed; do not call that a globally renormalized new population inference.
4. Report 5/50/95 posterior quantiles for both metrics and policies, median
   shifts, same-functional model contrasts, contributions from newly covered
   support and any concentrated weights. Compare grid doubling with the
   existing absolute 2e-6 quantile tolerance. No chosen pp threshold will turn
   a bank-replacement comparison into a physical robustness certification.
5. Independently recompute at least the mixture-density Jacobian integral,
   selection components, full event-product posterior and quantiles from HDF
   inputs and saved event integrals, without calling the producer's weighting
   or posterior routines. Match relative selection/density at 1e-9, quantiles
   at 2e-6; any failure remains in the record and needs explanation.
6. Derive the union support from the four official mass, spin and redshift
   proposals. Quantify remaining target mass; observed detections or zeros in
   the added region cannot certify absent likelihood/detection contribution
   outside the union or beneath the hopeless cut.

## Deliverables and limits

Keep raw HDF and any direct sample cache local/ignored. Publishable material is
original scripts, source checksums, derived aggregate results, an independently
computed validation record, a scientific figure and a concise bilingual account.
Add to the research roadmap after the checks; do not silently replace v0.3.
No new physical formation-channel mapping, historical pipeline execution,
campaign-level coverage or common-reference tilt claim follows from this step.
