# Reproducing the population research

Current snapshot (2026-10-25; research paused at the user's request): this guide also documents the completed 69-source
restoration, invariant effective-spin control, paired O3/event/joint finite-bank
diagnostics and conditional missing-support follow-up. The v0.2 paper/PDF and
historical workbench view retain earlier frozen results; the separate v0.3
LaTeX/PDF revision and initial population-view adapter incorporate the all-69 follow-ups. Earlier source-recovery counts
below are historical snapshots; the current catalog note and final commands in
this guide take precedence for the all69 follow-up. The adapter displays saved results only; its complete interaction and responsive
layout verification remains unfinished. This is an interim research snapshot.

This guide covers the two scientific stages in paper v0.2 and the separately
published stage-three through stage-five follow-up studies. The result is a
conditional, phenomenological spin-tilt analysis of a tutorial dataset, not an
LVK reproduction or a measurement of astrophysical formation-channel fractions.
See [stage one](population-stage1/README.md), [stage two](population-stage2/README.md)
and the [literature comparison](population-audit/literature-overlap.md).

Post-release follow-up: the [stage-three provenance pilot](population-stage3/README.md)
confirms the proposed 69-event identifier set and exactly matches two events to
single-waveform IMRPhenomXPHM samples. Their configurations record different
spin reference frequencies (10 and 20 Hz). This adds an explicit limitation to
physical interpretation; it does not change the frozen v0.2 numerical products.

The [stage-four follow-up](population-stage4/README.md) reconstructs all 504,649
stored prior values, runs a two-event conditional PE replacement test, and
identifies an influential source draw absent from a tutorial subset. Its raw
pilot files and scripts are separate from the frozen v0.2 reproduction path.

The [stage-five single-event study](population-stage5/README.md) matches the
cosmologically reweighted release back to 169,040 original-prior draws. It
checks prior-measure conversion and shows reduced single-draw influence when
the original draws are retained with their weights. This is a dependent-data
representation test, not an independent astrophysical observation.

The [stage-six follow-up](population-stage6/README.md) adds a weighted-event
interface, verifies it against the stage-five source rows, and screens all 69
stored events under six fixed settings. It separates representation-weight ESS,
population-integral contribution ESS, single-row influence and leave-event-out
influence. It identifies four follow-up candidates; their original-prior files
have only been located in source metadata, not yet downloaded or verified.

Run commands from the repository root. Use a fresh clone or copy if you want to
preserve the published outputs: the numerical scripts write their result files
and validation reports to the stage directories by default. Do not run `model.py`
as a prerequisite; `results/model.json` is the frozen, tracked input for these
stages, and `model.py` would overwrite it.

## Environment

The main calculations were run with Python 3.14.5, NumPy 2.4.6, SciPy 1.18.0,
pandas 2.3.3, Astropy 8.0.1, requests 2.34.2 and matplotlib 3.11.0. The separate
HDF provenance audit used h5py 3.16.0 in a Python 3.12 environment. The file below
collects the direct dependencies; it is not a complete transitive lockfile, and
this release does not claim a fresh installation of that combined environment
has been validated on every platform.

```text
python -m venv .population-venv
```

Activate that environment using the standard command for your shell, then run:

```text
python -m pip install -r research/population-requirements.txt
```

GWPopulation, Bilby, GWSurrogate and the visualization compiler are not needed to
run these standalone population calculations. Their methods and public tutorial
are credited in [NOTICE.md](../NOTICE.md). Installation of the browser workbench
is separate from this numerical workflow.

## 1. Offline stage-one reproduction and synthetic validation

The tracked `results/model.json` contains the original weight moments. Its
SHA256 is `6dcf3c9f4fd51b2007977d131c61f69c7d861afb986551ea82b70a7667a51083`.
These commands require no raw-data download:

```text
python -X utf8 research/population_functionals.py
python -X utf8 research/test_population_functional_math.py
python -X utf8 research/validate_population_functionals.py
python -X utf8 research/validate_population_mc_synthetic.py
python -X utf8 research/validate_population_mc_implementation.py
```

The first three reproduce and independently check the common-prior functional
analysis. The last two validate the empirical resampling mechanism and numerical
implementation on small synthetic inputs; they do not validate the real data
selection or population-model correctness. SVG figures are tracked. PNG files
are local rendering previews.

## 2. Restore the public tutorial inputs for stage two

The original download and restricted-conversion scripts are already public:
[download_data.py](../download_data.py) and [prepare_data.py](../prepare_data.py).
`download_data.py` retrieves the three GWOSC catalog JSON files and two tutorial
files, about 39 MB in total. It does not install packages or execute downloaded
code. The URLs are recorded in the script and original provenance manifests.

```text
python -X utf8 download_data.py
```

Before conversion, verify the two downloaded tutorial objects against the frozen
research inputs. Run this Python snippet from the repository root:

```python
from pathlib import Path
from hashlib import sha256

expected = {
    "gwtc-3-samples.pkl": "f0e2554a595630a4bea9c7db58dc9f83bef036ae800d4bc318cdd676655ba2dd",
    "gwtc-3-injections.pkl": "99b1ba450fd554f2d4b992c9eb45d389f460d224f5684691b114c59ed8d26acc",
}
for name, checksum in expected.items():
    assert sha256((Path("data") / name).read_bytes()).hexdigest() == checksum, name
```

If a checksum differs, stop this frozen-input reproduction and retain the changed
file as a separate data version. The tutorial Drive links are public but their
continued availability is controlled by their owner. Do not substitute arbitrary
pickle files. The conversion script permits only a small set of NumPy/pandas
reconstruction types; it does not use unrestricted `pickle.load`.

```text
python -X utf8 prepare_data.py
```

This produces `data/posteriors.npz`, `data/injections.npz` and
`data/event_summary.json`. The published NPZ SHA256 values are:

| Input | SHA256 |
|---|---|
| `data/posteriors.npz` | `3aec709780ac373e485682323c20cbe6ed0c14028c66b030bb73d848739489c4` |
| `data/injections.npz` | `00ce46283aaf5c1e3df1a8b479a2f3c5bf904e4aef6aa0cf641f93df54bf7a13` |

The stage-two script recomputes raw weights and verifies their moments against
the frozen model. Archive bytes can differ across dependency/serialization
versions even when arrays agree; do not silently replace a checksum in a frozen
result. Investigate such differences and record a new run separately.

## 3. Run the paired empirical diagnostic

```text
python -X utf8 research/population_mc_bootstrap.py
```

Defaults are two fixed seeds, 1024 repetitions each, and 1001 r-grid points. This
writes `paired-bootstrap-results.json`, the SVG/PNG comparison figure and a local
`bootstrap-means.npz` cache. It does not change raw inputs or the frozen model.
For a short implementation check, a separate output directory may be used:

```text
python -X utf8 research/population_mc_bootstrap.py --repeats-per-seed 16 --output research/local-population-smoke
```

A short run does not reproduce the published 2048-repetition results. Fixed RNG
seeds do not promise byte-identical outputs across all libraries and hardware;
reported covariance and quadrature checks specify the numerical tolerances.
The bootstrap repeats an explicitly defined empirical mechanism with an omitted
zero-weight bucket. It does not reproduce the original injection campaign's full
sampling design or add independent observations.

## 4. Rebuild the injection provenance and mixture-weight control

The HDF, local row map and third-party source snapshots are intentionally not
distributed. All can be retrieved or rebuilt using public URLs:

```text
python -X utf8 research/population-audit/provenance_stage2_fetch.py
python -X utf8 research/population-stage2/provenance_download_bbh.py
python -X utf8 research/population-stage2/provenance_match_injections.py
python -X utf8 research/population_mixture_weight_control.py
```

The first command retrieves small official metadata and tutorial documents.
Check its manifest for failed downloads; retry before continuing if the required
`zenodo-7890398.json` metadata is absent. Network/TLS errors were observed during
the original retrieval and are recorded rather than treated as valid evidence.
The HDF download then retrieves exactly one 63,149,928-byte public BBH file and
verifies the official MD5 `c21337a25b7318fdfd471cce549a547e` before finalizing it.

The match script independently checks the HDF checksum, unique source-row
identities, seven exact physical-coordinate arrays, the specified detection mask
and the draw-density Jacobian. It rebuilds
`research/population-stage2/provenance-injection-row-map.npz` with keys
`original_row_indices`, `run`, and `mixture_weight`. The last command uses that
map to reproduce the deterministic mixture-weight and normalization controls.

Optional documentation-history audit:

```text
python -X utf8 research/population-stage2/provenance_history_audit.py
```

It needs both fetched notebooks and both Zenodo metadata files. The original
notebook is commit-pinned; the current notebook and API history can change after
this release. Those documentary changes are separate from the HDF numerical
match. See the published [provenance report](population-stage2/provenance-README.md)
for the exact evidence boundary.

## 5. Follow-up stages three through five (2026-10-25)

These follow-ups do not regenerate the paper or overwrite `results/model.json`.
Restore the tutorial arrays using section 2, then run each stage's commands in
the order listed in its README:

1. [Stage three](population-stage3/README.md): fetch documentary sources and two
   official cosmo HDF files (64,794,885 bytes total), audit event membership,
   match the two PE labels and source rows, validate, and plot. The row-map NPZ
   is generated locally and is required by stage four.
2. [Stage four](population-stage4/README.md): fetch pinned source snapshots,
   reconstruct the stored prior, extract the pilot bank, run the six-setting
   sensitivity study, validate, and plot. This produces the local pilot-bank
   and checkpoint NPZ files needed by the follow-up calculations.
3. [Stage five](population-stage5/README.md): download the fixed original-prior
   HDF (335,235,368 bytes) and methods papers, extract/match representations,
   run the tail study, validate, and plot. Continuous-weight identities and
   independent full-product checks are included.

Use the requirements above; HDF extraction needs h5py and the prior/tail
calculations need Astropy. The original runs used separate Python environments
as recorded in each report. Raw HDFs, NPZ intermediates, PDF/HTML/source snapshots,
and PNG previews are ignored by Git; scripts, source URLs and hashes, derived
JSON summaries and scientific SVGs are public. The validators need the restored
local inputs, not only the tracked JSON summaries. A clean-machine end-to-end
download and installation has not been certified by this publication.

## 6. Weighted-event interface and catalog screening

After restoring the stage-three and stage-five local inputs, run the six
commands in the [stage-six report](population-stage6/README.md). They are
offline; only `validate_source_rows.py` needs h5py, while the numerical interface
and screen use NumPy, SciPy and Astropy. Intermediate NPZ files are local caches.
The source-row validator must run after `validate_interface.py` creates its
weighted-event cache. This stage does not execute a new posterior sampler or
download the four candidates listed in its follow-up source plan.

## Scope of reproducibility

### Follow-up stages seven and eight

The [stage-seven report](population-stage7/README.md) lists the ordered commands
for four additional event releases (eight HDF files, 1,658,552,960 bytes), their
exact row matches, joint replacements, independent verification and figure.
It requires the stage-five and stage-six caches. The
[selection study](population-selection-audit/README.md) additionally restores
the fixed full-O3 BBH injection file (88,842,240 bytes), verifies its exact link
to the combined-file subset, and propagates a run-weighted replacement and a
joint event/selection redshift-support restriction. A separate 512-replicate
benchmark calibrates only an explicitly specified artificial thinning design
conditional on that saved injection bank. Neither study overwrites the frozen
paper/results. Follow each report's commands; raw sources and NPZ caches are
local prerequisites, not tracked release files.

The [spin-reference audit](population-reference/README.md) restores azimuths and
released infinity tilts from five original-prior event files. Its precession
calculations require an isolated Python 3.12 runtime with the pinned NumPy,
SciPy and precession versions. Follow the extraction, diagnosis, held-out
validation and near-equal-mass refinement commands in order; retain both the
initial and refined diagnostics. The study verifies a numerical convention on
430 distinct rows, not a common-reference full-catalog population analysis.
The [search-policy control](population-selection-audit/pipeline-policy.md) then
uses the existing O3 bank and checkpoints to compare the tutorial search union
with the matched-filter-only rule described in LVK Appendix A. Both studies
preserve the frozen paper and original three input/result hashes.

The calculation, frozen summary input, original scripts, exact input hashes and
public download routes are supplied. Stage one is reproducible from tracked
inputs alone; stage two additionally needs the public tutorial arrays. The
injection-side source rows, finite-file mask and measure conversion have been
recovered. The follow-ups exactly identify two single-waveform PE pilot arrays,
recover the tutorial prior formula on all stored rows, and compare dependent
representations of five events using their matching prior densities. The remaining
event identities, full PE support, RIFT prior provenance, spin-reference
consistency and selection-consistency checks remain unfinished. Neither these
instructions nor the paper claim that those remaining scientific uncertainties
have been resolved.

The latest local [catalog expansion](population-catalog/README.md) inventories
reference frequencies, spin fields and analytic prior declarations in all 69
fixed candidate files. Its saved 20-source snapshot replaces 15 additional
full original-prior event posteriors beyond the five-source baseline, retaining
49 tutorial representations in the fixed 69-event likelihood. Six missing
analytic-prior strings are investigated using released LALInference
configurations and pinned historical defaults; full-row checks are incremental.
It also records a label-dependent infinity-tilt discrepancy, the matching
historical PESummary source, and a retained failed 1e-6 release-agreement gate.
Historical Python-default replay on the same 70 rows explains the small
numerical remainder, without reclassifying the original gate or claiming a
reconstructed historical runtime.
Run its source, extraction, evolution, comparison and independent-validation
commands in the documented environments. Metadata availability is complete;
the remaining source connections and common-reference selection
protocol are not. Raw files, range caches and third-party source snapshots
remain local; public retrieval URLs and checksums are recorded.

The completed69-source expansion contains11,180,578 original-prior posterior
rows. Follow both label-resolution commands before extraction: index48 needs
the previously checksum-verified cosmo pilot file and a three-way connection;
index72 needs an exact decomposition of its Mixed file. Containing labels need
not be unique, and Mixed is not silently treated as Phenom. All six recovered
configuration-prior families now have full-row range checks. The full native
source replacement remains distinct from a common-reference tilt analysis.

The supplementary control in `population-invariant/` requires this completed,
independently validated69-source gate. In order, run `validate_angular.py`,
`calibrate_sign_experiment.py`, `run_control.py`, `validate_control.py`, and
`plot_control.py` in the main NumPy/SciPy/Astropy environment. Its pre-recorded
design is `design.json`; it measures the binary negative-effective-spin
fraction with an explicit2PN/isotropic-measure approximation. The known
selected/noisy-sign experiment is implementation calibration, not a real-GW
coverage study. Do not reuse the older five-event thinning intervals as
uncertainty for either new69-event result.

The all69 paired finite-bank follow-up is now complete. After the validated
full-source invariant control, run the following in order in the main
NumPy/SciPy/Astropy environment (Matplotlib is needed for plots):

```text
python research/population-invariant/paired_o3.py
python research/population-invariant/validate_paired_o3.py
python research/population-invariant/paired_pe.py
python research/population-invariant/diagnose_pe_derivative.py
python research/population-invariant/validate_paired_pe.py
python research/population-invariant/plot_finite_banks.py
python research/population-invariant/support_sensitivity.py
python research/population-invariant/validate_support.py
python research/population-invariant/plot_support.py
python research/population-invariant/plan_tail_validation.py
python research/population-invariant/summarize_followup.py
```

These scripts do not execute root `model.py`; only established functions are
imported. Do not run root `model.py`, which would replace frozen early inputs
and results. Do not modify `run_control.py` casually: its bytes are part of
the saved event-integral protocol. No source download or all-event refit is
necessary to reproduce support sensitivity from the already validated control
checkpoints.

The paired event experiment retains6913 rows per source and512 replicates;
it keeps69 observations. Shared O3 draws retain128454 of284672 accepted rows,
with O1/O2 fixed. The retained original derivative-step failure is documented
in `pe-derivative-diagnosis.json`: an unchanged relative gate passes at a
larger real step and independent complex derivatives validate the weak score.
Do not relabel the original step as passing or treat this as a modification
of production posteriors. Both paired validators use independent products,
source calculations and finite-population covariance formulas.

The support follow-up was designed after observing a zero event component.
It derives conditional CDF/quantile envelopes for an *assumed* log-correction
budget, and evaluates nine declared completions of one event integral.
It cannot estimate missed likelihood from a finite-bank zero. Its omitted-O3
domain bound uses only detection probability<=1; it is intentionally too broad
to certify negligible missing selection. `tail-validation-plan.json` is an
unexecuted direct-draw experiment/resource plan, not an injection dataset.
Its binomial bounds require independent direct target draws and a validated
historical detection law; weighted bank samples do not automatically qualify.

All three full-source frozen hashes remain unchanged. Older manuscript/PDF
v0.2 versions describe earlier snapshots and must not be cited as the report of
these all69 follow-ups. The new v0.3 manuscripts are separately versioned.
Current scientific results and their limits are
recorded in [the catalog research note](population-catalog/README.md).
The final summary command verifies all result/validation digest links, the
unchanged event-integral protocol, the three frozen hashes and preservation
of the original derivative-step failure; it writes `finite-bank-milestone.json`.

The interim bilingual paper revision is in `paper/v0.3/latex-en` and
`paper/v0.3/latex-zh`, each with its own complete source tree and five vector
figure PDFs. Both 14-page exports are in `output/pdf/`. From the repository root,
`python research/paper/v0.3/build_evidence.py` regenerates shared facts/tables
and publication figures from validated results, while
`python research/paper/v0.3/compile_paper.py --render` compiles both projects
and renders all pages with Poppler. The PDF environment additionally needs
pypdf, pdfplumber and Pillow for `inspect_paper.py`. After actual visual and
bilingual self-review, `validate_release.py` checks the new artifact chain;
it does not recompile or independently review the science. Rebuilds invalidate
old PDF/review hashes and require inspection again. The earlier finite-bank
milestone retains the status at its creation; the v0.3 build/review manifests
record this later manuscript step. No current-cycle Git push is implied.

Regenerating publication figures also requires the validated local
`population-invariant/control-checkpoints.npz` and `support-checkpoints.npz`,
plus the frozen tutorial NPZ inputs checked by the evidence builder. These
derived caches are rebuilt by the documented research pipeline and are not
distributed in Git. Compiling the already generated self-contained LaTeX
projects with their supplied figure PDFs does not require those raw/cache
files; the Python build wrapper additionally checks tracked research JSONs.
