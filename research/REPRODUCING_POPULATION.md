# Reproducing the population research

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

## Scope of reproducibility

The calculation, frozen summary input, original scripts, exact input hashes and
public download routes are supplied. Stage one is reproducible from tracked
inputs alone; stage two additionally needs the public tutorial arrays. The
injection-side source rows, finite-file mask and measure conversion have been
recovered. The follow-ups exactly identify two single-waveform PE pilot arrays,
recover the tutorial prior formula on all stored rows, and compare two dependent
representations of one event using their matching prior densities. The remaining
event identities, full PE support, RIFT prior provenance, spin-reference
consistency and selection-consistency checks remain unfinished. Neither these
instructions nor the paper claim that those remaining scientific uncertainties
have been resolved.
