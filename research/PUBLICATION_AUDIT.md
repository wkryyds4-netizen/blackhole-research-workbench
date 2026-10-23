# Scientific release audit

Date: 2026-10-24. This is the working-tree publication audit for the population
research release. It does not certify the entire Git history or future additions.
The release coordinator remains responsible for checking the final staged tree.

## Files suitable for this release

| Content | Paths |
|---|---|
| Scientific roadmap and literature/provenance analysis | `research/SCIENTIFIC_ROADMAP.md`, `research/population-audit/*.md`, authored audit JSON |
| Original numerical methods and independent checks | `research/population_*.py`, `research/test_population_functional_math.py`, `research/validate_population_*.py`, `research/audit_population_inputs.py` |
| Original provenance acquisition/check scripts | `research/population-audit/*.py`, `research/population-stage2/provenance_*.py` |
| Stage-one calculated summaries and figure | `research/population-stage1/*.json`, `README.md`, `prior-support-control.svg` |
| Stage-two calculated summaries and figure | `research/population-stage2/*.json`, `*.md`, `paired-mc-comparison.svg` |
| Small authored source manifests/range audit | `population-audit/provenance-sources/manifest.json`, `population-stage2/provenance-sources/fetch-manifest.json`, `remote-headers.json` |
| Reproduction dependencies and documentation | `research/population-requirements.txt`, `research/REPRODUCING_POPULATION.md`, this report and `publication-audit.json` |
| Notices | updated root `NOTICE.md`, `CONTRIBUTING.md` |

The paper/source manuscript and any additional final artifacts must be added to
the coordinator's final scan; they were being prepared concurrently. These path
patterns are an allowlist by purpose, not an instruction to add every future
file matching a wildcard.

## Exclusions

- Raw posterior/injection pickle, HDF/HDF5 and direct NPZ array conversions,
  including the 63,149,928-byte official injection HDF and source-row map.
- The local bootstrap mean cache, environments, package caches, logs, private
  configuration, access credentials and raw working directories.
- Full fetched third-party notebooks, HTML source, GitHub API history/tree and
  Zenodo record snapshots. Public links, hashes and fetch scripts provide the
  provenance without needing those documentary copies in the repository.
- SXS2327 raw data and direct derived trajectories, as already excluded by the
  previous public release; no new redistribution license has been established.

The injection HDF directory has a local `.gitignore` for `*.hdf5` and `*.part`;
`git check-ignore` confirmed that the downloaded file is excluded. NPZ/pickle
data and `/data/` are covered by root ignore rules. The release coordinator was
asked to additionally exclude the new `.population-venv/`, smoke-output directory
and bulk stage-two source snapshots while preserving the authored manifests.

## Attribution review

`provenance-history.json` contains limited verbatim tutorial source cells. The
fixed 2024 GWPopulation repository license was retrieved and checked: MIT,
Copyright (c) 2019 Colm Talbot. `NOTICE.md` now identifies those excerpts and
includes that full notice. This retains their upstream provenance without
changing or relabeling the research output. The tutorial code license is not
asserted to license the separately hosted data files.

All raw scientific inputs remain separate downloads. The scientific SVGs are
original plots of this project's calculated diagnostics. Reports distinguish
the source data, cited numerical methods and this project's implementation;
they do not claim LVK endorsement or a novel formation-channel measurement.

## Local path and credential check

The heuristic scan considered all tracked working-tree paths and nonignored
untracked files at scan time: **241 text files**, no binary files and no file
larger than 5 MB. It found no credential-like tokens, private-key markers,
credential assignments or personal absolute filesystem paths. Only filenames,
categories and counts were retained; candidate credential values would not
have been printed. The exact scope and limitations are in
[publication-audit.json](publication-audit.json).

An earlier broad path regex matched public `https://api.github.com/users/...`
URLs in a fetched GitHub snapshot. Inspection established these were public
account API links, not local user paths. That false positive was removed by
distinguishing URL host paths from filesystem paths.

This scan does not prove the absence of all secrets. It also does not inspect
old commits, ignored local files or artifacts created after the scan.

## Reproduction review

The existing public repository already tracks `download_data.py`,
`prepare_data.py`, `model.py`, `validate_model.py` and the frozen
`results/model.json`. The raw acquisition path is therefore public. The new
[reproduction guide](REPRODUCING_POPULATION.md) supplies:

1. Offline stage-one and synthetic checks using tracked inputs.
2. Exact tutorial input SHA256 checks before the restricted pickle conversion.
3. Full stage-two commands, seeds and repetition counts.
4. Official HDF retrieval with MD5 validation and reconstruction of the excluded
   row map before the mixture-weight control.
5. Explicit limits: externally hosted tutorial objects may become unavailable,
   snapshot metadata can change, and event-level PE provenance remains incomplete.

The root requirements alone omit matplotlib and h5py. The separate
`population-requirements.txt` adds the observed versions and inherits the root
requirements. The main science environment was Python 3.14.5; provenance HDF
checks used Python 3.12 with h5py 3.16.0. A cross-platform fresh installation of
the combined requirements has not been claimed as tested in this audit.

No scientific numerical result, original input, model parameter or hypothesis
was changed for this publication audit. No files were staged, committed or
pushed by this audit task.
