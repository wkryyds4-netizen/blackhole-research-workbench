# Code and data attribution

The root MIT license covers this project's original code and documentation.
It does not relicense third-party code, scientific datasets, model weights, or
figures. Original notices remain authoritative.

## Rendering and scientific software

- The renderer adapts **algometrix/blackhole-sim**, MIT, and uses **Three.js**, MIT.
  See `simulator/THIRD_PARTY.md`, `simulator/renderer-src/vendor/blackhole-sim/`,
  and `simulator/public/licenses/THIRD-PARTY.txt` for notices and provenance.
- **GWSurrogate** uses the MIT/X Consortium license. The portable Windows build
  and compatibility patch are documented in `research/surrogate-results/`.
  The build cache and third-party source archive are not redistributed here.
- **Ripple**, JAX, NumPy, SciPy, SXS and other installed packages retain their
  own licenses. Python environments and downloaded packages are excluded.

## Numerical relativity and surrogate results

SXS Collaboration's historical **0001, 0003 and 0389** snapshots are CC BY 4.0.
Sources, contributors, DOI links and transformations are listed in
`simulator/public/licenses/SXS-DATA.txt`. These JSON files are derived snapshots;
0003 is explicitly retained as a deprecated historical demonstration.

NRSur7dq4 model attribution, model DOI and CC BY 4.0 terms are recorded in
`simulator/public/licenses/SURROGATE-DATA.txt`. Published waveform and dynamics
JSON retain provenance and are derived from this model, rather than observations.

**SXS2327 is excluded from this public distribution** because its source record
did not declare a redistribution license at retrieval. Local research reports,
URLs and file hashes do not constitute permission to redistribute its samples.
Do not add locally downloaded 2327 data or direct derived trajectories to a PR.

## Observations and population analysis

GW150914 detector data come from **GWOSC / LIGO Scientific Collaboration and
Virgo Collaboration**. GWOSC states CC BY 4.0 for its data:
https://gwosc.org/about/ . Event source:
https://gwosc.org/eventapi/html/GWTC-1-confident/GW150914/v3/ .
The local preparation script whitens, filters, crops and resamples the data;
small processed event JSON and download manifests describe those changes.

Population analysis follows **Colm Talbot's GWPopulation GWTC-3 tutorial**:
https://colmtalbot.github.io/gwpopulation/examples/GWTC3.html .
The tutorial uses preprocessed LVK posterior samples and search injections from:

- https://doi.org/10.5281/zenodo.5546663
- https://doi.org/10.5281/zenodo.6513631
- https://doi.org/10.5281/zenodo.7890398

This project computes additional conditional models, importance weights,
selection corrections and numerical diagnostics. Downloaded tutorial sample
files and direct NPZ conversions are not redistributed. The tutorial software
license is not asserted to license those sample files. Public summary outputs
are exploratory calculations, not an official LVK result or formation-channel
measurement. No endorsement by the cited authors or collaborations is implied.
