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

### Scientific population release and documentary excerpts

The stage-one and stage-two reports, scientific SVG figures, diagnostic JSON,
and standalone numerical scripts are this project's analyses. Their source
observations and underlying methods retain the attribution above. Reproduction
instructions and the boundary between public code and separately retrieved data
are in `research/REPRODUCING_POPULATION.md`.

The official BBH injection file used for the stage-two row match is
`o1+o2+o3_bbhpop_real+semianalytic-LIGO-T2100377-v2.hdf5`, published in
https://doi.org/10.5281/zenodo.5636816 and
https://doi.org/10.5281/zenodo.7890398 with the same MD5. It and the direct row-map
NPZ are not redistributed. Download/checksum scripts and aggregate audit results
are supplied. The weighting prescription is credited to Reed Essick (2021),
https://doi.org/10.3847/2515-5172/ac2ba7 .

`research/population-stage2/provenance-history.json` contains short documentary
source-cell excerpts from Colm Talbot's `GWPopulation` tutorial. These excerpts
retain the upstream MIT license below. They are not original project prose, and
the upstream software license does not grant rights to the separately hosted
tutorial samples or LVK data. The checked historical license is:
https://github.com/ColmTalbot/gwpopulation/blob/2daa71e8c56362f2f64e2166a0055478557f239a/LICENSE.md

```text
MIT License

Copyright (c) 2019 Colm Talbot

Permission is hereby granted, free of charge, to any person obtaining a copy
of this software and associated documentation files (the "Software"), to deal
in the Software without restriction, including without limitation the rights
to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is
furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all
copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
SOFTWARE.
```
