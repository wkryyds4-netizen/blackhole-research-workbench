"""Record what the frozen injection HDF and official definition actually expose."""
from pathlib import Path
import hashlib
import json
import h5py

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]


def main():
    original=json.loads((ROOT/'research/population-stage2/provenance-injection-match.json').read_text())
    path=ROOT/'research/population-stage2/source-data/o1+o2+o3_bbhpop_real+semianalytic-LIGO-T2100377-v2.hdf5'
    assert hashlib.sha256(path.read_bytes()).hexdigest()==original['source_sha256']
    manifest=json.loads((BASE/'sources/manifest.json').read_text())[0]
    doc=BASE/'sources'/manifest['file'];assert hashlib.sha256(doc.read_bytes()).hexdigest()==manifest['sha256']
    with h5py.File(path,'r') as h:
        g=h['injections'];fields=list(g)
        attributes={name:{k:float(v) for k,v in group.attrs.items()} for name,group in [('root',h),('injections',g)]}
        reference_fields=[k for k in fields if any(token in k.lower() for token in ['freq','wave','approxim','reference'])]
        result=dict(schema='blackhole-injection-metadata-followup/1',source_sha256=original['source_sha256'],
            official_definition_sha256=manifest['sha256'],attributes=attributes,fields=fields,
            candidate_reference_or_waveform_fields=reference_fields,
            hdf_spin_reference_frequency_established=False,
            hdf_complete_trial_design_established=False,
            interpretation='The frozen summary HDF has no frequency/waveform/reference fields. The official definition identifies semianalytic O1/O2 and real-search O3, but explicitly omits individual draw-distribution details. This is evidence of limits in these two artifacts, not proof that no other public record exists.')
    (BASE/'metadata-audit.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='fields'},indent=2))


if __name__=='__main__':main()
