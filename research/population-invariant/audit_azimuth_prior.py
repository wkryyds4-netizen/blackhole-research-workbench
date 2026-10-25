"""Check the azimuth-prior premise separately from PN measure preservation."""
from pathlib import Path
import hashlib
import json
import re
import numpy as np

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
CAT=ROOT/'research/population-catalog'
SOURCES=ROOT/'research/population-selection-audit/sources'


def main():
    metadata=json.loads((CAT/'metadata-results.json').read_text());rows=[];missing=[]
    for event in metadata:
        if event.get('mode')!='remote':continue
        analysis=next(r for r in event['labels'] if r['label']=='C01:IMRPhenomXPHM')
        declaration=analysis.get('analytic_prior_strings',{}).get('phi_12')
        if not declaration:missing.append(event['index']);continue
        text=declaration[0];bounds={k:float(v) for k,v in re.findall(r'(minimum|maximum)=([\d.eE+\-]+)',text)}
        assert text.startswith('Uniform(') and bounds['minimum']==0 and abs(bounds['maximum']-2*np.pi)<1e-12
        rows.append(dict(index=event['index'],mode='released_analytic_uniform_azimuth',declaration=text))
    configs=json.loads((CAT/'missing-prior-configurations.json').read_text())
    assert set(missing)=={r['index'] for r in configs}
    manifest=json.loads((SOURCES/'manifest.json').read_text())
    name='gitlab-1438-af4bc21a9f92-lalinference--lib--LALInferenceInitCBC.c';source=SOURCES/name
    raw=source.read_bytes();entry=next(r for r in manifest if r['file']==name)
    assert hashlib.sha256(raw).hexdigest()==entry['sha256']
    text=raw.decode();assert 'REAL8 phi12min=0.0,phi12max=LAL_TWOPI;' in text
    assert 'LALInferenceRegisterUniformVariableREAL8(state, model->params, "phi12", 0.0, phi12min,phi12max,phi12Vary);' in text
    for row in configs:
        flags=[p for p in row['configuration_paths'] if p.startswith('engine/') and any(t in p.lower() for t in ['phi','azimuth','fix'])]
        assert not flags,(row['index'],flags)
        rows.append(dict(index=row['index'],mode='released_configuration_and_historical_uniform_default',configured_azimuth_or_fixed_engine_overrides=flags))
    assert len(rows)==69
    report=dict(schema='blackhole-relative-azimuth-prior-audit/1',analyses=69,analytic_declarations=63,configuration_default_families=6,rows=rows,
        default_source=dict(file=name,url=entry['url'],sha256=entry['sha256']),
        scope='Relative-azimuth prior declarations/default family for each selected Phenom-family analysis. Does not prove PN preservation of isotropic measure, identify unrecorded runtime overrides, or recover injection reference metadata.')
    (BASE/'azimuth-prior-audit.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:report[k] for k in ['analyses','analytic_declarations','configuration_default_families']},indent=2))


if __name__=='__main__':main()
