"""Check restored rows against recovered LALInference bounds and coordinate density."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
sys.path.insert(0,str(ROOT/'research/population-stage4'))
from recover_prior import LAL


def main():
    recovery=json.loads((BASE/'configuration-prior-recovery.json').read_text());provenance=json.loads((BASE/'event-provenance.json').read_text())
    output=BASE/'configuration-prior-row-audit.json';rows=[];rng=np.random.default_rng(20261025)
    for event in recovery['events']:
        info=next((r for r in provenance if r['index']==event['index'] and r.get('native_source_identity_resolved')),None)
        if info is None:continue
        path=BASE/f'event-{event["index"]}.npz';assert hashlib.sha256(path.read_bytes()).hexdigest()==info['cache_sha256']
        with np.load(path,allow_pickle=False) as bank:d={k:bank[k] for k in bank.files}
        violations={};ranges={}
        for name,bounds in event['bounds'].items():
            key='mass_1_detector' if name=='mass_1' else name
            values=d[key] if key in d else np.arccos(d['cos_'+name])
            violations[name]=int(np.count_nonzero((values<bounds['minimum']-1e-9)|(values>bounds['maximum']+1e-9)))
            ranges[name]=[float(values.min()),float(values.max())]
        assert not any(violations.values()),(event['event'],violations)
        # Numerical Jacobian of source coordinates -> historical sampling coordinates.
        def coordinates(v):
            m,q,z,c1,c2=v;mc=m*(1+z)*q**.6/(1+q)**.2
            return np.array([mc,q,np.log(LAL.luminosity_distance(z).value),np.arccos(c1),np.arccos(c2)])
        errors=[]
        for i in rng.choice(len(d['mass_1']),32,replace=False):
            v=np.array([d['mass_1'][i],d['mass_ratio'][i],d['redshift'][i],d['cos_tilt_1'][i],d['cos_tilt_2'][i]])
            steps=1e-5*np.maximum(abs(v),.01);steps[2]=min(steps[2],v[2]/10)
            steps[3:]=np.minimum(steps[3:],(1-abs(v[3:]))/10)
            jac=np.empty((5,5))
            for j,h in enumerate(steps):
                up=v.copy();down=v.copy();up[j]+=h;down[j]-=h;jac[:,j]=(coordinates(up)-coordinates(down))/(2*h)
            mc,q,logdl,t1,t2=coordinates(v);dl=np.exp(logdl);mdet=v[0]*(1+v[2])
            historical=(mdet**2/mc)*dl**3*np.sin(t1)*np.sin(t2)
            ddz=dl/(1+v[2])+(1+v[2])*LAL.hubble_distance.value/LAL.efunc(v[2])
            direct=dl**2*ddz*v[0]*(1+v[2])**2
            errors.append(abs(historical*abs(np.linalg.det(jac))/direct-1))
        error=float(max(errors));assert error<1e-5
        rows.append(dict(index=event['index'],event=event['event'],rows=len(d['mass_1']),bounds_violations=violations,observed_ranges=ranges,
            numerical_sampling_coordinate_jacobian_max_relative_error=error,source_cache_sha256=info['cache_sha256'],
            source_supported_default_family=True,scope='Released configuration plus documented defaults; rows satisfy recorded support. No bound on likelihood outside it.'))
    report=dict(schema='blackhole-configuration-prior-row-audit/1',recovery_sha256=hashlib.sha256((BASE/'configuration-prior-recovery.json').read_bytes()).hexdigest(),events=rows,
        pending=[e['index'] for e in recovery['events'] if e['index'] not in {r['index'] for r in rows}],
        density_comparison='Historical (Mchirp,q,logDL,tilt1,tilt2) branch times numerical five-dimensional Jacobian versus direct source-coordinate density, up to fixed normalizations.')
    output.write_text(json.dumps(report,indent=2),encoding='utf-8');print(json.dumps(dict(completed=[r['index'] for r in rows],pending=report['pending'],errors=[r['numerical_sampling_coordinate_jacobian_max_relative_error'] for r in rows]),indent=2))


if __name__=='__main__':main()
