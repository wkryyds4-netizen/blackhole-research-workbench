"""Controlled initial-L variants; exploratory diagnosis, not a fitted correction.

PN angular momentum coefficients: LALSuite 7.26 documentation, citing
Bohé et al. 2013 arXiv:1212.5520 Eq. 4.7 and Kidder 1995 arXiv:gr-qc/9506022.
Independent formula evaluation; no upstream source code is vendored/executed.
"""
from pathlib import Path
import json
import hashlib
import warnings
import numpy as np
from transport import to_infinity

BASE=Path(__file__).resolve().parent
MTSUN=4.9254909476412675e-6  # SI seconds per solar mass (LAL convention)


def initial_state(inputs,kind):
    k=dict(inputs);q=k['q'];x=np.array([1.,q])/(1+q);eta=x.prod()
    v=(np.pi*MTSUN*k['total_detector_mass']*k['reference_hz'])**(1/3)
    t1,t2,phi=k['theta1'],k['theta2'],k['phi12']
    unit=np.array([[np.sin(t1),0,np.cos(t1)],
                   [np.sin(t2)*np.cos(phi),np.sin(t2)*np.sin(phi),np.cos(t2)]])
    spin=unit*(x*x*np.array([k['a1'],k['a2']]))[:,None]
    correction=1.
    if kind!='newtonian_frequency':correction+=(1.5+eta/6)*v*v
    if kind in ('L_2PN_nonspin','L_2PN_spin_norm','L_2PN_spin_vector'):
        correction+=(27/8-19*eta/8+eta*eta/24)*v**4
    direction=np.array([0.,0.,correction])
    if kind in ('L_2PN_spin_norm','L_2PN_spin_vector'):
        direction-=v**3*np.sum(((1+3/x)/4)[:,None]*spin,axis=0)
        direction[2]-=v**3*np.sum((1+27/x)*spin[:,2])/12
    L=eta/v*np.linalg.norm(direction)
    if kind=='L_2PN_spin_vector':
        axis=direction/np.linalg.norm(direction)
        cosine=unit@axis
        inplane=unit-cosine[:,None]*axis
        cosphi=np.dot(inplane[0],inplane[1])/np.prod(np.linalg.norm(inplane,axis=1))
        k.update(theta1=float(np.arccos(np.clip(cosine[0],-1,1))),
                 theta2=float(np.arccos(np.clip(cosine[1],-1,1))),
                 phi12=float(np.arccos(np.clip(cosphi,-1,1))))
    k['initial_separation']=float((L/eta)**2)
    return k


def main():
    source_path=BASE/'reference-results.json'
    source=json.loads(source_path.read_text());rows=[];failures=[];warning_log=[]
    kinds=['newtonian_frequency','L_1PN_nonspin','L_2PN_nonspin','L_2PN_spin_norm','L_2PN_spin_vector']
    for n,record in enumerate(source['transport_comparison']):
        for kind in kinds:
            try:
                with warnings.catch_warnings(record=True) as caught:
                    warnings.simplefilter('always')
                    value=to_infinity(**initial_state(record['inputs'],kind),solver='DOP853',rtol=1e-11,atol=1e-13)
                meaningful=sorted({str(w.message) for w in caught if not issubclass(w.category,DeprecationWarning)})
                if meaningful:warning_log.append(dict(index=record['index'],source_row=record['source_row'],kind=kind,warnings=meaningful))
                value.update(index=record['index'],source_row=record['source_row'],kind=kind,
                    difference=[value[f'cos_tilt_{j}_infinity']-record['released_cosines'][j-1] for j in [1,2]],
                    chieff_difference_from_native=value['chieff']-record['chieff'])
                rows.append(value)
            except Exception as exc:
                failures.append(dict(index=record['index'],source_row=record['source_row'],kind=kind,error=str(exc)))
        if (n+1)%16==0:print(f'{n+1}/80 source rows diagnosed',flush=True)
    summaries=[]
    for kind in kinds:
        group=[r for r in rows if r['kind']==kind]
        summaries.append(dict(kind=kind,rows=len(group),abs_difference_quantiles=np.quantile(abs(np.array([r['difference'] for r in group])),[.5,.9,.99,1]).tolist(),
            max_chieff_difference_from_native=max(abs(r['chieff_difference_from_native']) for r in group)))
    result=dict(schema='blackhole-spin-initialization-diagnosis/1',input_sha256=hashlib.sha256(source_path.read_bytes()).hexdigest(),
        purpose='Exploratory controlled variants after observing release discrepancy. No historical setting is inferred solely from a smaller residual.',
        summaries=summaries,rows=rows,failures=failures,warnings=warning_log)
    (BASE/'initialization-results.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps(summaries,indent=2),flush=True)


if __name__=='__main__':main()
