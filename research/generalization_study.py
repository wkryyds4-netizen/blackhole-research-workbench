"""Extend the existing conditional analysis without changing its data or model.

Evaluate Monte Carlo screening throughout f's posterior, not just at its median.
The result is an integration diagnostic, not a probability that a model is true.
"""
from pathlib import Path
import json,hashlib
import numpy as np

ROOT=Path(__file__).resolve().parents[1]

def analyze(rec,f):
    em=np.array(rec['event_moments']);n=len(em);counts=np.array(rec['sample_counts'])
    a,b,aa,bb,ab=rec['selection']
    e=f[:,None]*em[:,0]+(1-f[:,None])*em[:,1]
    e2=f[:,None]**2*em[:,2]+(1-f[:,None])**2*em[:,3]+2*f[:,None]*(1-f[:,None])*em[:,4]
    ess=counts*e**2/e2
    xi=f*a+(1-f)*b
    second=f*f*aa+(1-f)**2*bb+2*f*(1-f)*ab
    injess=xi**2/second
    var_events=((e2/e**2-1)/counts).sum(axis=1)
    mc=np.sqrt(np.maximum(0,var_events)+n*n*np.maximum(0,second/xi**2-1/rec['total_generated']))
    passed=(ess.min(axis=1)>=50)&(injess>=4*n)&(mc<=1)
    eventlog=np.log(e).sum(axis=1);logp=eventlog-n*np.log(xi);p=np.exp(logp-logp.max());p/=p.sum()
    score=-n*f*a/xi
    derivative=float(np.sum(p*(f-np.sum(p*f))*(score-np.sum(p*score))))
    def mean_at(k):
        lp=eventlog-n*np.log(f*a*k+(1-f)*b);pp=np.exp(lp-lp.max());pp/=pp.sum();return float(pp@f)
    eps=1e-5;numeric=(mean_at(np.exp(eps))-mean_at(np.exp(-eps)))/(2*eps)
    assert abs(derivative-numeric)<1e-7,(derivative,numeric)
    return dict(event_log_likelihood=(eventlog-eventlog.max()).round(9).tolist(),
                selection=rec['selection'],n_events=n,
                min_event_ess=ess.min(axis=1).round(3).tolist(),
                injection_ess=injess.round(3).tolist(),mc_log_sd=mc.round(6).tolist(),
                mc_event_var=var_events.round(9).tolist(),total_generated=rec['total_generated'],
                mc_pass=passed.tolist(),posterior_mass_passing=float(p[passed].sum()),
                mean_fraction=float(p@f),mean_sensitivity_log_efficiency=derivative,
                derivative_check_error=abs(derivative-numeric),median_interval=rec['corrected'])

def main():
    file=ROOT/'results/model.json';data=json.loads(file.read_text(encoding='utf-8'));f=np.array(data['f_grid'])
    scenarios={k:analyze(r,f) for k,r in data['scenarios'].items()}
    base=scenarios[data['default']]
    result=dict(schema='population-diagnostics/1',generated_date='2026-10-23',
        input_sha256=hashlib.sha256(file.read_bytes()).hexdigest(),default=data['default'],f_grid=f.tolist(),
        event_count=len(data['events']),scenarios=scenarios,
        interpretation='Conditional phenomenological tilt mixture; MC screening is diagnostic only. No averaging of the 81 models.',
        selection_audit=data['caveat'])
    for p in [ROOT/'results/generalization.json',ROOT/'simulator/public/research-data.json']:
        p.write_text(json.dumps(result,ensure_ascii=False,separators=(',',':'),allow_nan=False),encoding='utf-8')
    summary=dict(baseline_posterior_mass_passing=base['posterior_mass_passing'],
       baseline_mean_sensitivity_log_efficiency=base['mean_sensitivity_log_efficiency'],
       scenarios_passing_at_least_90_percent=sum(r['posterior_mass_passing']>=.9 for r in scenarios.values()),
       scenarios=81,max_derivative_check_error=max(r['derivative_check_error'] for r in scenarios.values()))
    (ROOT/'results/generalization-summary.json').write_text(json.dumps(summary,indent=2),encoding='utf-8')
    print(json.dumps(summary,indent=2))

if __name__=='__main__':main()
