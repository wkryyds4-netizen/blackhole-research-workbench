"""Reproducible conditional two-component population inference, GWTC-3.

The components are phenomenological spin-tilt populations, NOT identified
astrophysical formation channels. Full tutorial posterior draws and injection
draws are importance reweighted. N is conditioned on; no rate inference.
"""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.special import betaln, ndtr, logsumexp
from scipy.optimize import linear_sum_assignment
from astropy.cosmology import Planck15

ROOT = Path(__file__).resolve().parent
SIGMAS = [.3, .6, 1.2]
SPINS = [(1.5,8.5),(2,5),(1,1)]
MASSES = [(2,.1),(3,.1),(2,0)]
LAMBDAS = [0,2,4]
F = (np.arange(501)+.5)/501

def power_pdf(x, slope, low, high):
    norm=(high**(1-slope)-low**(1-slope))/(1-slope)
    return np.where((x>=low)&(x<=high),x**(-slope)/norm,0)

def magnitude_pdf(x,a,b):
    v=np.clip(x,1e-12,1-1e-12)
    return np.exp((a-1)*np.log(v)+(b-1)*np.log1p(-v)-betaln(a,b))

def aligned_tilt(x,sigma):
    norm=sigma*np.sqrt(2*np.pi)*(ndtr(0)-ndtr(-2/sigma))
    return np.exp(-.5*((x-1)/sigma)**2)/norm

def identify_events():
    summaries=json.loads((ROOT/'data/event_summary.json').read_text())
    catalog={e['commonName']:e for p in sorted((ROOT/'data').glob('GWTC*.json'))
             for e in json.loads(p.read_text())['events'].values() if e['mass_1_source'] is not None}
    events=list(catalog.values())
    def cost(s,e):
        return sum(((s[k][1]-e[c])/max((s[k][2]-s[k][0])/2,floor))**2
                   for k,c,floor in [('mass_1','mass_1_source',1),('mass_2','mass_2_source',1),('chi_eff','chi_eff',.05),('redshift','redshift',.05)])
    scores=np.array([[cost(s,e) for e in events] for s in summaries])
    rows,cols=linear_sum_assignment(scores)
    out=[]
    for r,c in zip(rows,cols):
        s=summaries[r]; e=events[c]
        s.update(name=e['commonName'],match_cost=float(scores[r,c]),far=e['far'],gps=e['GPS'],catalog=e['catalog.shortName'])
        assert s['match_cost']<1, 'Ambiguous catalog match: manual audit needed'
        s['included']=s['mass_2'][1]>3 and (e['far'] is None or e['far']<1)
        s['exclusion']=None if s['included'] else ('secondary median < 3 solar masses' if s['mass_2'][1]<3 else 'catalog FAR >= 1/year')
        out.append(s)
    # Independently verify the matching retains the chronology of all 73 inputs.
    assert np.all(np.diff([s['gps'] for s in out])>0)
    assert sum(s['included'] for s in out)==69
    (ROOT/'data/sample_audit.json').write_text(json.dumps(out,ensure_ascii=False,indent=2),encoding='utf-8')
    return out

def base_density(d,spin,mass,lam):
    m=d['mass_1']; q=d['mass_ratio']; z=d['redshift']
    alpha,peak=MASSES[mass]
    mmin,mmax=3.,100.
    gaussian=np.exp(-.5*((m-35)/5)**2)/(5*np.sqrt(2*np.pi)*(ndtr((mmax-35)/5)-ndtr((mmin-35)/5)))
    pm=((1-peak)*power_pdf(m,alpha,mmin,mmax)+peak*gaussian)*((m>=mmin)&(m<=mmax))
    qmin=mmin/m
    pq=np.where((q>=qmin)&(q<=1)&(qmin<1),2*q/np.maximum(1-qmin**2,1e-12),0)
    a,b=SPINS[spin]
    pz=np.interp(z,ZGRID,ZDV)*(1+z)**(lam-1)/ZNORM[lam]
    pz*=((z>=0)&(z<=2.3))
    return pm*pq*magnitude_pdf(d['a_1'],a,b)*magnitude_pdf(d['a_2'],a,b)*pz

def component_weights(d,spin,mass,lam,sigma):
    base=base_density(d,spin,mass,lam)/d['prior']
    return base*aligned_tilt(d['cos_tilt_1'],sigma)*aligned_tilt(d['cos_tilt_2'],sigma),base*.25

def reduce_weights(a,b):
    return np.array([a.mean(),b.mean(),np.mean(a*a),np.mean(b*b),np.mean(a*b)])

def density(rec,correct=True,prior='uniform',drop=None,calibration=1):
    moments=np.array(rec['event_moments'])
    if drop is not None:
        moments=np.delete(moments,drop,axis=0)
    mix=F[:,None]*moments[:,0]+(1-F[:,None])*moments[:,1]
    logs=np.log(mix).sum(axis=1)
    if correct:
        a,b=rec['selection'][:2]
        logs-=len(moments)*np.log(F*a*calibration+(1-F)*b)
    if prior=='central':
        logs+=np.log(F*(1-F))
    elif prior=='edge':
        logs-=.5*np.log(F*(1-F))
    p=np.exp(logs-logsumexp(logs))*len(F)
    return p

def interval(p):
    c=np.cumsum(p)/np.sum(p)
    return np.interp([.05,.5,.95],c,F).tolist()

def diagnostics(rec,f):
    em=np.array(rec['event_moments']); n=len(em)
    e=f*em[:,0]+(1-f)*em[:,1]
    e2=f*f*em[:,2]+(1-f)**2*em[:,3]+2*f*(1-f)*em[:,4]
    counts=np.array(rec['sample_counts'])
    ess=counts*e*e/e2
    a,b,aa,bb,ab=rec['selection']
    xi=f*a+(1-f)*b
    second=f*f*aa+(1-f)**2*bb+2*f*(1-f)*ab
    injess=xi*xi/second
    var_events=np.sum((e2/e**2-1)/counts)
    var_sel=n*n*max(0,second/xi**2-1/rec['total_generated'])
    return dict(min_event_ess=float(ess.min()),injection_ess=float(injess),log_likelihood_mc_sd=float(np.sqrt(var_events+var_sel)))

def main():
    global ZGRID,ZDV,ZNORM
    ZGRID=np.linspace(0,2.3,3001)
    ZDV=Planck15.differential_comoving_volume(ZGRID).value
    ZNORM={lam:np.trapezoid(ZDV*(1+ZGRID)**(lam-1),ZGRID) for lam in LAMBDAS}
    audit=identify_events(); selected=[s for s in audit if s['included']]
    npz=np.load(ROOT/'data/posteriors.npz')
    keys=['mass_1','mass_ratio','redshift','a_1','a_2','cos_tilt_1','cos_tilt_2','prior']
    events=[{k:npz[f'{s["index"]}_{k}'] for k in keys} for s in selected]
    injections=dict(np.load(ROOT/'data/injections.npz'))
    total=float(injections['total_generated'])
    scenarios={}
    for si,sigma in enumerate(SIGMAS):
        for ai in range(3):
            for mi in range(3):
                for zi,lam in enumerate(LAMBDAS):
                    key=f'{si}-{ai}-{mi}-{zi}'
                    weights=[component_weights(d,ai,mi,lam,sigma) for d in events]
                    event_moments=[reduce_weights(a,b).tolist() for a,b in weights]
                    a,b=component_weights(injections,ai,mi,lam,sigma)
                    selection=[a.sum()/total,b.sum()/total,(a*a).sum()/total**2,(b*b).sum()/total**2,(a*b).sum()/total**2]
                    rec=dict(event_moments=event_moments,sample_counts=[len(d['prior']) for d in events],selection=selection,total_generated=total)
                    rec['corrected']=interval(density(rec))
                    rec['uncorrected']=interval(density(rec,False))
                    rec['diagnostics']=diagnostics(rec,rec['corrected'][1])
                    scenarios[key]=rec
        print(f'Completed tilt width {sigma}',flush=True)
    data=dict(events=selected,scenarios=scenarios,default='1-1-0-1',f_grid=F.tolist(),
              n_found=len(injections['prior']),n_generated=int(total),analysis_time_years=float(injections['analysis_time']),
              sources=json.loads((ROOT/'data/manifest.json').read_text()),
              tutorial='https://colmtalbot.github.io/gwpopulation/examples/GWTC3.html',
              caveat='Conditional phenomenological tilt-mixture analysis; preprocessed tutorial inputs, not an LVK reproduction. The found-injection threshold and waveform identities are not retained in the tutorial files.')
    (ROOT/'results').mkdir(exist_ok=True)
    (ROOT/'results/model.json').write_text(json.dumps(data,separators=(',',':'),ensure_ascii=False),encoding='utf-8')
    default=scenarios[data['default']]
    loo=[]
    for i,e in enumerate(selected):
        ci=interval(density(default,drop=i)); loo.append(dict(name=e['name'],interval=ci,delta=ci[1]-default['corrected'][1]))
    report=dict(baseline=default,scenario_median_range=[min(r['corrected'][1] for r in scenarios.values()),max(r['corrected'][1] for r in scenarios.values())],leave_one_out=sorted(loo,key=lambda d:abs(d['delta']),reverse=True),
                worst_diagnostics={k:max(r['diagnostics'][k] for r in scenarios.values()) if k=='log_likelihood_mc_sd' else min(r['diagnostics'][k] for r in scenarios.values()) for k in default['diagnostics']})
    (ROOT/'results/summary.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in report.items() if k not in ['baseline','leave_one_out']},indent=2))
    print('Baseline:',default['corrected'],default['uncorrected'],default['diagnostics'])
    print('Top LOO:',report['leave_one_out'][:3])

if __name__=='__main__': main()
