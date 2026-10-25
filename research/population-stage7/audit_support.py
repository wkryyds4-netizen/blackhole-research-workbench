"""Exploratory retained-row boundary diagnostics, explicitly not missing-tail bounds."""
from pathlib import Path
import sys
import json
import argparse
import numpy as np
from scipy.special import betainc,betaincc
from scipy.integrate import quad
from scipy.stats import beta

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
sys.path.insert(0,str(ROOT/'research/population-stage5'))
from run_tail_study import density_pair,init_model,weights,coordinates,KEYS,PARAMS


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--indices',nargs='+',type=int,default=[19,20,25,50]);args=parser.parse_args()
    init_model();provenance=json.loads((BASE/'event-provenance.json').read_text())
    reference=np.array(json.loads((ROOT/'research/population-stage6/catalog-screen.json').read_text())['baseline_quantiles'])
    r,fs=coordinates();out=[];tail_checks=[]
    for ix in args.indices:
        rec=next(x for x in provenance if x['index']==ix);bounds=rec['prior_audit']['bounds']
        with np.load(BASE/f'event-{ix}.npz',allow_pickle=False) as bank:
            d={k:bank[f'nocosmo_{k}'] for k in PARAMS};_,pn,_=density_pair(d)
            w=weights(dict(d,prior=pn));fractions=[]
            for k in range(6):
                f=float(np.interp(reference[k,1],r,fs[k]));raw=f*w[:,k,0]+(1-f)*w[:,k,1]
                fractions.append(raw/raw.sum())
            fields={}
            for name,bound in bounds.items():
                if name.startswith('tilt_'):values=np.arccos(d['cos_'+name])
                else:values=bank['nocosmo_'+('mass_1_detector' if name=='mass_1' else name)]
                lo,hi=bound['minimum'],bound['maximum'];width=hi-lo
                masks={'lower':values<=lo+.01*width,'upper':values>=hi-.01*width}
                fields[name]=dict(bounds=bound,observed_min=float(values.min()),observed_max=float(values.max()),
                    bands={side:dict(rows=int(mask.sum()),row_fraction=float(mask.mean()),
                                    integral_fraction_at_original_medians=[float(v[mask].sum()) for v in fractions])
                           for side,mask in masks.items()})
            spin_tails=[]
            for a,b in [(1.5,8.5),(2,5)]:
                per=[]
                for name in ['a_1','a_2']:
                    lo,hi=bounds[name]['minimum'],bounds[name]['maximum']
                    tail=float(betainc(a,b,lo)+betaincc(a,b,hi))
                    numerical=quad(lambda x:beta.pdf(x,a,b),0,lo,epsabs=1e-24)[0]+quad(lambda x:beta.pdf(x,a,b),hi,1,epsabs=1e-24)[0]
                    assert abs(numerical/tail-1)<1e-8
                    tail_checks.append(abs(numerical/tail-1));per.append(tail)
                spin_tails.append(dict(beta=[a,b],population_probability_outside_either_spin_support=per[0]+per[1]-per[0]*per[1]))
            out.append(dict(index=ix,event=rec['event'],n_rows=len(pn),fields=fields,spin_support_population_mass=spin_tails))
    result=dict(schema='blackhole-retained-support-diagnostic/1',boundary_band_fraction=.01,keys=KEYS,events=out,
        beta_tail_quadrature_max_relative_error=max(tail_checks),
        design='Exploratory follow-up defined after two main comparisons; see design.md.',
        caution='Observed-row diagnostics cannot bound unobserved likelihood tails. The Beta probability is population mass, not posterior mass or an event-evidence error bound.')
    suffix='' if args.indices==[19,20,25,50] else 'partial-'+'-'.join(map(str,args.indices))+'-'
    (BASE/(suffix+'support-audit.json')).write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps([{k:e[k] for k in ['event','n_rows','spin_support_population_mass']} for e in out],indent=2))


if __name__=='__main__':main()
