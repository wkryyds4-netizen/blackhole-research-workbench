"""Independent full products and accepted-bank selection normalization checks."""
from pathlib import Path
import hashlib
import json
import sys
import numpy as np
from scipy.integrate import cumulative_trapezoid

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]
sys.path.insert(0,str(BASE))
from run_selection_study import physical,init_model,weights


def main():
    result=json.loads((BASE/'selection-results.json').read_text());source=json.loads((BASE/'o3-provenance.json').read_text())
    for path,sha in result['frozen_hashes'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha
    assert hashlib.sha256((BASE/'o3-bank.npz').read_bytes()).hexdigest()==result['o3_bank_sha256']
    errors=[]
    with np.load(BASE/'selection-checkpoints.npz',allow_pickle=False) as cp:
        r=cp['r'];fs=cp['fs']
        for si,name in enumerate(cp['names'].tolist()):
            for k in range(6):
                means=cp['means'][si,:,k,:];xi=cp['selection'][si,k]
                # Ratios remove large hyperparameter-independent constants before multiplication.
                integrals=1+fs[k,None,:]*(means[:,0,None]/means[:,1,None]-1)
                efficiency=1+fs[k]*(xi[0]/xi[1]-1)
                p=np.prod(integrals.astype(np.longdouble),axis=0)/efficiency.astype(np.longdouble)**len(means)
                c=cumulative_trapezoid(np.asarray(p/p.max(),dtype=float),r,initial=0);c/=c[-1]
                q=np.interp([.05,.5,.95],c,r)
                errors.append(float(np.max(abs(q-result['scenarios'][name]['quantiles'][k]))))
    assert max(errors)<1e-11
    init_model()
    with np.load(BASE/'o3-bank.npz',allow_pickle=False) as cp:
        raw={k[5:]:cp[k] for k in cp.files if k.startswith('full_')}
        mapping=cp['combined_o3_to_full']
    # Sum on the entire accepted bank with an explicit indicator, including failed detections.
    w=weights(physical(raw));found=np.maximum.reduce([raw[k] for k in raw if k.startswith('ifar_')])>1
    accepted_sum=np.einsum('n,nki->ki',found.astype(float),w)
    attrs=source['attributes'];old=source['combined_attributes']
    c=attrs['analysis_time_s']/old['analysis_time_s']/attrs['total_generated']
    target=np.array(result['selection_components']['full_o3'])-np.array(result['selection_components']['o1'])-np.array(result['selection_components']['o2'])
    full_error=float(np.max(abs(accepted_sum*c/target-1)));assert full_error<1e-10
    fraction=len(mapping)/len(w);mw=source['mixture_weight_combined_o3_range'][0]
    identity_error=abs(mw*fraction/old['total_generated']/c-1);assert identity_error<1e-12
    subset=np.einsum('n,nki->ki',found[mapping].astype(float),w[mapping])
    recovered=subset*mw/old['total_generated']
    subset_error=float(np.max(abs(recovered/np.array(result['selection_components']['o3'])-1)));assert subset_error<1e-10
    report=dict(schema='blackhole-full-o3-validation/1',status='pass',
        full_product_max_quantile_error=max(errors),all_accepted_selection_sum_max_relative_error=full_error,
        recovered_combined_o3_max_relative_error=subset_error,
        time_trial_mixture_thinning_identity_relative_error=identity_error,
        redshift_cap_normalization_quantile_error=result['normalized_cap_quantile_error'],
        grid_max_quantile_error=result['max_grid_quantile_error'],frozen_hashes_verified=result['frozen_hashes'],
        interpretation='Validates the declared finite-sample estimator and exact source subset; does not establish actual uniform random thinning, repeated-campaign coverage, or absent high-redshift sensitivity.')
    (BASE/'selection-validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
