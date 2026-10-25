"""Independent indicator sums and direct products for the search-policy control."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.integrate import cumulative_trapezoid
from run_selection_study import init_model,physical,weights

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]


def main():
    result=json.loads((BASE/'pipeline-results.json').read_text())
    source=json.loads((BASE/'o3-provenance.json').read_text())
    study=json.loads((BASE/'selection-results.json').read_text())
    assert hashlib.sha256((BASE/'o3-bank.npz').read_bytes()).hexdigest()==result['inputs']['o3_bank']
    for path,sha in study['frozen_hashes'].items():assert hashlib.sha256((ROOT/path).read_bytes()).hexdigest()==sha
    init_model()
    with np.load(BASE/'o3-bank.npz',allow_pickle=False) as bank:raw={k[5:]:bank[k] for k in bank.files if k.startswith('full_')}
    w=weights(physical(raw))
    # Explicit, documented pipeline set; no substring-based union construction.
    found=np.maximum.reduce([raw['ifar_gstlal'],raw['ifar_mbta'],raw['ifar_pycbc_bbh'],raw['ifar_pycbc_hyperbank']])>1
    assert int(found.sum())==result['counts']['full_matched_filter']
    factor=source['attributes']['analysis_time_s']/source['combined_attributes']['analysis_time_s']/source['attributes']['total_generated']
    reference=np.array(study['selection_components']['o1'])+np.array(study['selection_components']['o2'])
    direct=reference+factor*np.einsum('n,nki->ki',found.astype(float),w)
    error=float(np.max(abs(direct/np.array(result['selection_components']['matched_filter'])-1)))
    assert error<1e-10
    errors=[]
    with np.load(BASE/'pipeline-checkpoints.npz',allow_pickle=False) as bank:
        r=bank['r'];fs=bank['fs'];means=bank['means']
        for name in ['union','matched_filter']:
            for k in range(6):
                x=bank['selection_'+name][k]
                integrals=1+fs[k,None,:]*(means[:,k,0,None]/means[:,k,1,None]-1)
                efficiency=1+fs[k]*(x[0]/x[1]-1)
                product=np.prod(integrals.astype(np.longdouble),axis=0)/efficiency.astype(np.longdouble)**len(means)
                cdf=cumulative_trapezoid(np.asarray(product/product.max(),float),r,initial=0);cdf/=cdf[-1]
                q=np.interp([.05,.5,.95],cdf,r)
                errors.append(float(np.max(abs(q-np.array(result[name+'_quantiles'])[k]))))
    assert max(errors)<1e-11
    report=dict(schema='blackhole-pipeline-policy-validation/1',status='pass',
        all_accepted_indicator_relative_error=error,full_product_quantile_max_error=max(errors),
        frozen_hashes_verified=study['frozen_hashes'],
        scope='Estimator and pipeline-indicator validation only; event-search compatibility is documentary, not an independent FAR-log audit.')
    (BASE/'pipeline-validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print(json.dumps(report,indent=2))


if __name__=='__main__':main()
