"""Read only whitelisted numeric/pandas pickle objects; export safe NPZ arrays."""
from pathlib import Path
import importlib
import json
import pickle
import numpy as np

ROOT = Path(__file__).resolve().parent
ALLOWED = {
 ('pandas.core.frame','DataFrame'),
 ('pandas.core.internals.managers','BlockManager'),
 ('pandas._libs.internals','_unpickle_block'),
 ('pandas.core.indexes.base','_new_Index'),
 ('pandas.core.indexes.base','Index'),
 ('numpy.core.numeric','_frombuffer'),
 ('numpy.core.multiarray','_reconstruct'),
 ('numpy.core.multiarray','scalar'),
 ('numpy','dtype'), ('numpy','ndarray'), ('builtins','slice'),
}

def numeric_array(reconstruct, args, state, extra=None):
    if reconstruct is not np._core.multiarray._reconstruct or args[0] is not np.ndarray:
        raise pickle.UnpicklingError('Not a numeric NumPy reconstruction')
    if extra is not None or np.dtype(state[2]).hasobject:
        raise pickle.UnpicklingError('Unexpected object array')
    array = reconstruct(*args)
    array.__setstate__(state)
    return array

class NumericUnpickler(pickle.Unpickler):
    def find_class(self, module, name):
        if (module,name) == ('dill._dill','_create_array'):
            return numeric_array
        if (module,name) not in ALLOWED:
            raise pickle.UnpicklingError(f'Forbidden global: {module}.{name}')
        return getattr(importlib.import_module(module.replace('numpy.core.', 'numpy._core.')), name)

def load(path):
    with path.open('rb') as f:
        return NumericUnpickler(f).load()

if __name__ == '__main__':
    events = load(ROOT/'data/gwtc-3-samples.pkl')
    inj = load(ROOT/'data/gwtc-3-injections.pkl')
    print('Event count:',len(events),'columns:',list(events[0].columns))
    print('Injection metadata:', {k:(np.shape(v) if np.ndim(v) else float(v)) for k,v in inj.items()})
    summaries=[]
    arrays={}
    for i,df in enumerate(events):
        d={k:df[k].to_numpy(dtype=float) for k in df.columns}
        d['chi_eff']=(d['a_1']*d['cos_tilt_1']+d['mass_ratio']*d['a_2']*d['cos_tilt_2'])/(1+d['mass_ratio'])
        d['mass_2']=d['mass_1']*d['mass_ratio']
        summary={'id':f'E{i+1:02d}','index':i,'samples':len(df)}
        for k in ['mass_1','mass_2','mass_ratio','chi_eff','redshift']:
            summary[k]=np.quantile(d[k],[.05,.5,.95]).tolist()
        summaries.append(summary)
        arrays.update({f'{i}_{k}':v for k,v in d.items()})
    np.savez_compressed(ROOT/'data/posteriors.npz',**arrays)
    np.savez_compressed(ROOT/'data/injections.npz',**inj)
    (ROOT/'data/event_summary.json').write_text(json.dumps(summaries,indent=2),encoding='utf-8')
    print('Summary:',[(s['id'],round(s['mass_1'][1],1),round(s['mass_2'][1],1),round(s['chi_eff'][1],2)) for s in summaries])
