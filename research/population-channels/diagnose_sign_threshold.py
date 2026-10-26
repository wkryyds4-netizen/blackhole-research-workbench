"""Post-preview diagnostic of a discontinuous functional near spin zero.

Empirical dead zones do not estimate the true atom at zero or channel weights.
They expose sensitivity of a smoothed emulator to the strict-sign threshold.
"""
from pathlib import Path
import hashlib
import json
import numpy as np

BASE=Path(__file__).resolve().parent


def sha(path):return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    record=json.loads((BASE/'channel-preview.json').read_text())
    path=BASE/'channel-preview-samples.npz';assert sha(path)==record['sampling_cache_sha256']
    epsilons=[0.,1e-6,1e-4,1e-3,.01,.05]
    rows=[]
    with np.load(path,allow_pickle=False) as bank:
        for case in record['cases']:
            key=f'{case["channel"]}_birth{case["birth_spin"]:g}'
            data=bank[key];x=data[:,2]
            rows.append(dict(channel=case['channel'],birth_spin=case['birth_spin'],rows=len(x),
                             thresholds=[dict(epsilon=epsilon,negative=float((x < -epsilon).mean()),
                                              neutral=float((abs(x)<=epsilon).mean()),positive=float((x > epsilon).mean())) for epsilon in epsilons]))
    report=dict(schema='blackhole-channel-emulator-sign-threshold/1',
                timing='Post-preview exploratory source diagnostic; fixed existing draws, no GW fit.',
                input_sha256=sha(BASE/'channel-preview.json'),tool_sha256=sha(Path(__file__)),cases=rows,
                scope='Dead-zone sensitivity of learned continuous densities only. A neutral band does not identify a true zero-spin atom; strict negative fractions from a smoothed fixed-birth-spin model cannot be assumed to equal exact population-synthesis probabilities.')
    (BASE/'sign-threshold-results.json').write_text(json.dumps(report,indent=2)+'\n',encoding='utf-8')
    print(json.dumps([row for row in rows if row['birth_spin']==0],indent=2))


if __name__=='__main__':main()
