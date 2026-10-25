"""Check ranked rows directly against the checksum-verified official O3 HDF."""
from pathlib import Path
import json
import hashlib
import h5py

BASE=Path(__file__).resolve().parent


def main():
    result=json.loads((BASE/'policy-influence-results.json').read_text())
    source=json.loads((BASE/'o3-provenance.json').read_text())['source']
    path=BASE/'sources/endo3_bbhpop-LIGO-T2100113-v12.hdf5'
    assert hashlib.sha256(path.read_bytes()).hexdigest()==source['sha256']
    checks=[]
    with h5py.File(path) as h:
        g=h['injections']
        for family in result['families']:
            for row in family['ranked_rows']:
                i=row['source_row'];expected=dict(row['parameters'],**row['ifar'])
                assert all(float(g[k][i])==v for k,v in expected.items())
                assert float(g['ifar_cwb'][i])>1 and max(float(g[k][i]) for k in row['ifar'] if k!='ifar_cwb')<=1
                checks.append(dict(family=family['label'],source_row=i,all_stored_values_exact=True,cwb_only=True))
    report=dict(status='pass',source_sha256=source['sha256'],checked_rows=checks,
        source_row_convention='zero based within full O3 injections datasets')
    (BASE/'policy-influence-source-validation.json').write_text(json.dumps(report,indent=2),encoding='utf-8')
    print('Verified all ranked-row records directly against official HDF.')


if __name__=='__main__':main()
