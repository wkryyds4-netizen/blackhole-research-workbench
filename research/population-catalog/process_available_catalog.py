"""Run the dependent catalog stages with fail-fast validation.

An optional finite follow mode processes newly verified downloads in batches.
It stops after the fixed 69-file plan, a download failure, any scientific check
failure, or a bounded timeout. It does not download, publish, or schedule work.
"""
from pathlib import Path
import argparse
from datetime import datetime,timezone
import hashlib
import json
import os
import subprocess
import sys
import time

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]


def atomic(path,value):
    temporary=path.with_suffix('.tmp')
    temporary.write_text(json.dumps(value,indent=2),encoding='utf-8')
    temporary.replace(path)


def main():
    parser=argparse.ArgumentParser()
    parser.add_argument('--follow-downloads',action='store_true')
    parser.add_argument('--maximum-minutes',type=float,default=90)
    args=parser.parse_args();assert 0<args.maximum_minutes<=120
    lock=BASE/'processing.lock'
    descriptor=os.open(lock,os.O_CREAT|os.O_EXCL|os.O_WRONLY)
    os.write(descriptor,str(os.getpid()).encode());os.close(descriptor)
    started=time.time();last=set();status_path=BASE/'processing-status.json'
    adapter=ROOT/'simulator/.adapters-venv/Scripts/python.exe'
    assert adapter.exists()
    state=dict(schema='blackhole-native-catalog-processing/1',started_at=datetime.now(timezone.utc).isoformat(),
        follow_downloads=args.follow_downloads,planned_source_files=69,cycles=[])
    def save(phase,**extra):
        state.update(phase=phase,updated_at=datetime.now(timezone.utc).isoformat(),**extra);atomic(status_path,state)
    try:
        while True:
            registry=json.loads((BASE/'download-results.json').read_text())
            available={r['index'] for r in registry if r['whole_file_md5_verified']}
            batch_path=BASE/'download-batch-status.json'
            terminal=None
            if batch_path.exists() and batch_path.stat().st_mtime>=started:
                terminal=json.loads(batch_path.read_text())
            if not last or len(available-last)>=4 or len(available)==69 or terminal:
                if available!=last:
                    save('processing',verified_downloads=len(available),active_step=None)
                    scripts=[(str(adapter),'extract_verified_events.py',['--available'])]
                    scripts += [(sys.executable,name,[]) for name in [
                        'audit_configuration_rows.py','run_native_catalog.py','validate_native_catalog.py',
                        'summarize_catalog.py','plot_native_catalog.py','plot_catalog.py']]
                    for interpreter,name,options in scripts:
                        save('processing',active_step=name)
                        command=[interpreter,str(BASE/name),*options]
                        with (BASE/'processing-step.log').open('w',encoding='utf-8') as output:
                            completed=subprocess.run(command,cwd=ROOT,stdout=output,stderr=subprocess.STDOUT)
                        if completed.returncode:
                            tail=(BASE/'processing-step.log').read_text(encoding='utf-8',errors='replace')[-5000:]
                            raise RuntimeError(f'{name} exited {completed.returncode}:\n{tail}')
                    result_path=BASE/'native-catalog-results.json';result=json.loads(result_path.read_text())
                    validation=json.loads((BASE/'native-catalog-validation.json').read_text())
                    assert validation['result_sha256']==hashlib.sha256(result_path.read_bytes()).hexdigest()
                    cycle=dict(verified_downloads=len(available),full_sources_in_likelihood=result['completed_full_source_events'],
                        width_gaps_pp=result['width_gaps_pp'],median_shift_pp=result['median_shift_pp'],
                        independent_quantile_error=validation['full_product_quantile_error'],result_sha256=validation['result_sha256'])
                    snapshots=BASE/'native-snapshots';snapshots.mkdir(exist_ok=True)
                    count=cycle['full_sources_in_likelihood']
                    for name in ['native-catalog-results.json','native-catalog-validation.json']:
                        target=snapshots/f'{count}-{name}';content=(BASE/name).read_bytes()
                        if target.exists():assert target.read_bytes()==content,'Snapshot collision: preserve and inspect before replacing.'
                        else:target.write_bytes(content)
                    state['cycles'].append(cycle);print(json.dumps(cycle),flush=True)
                    last=available
                if len(available)==69:
                    result=json.loads((BASE/'native-catalog-results.json').read_text())
                    assert result['completed_full_source_events']==69,'All files available, but source/prior resolution is incomplete.'
                    save('complete',active_step=None);break
                if terminal and terminal['missing']:
                    raise RuntimeError(f'Download batch ended with missing source indices: {terminal["missing"]}')
                if not args.follow_downloads:
                    save('available_batch_complete',active_step=None);break
            if time.time()-started>args.maximum_minutes*60:
                save('bounded_timeout',active_step=None);raise TimeoutError('Follow duration ended; downloads have not been stopped. Resume from saved state.')
            save('waiting_for_verified_sources',active_step=None)
            time.sleep(15)
    except Exception as error:
        save('failed',error=type(error).__name__+': '+str(error));raise
    finally:
        lock.unlink()


if __name__=='__main__':main()
