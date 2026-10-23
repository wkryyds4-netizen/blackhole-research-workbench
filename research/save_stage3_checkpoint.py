"""Record local stage-three artifacts without including environments/build caches."""
from pathlib import Path
import hashlib,json
root=Path(__file__).resolve().parent.parent
paths={root/p for p in [
    'README.md','.gitignore','simulator/README.md',
    'research/stage3-plan.json','research/compare_nr_surrogate.py',
    'research/waveform_metrics.py','research/waveform_reference.py',
    'research/test_waveform_metrics.py','research/test_waveform_reference.py',
    'research/test_nr_comparison.py','research/download_nr_waveform.py',
    'research/load_nr_waveform.py','research/validate_nr_waveform.py',
    'research/waveform-metrics-method.md','research/waveform-metrics-validation.json',
    'research/save_stage3_checkpoint.py','simulator/public/workbench.html',
    'simulator/public/workbench.js','simulator/public/workbench.css',
    'simulator/public/waveform-view.mjs','simulator/public/projects.json',
    'simulator/public/datasets/waveform-comparison.json',
    'simulator/public/licenses/SXS-DATA.txt',
    'simulator/qa-waveform-stage.cjs','simulator/test-waveform-view.mjs',
]}
for directory in ['research/stage3-results','research/wave-source-results','research/waveform-alignment-results']:
    paths.update(p for p in (root/directory).iterdir() if p.is_file() and p.name!='CHECKPOINT.json')
report={'date':'2026-10-23','status':'evaluated-integrated-verified',
    'checks':{'synthetic_metrics':15,'reference_parameters_and_rotation':8,
        'independent_real_nr_checks':8,'new_frontend_contracts':12,'new_browser_checks':17,
        'existing_frontend_contracts':15,'existing_research_browser_checks':22,
        'source_reader':'0389 all3x77 modes exactly equal HDF5;2327 two levels rotation power invariant <=3.58e-15'},
    'limits':['External archive comparison; training membership unverified',
        'Complex h22 finite-record comparison, not detector-network or formation inference',
        '2327 decoded only; memory, normalization and radiation frame calibration pending'],
    'files':[{'path':p.relative_to(root).as_posix(),'bytes':p.stat().st_size,'sha256':hashlib.sha256(p.read_bytes()).hexdigest()} for p in sorted(paths)]}
out=root/'research/stage3-results/CHECKPOINT.json'
out.write_text(json.dumps(report,indent=2,ensure_ascii=False),encoding='utf-8')
print(f'Saved {len(report["files"])} file hashes to {out}')
