"""Record reference-frame inputs without claiming a waveform-clock match."""
from pathlib import Path
import hashlib,json,sys,unittest
import h5py
from waveform_reference import interpolate_horizon_reference

root=Path(__file__).resolve().parent
out=root/'waveform-alignment-results'
out.mkdir(exist_ok=True)
rows=[]
for identifier,level in [('0389','Lev3'),('0389','Lev2'),('0001','Lev5'),('2327','Lev3')]:
    folder=root/'sxs-data'/identifier/level
    meta=json.loads((folder/'metadata.json').read_text())
    source_time=meta.get('reference_time',meta.get('relaxed_measurement_time'))
    with h5py.File(folder/'Horizons.h5') as f:
        result=interpolate_horizon_reference(meta,f,source_time)
    rows.append({'simulation':'SXS:BBH:'+identifier,'level':level,
        'purpose':'coordinate reference audit, not surrogate evaluation parameters',
        'metadata_sha256':hashlib.sha256((folder/'metadata.json').read_bytes()).hexdigest(),
        'horizons_sha256':hashlib.sha256((folder/'Horizons.h5').read_bytes()).hexdigest(),
        'parameters':result,'training_membership':'unverified',
        'radiation_frame_alignment':'not performed','independent_validation':False})
report={'date':'2026-10-23','records':rows,'implementation':{p:hashlib.sha256((root/p).read_bytes()).hexdigest() for p in ['waveform_reference.py','test_waveform_reference.py']}}
(out/'reference-audit.json').write_text(json.dumps(report,ensure_ascii=False,indent=2,allow_nan=False),encoding='utf-8')
suite=unittest.defaultTestLoader.discover(str(root),pattern='test_waveform_reference.py')
result=unittest.TextTestRunner(verbosity=2).run(suite)
(out/'tests.json').write_text(json.dumps({'tests_run':result.testsRun,'errors':len(result.errors),'failures':len(result.failures),'passed':result.wasSuccessful()},indent=2),encoding='utf-8')
sys.exit(0 if result.wasSuccessful() else 1)
