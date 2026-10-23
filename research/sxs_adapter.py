"""Reproducible SpEC Horizons.h5 -> blackhole-scene/1 adapter.

No SXS package code is vendored. Reads the documented time/value datasets
with h5py. Coordinates are gauge dependent; these are not observed paths.
"""
from pathlib import Path
import hashlib, json
import h5py
import numpy as np

ROOT = Path(__file__).resolve().parent
PUBLIC = ROOT.parent/'simulator'/'public'/'datasets'
FIELDS = {'position':'CoordCenterInertial.dat', 'mass':'ChristodoulouMass.dat', 'spin':'chiInertial.dat'}
ROTATION = np.array([[1,0,0],[0,0,-1],[0,1,0]])  # proper Rx(+pi/2), det=+1
TOL = np.array([.003,.003,.003,1e-5,1e-4,1e-4,1e-4])
DATASETS = [('2327','Lev3'),('2327','Lev2'),('2327','Lev1'),('0001','Lev5'),('0003','Lev5'),('0389','Lev3'),('0389','Lev2')]

def write(path, data):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(data, ensure_ascii=False, separators=(',', ':'), allow_nan=False), encoding='utf-8')

def read(directory):
    meta = json.loads((directory/'metadata.json').read_text())
    manifest = json.loads((directory/'manifest.json').read_text())
    for item in manifest['files']:
        assert hashlib.sha256((directory/item['file']).read_bytes()).hexdigest() == item['sha256']
    prefix = 'reference' if 'reference_mass1' in meta else 'relaxed'
    mass = meta[prefix+'_mass1'] + meta[prefix+'_mass2']
    start = meta.get('reference_time', meta.get('relaxed_measurement_time'))
    assert mass > 0 and start is not None
    groups = {}
    with h5py.File(directory/'Horizons.h5') as file:
        for name in 'ABC':
            groups[name] = {}
            for field, dataset in FIELDS.items():
                a = np.asarray(file[f'Ah{name}.dir/{dataset}'])
                assert np.isfinite(a).all() and np.all(np.diff(a[:,0]) > 0)
                assert a.shape[1] == (2 if field == 'mass' else 4)
                groups[name][field] = a
    return meta, manifest, groups, mass, start

def bounds(group):
    return max(a[0,0] for a in group.values()), min(a[-1,0] for a in group.values())

def evaluate(group, times, mass):
    values = {}
    for field, a in group.items():
        assert times[0] >= a[0,0] and times[-1] <= a[-1,0], 'Extrapolation forbidden'
        v = np.column_stack([np.interp(times, a[:,0], a[:,j]) for j in range(1,a.shape[1])])
        values[field] = v@ROTATION.T if field in ('position','spin') else v
        if field != 'spin': values[field] /= mass
    return np.column_stack([values['position'], values['mass'], values['spin']])

def compress(t, y, tol):
    """Piecewise-linear simplification with a componentwise native-grid bound."""
    keep = {0, len(t)-1}; stack = [(0,len(t)-1)]
    while stack:
        a,b = stack.pop()
        if b-a < 2: continue
        u = (t[a+1:b]-t[a])/(t[b]-t[a])
        prediction = y[a] + u[:,None]*(y[b]-y[a])
        errors = np.max(np.abs(y[a+1:b]-prediction)/tol, axis=1)
        k = int(np.argmax(errors))+a+1
        if errors[k-a-1] > 1:
            keep.add(k); stack.extend([(a,k),(k,b)])
    return np.array(sorted(keep))

def adapt(directory):
    meta, manifest, groups, mass, reference_time = read(directory)
    common = bounds(groups['C'])[0]
    assert abs(common-meta['common_horizon_time']) < 1e-5
    binary_end = min(*(bounds(groups[k])[1] for k in 'AB'), np.nextafter(common,-np.inf))
    gap = common-binary_end
    start = max(reference_time, *(bounds(groups[k])[0] for k in 'AB'))
    frames, stats = [], []
    for ids, lo, hi in [('AB',start,binary_end), ('C',common,bounds(groups['C'])[1])]:
        # Include every source knot so error bounds hold between knots as well.
        t = np.unique(np.concatenate([[lo,hi]]+[a[:,0][(a[:,0]>=lo)&(a[:,0]<=hi)] for k in ids for a in groups[k].values()]))
        y = np.column_stack([evaluate(groups[k],t,mass) for k in ids])
        selected = compress(t,y,np.tile(TOL,len(ids)))
        reconstructed = np.column_stack([np.interp(t,t[selected],y[selected,j]) for j in range(y.shape[1])])
        error = np.abs(y-reconstructed).reshape(len(t),len(ids),7)
        stats.append({'ids':ids, 'native_union_samples':len(t), 'display_samples':len(selected),
                      'max_position_component_error_M':float(error[:,:,:3].max()),
                      'max_mass_error_fraction':float(error[:,:,3].max()),
                      'max_spin_component_error':float(error[:,:,4:].max())})
        # The AB endpoint at common is used to bound error but is not displayed:
        # IDs switch discretely to C, never a morph or extrapolated trajectory.
        if ids == 'AB': selected = selected[t[selected] < common]
        for i in selected:
            bodies = []
            for j,k in enumerate(ids):
                v = y[i,j*7:(j+1)*7]
                assert v[3] > 0 and np.linalg.norm(v[4:]) <= 1
                bodies.append({'id':k,'position':v[:3].tolist(),'mass':float(v[3]),'spin':v[4:].tolist()})
            frames.append({'t':float((t[i]-common)/mass),'bodies':bodies})
    assert 0 < len(frames) <= 30000
    assert all(b['t'] > a['t'] for a,b in zip(frames,frames[1:]))
    id_ = manifest['simulation'].split(':')[-1]
    prefix='reference' if 'reference_mass1' in meta else 'relaxed'
    omega=np.asarray(meta[prefix+'_orbital_frequency']); axis=omega/np.linalg.norm(omega)
    transverse=max(np.linalg.norm(np.cross(np.asarray(meta[prefix+'_dimensionless_spin'+str(i)]),axis)) for i in [1,2])
    classification='进动双黑洞' if transverse>1e-3 else '近无自旋双黑洞' if max(np.linalg.norm(meta[prefix+'_dimensionless_spin'+str(i)]) for i in [1,2])<1e-3 else '近对齐自旋双黑洞'
    deprecated='deprecated' in meta.get('keywords',[])
    title = f"SXS:BBH:{id_} · {manifest['level']} · {classification}"+('（旧版）' if deprecated else '')
    repository=manifest.get('repository','Zenodo')
    license_id=manifest.get('license',{}).get('id','unspecified')
    provenance = {'kind':'numerical','source':f"SXS Collaboration / SpEC Horizons; https://doi.org/{manifest['doi']}",
        'version':f"{repository} record {manifest['record_id']}; {manifest['level']}; sxs-adapter/2",
        'frame':'SpEC inertial coordinate centers; proper rotation (x,y,z) → (x,−z,y); no COM subtraction',
        'license':'CC-BY-4.0' if license_id=='cc-by-4.0' else '记录未声明许可；仅本地研究快照',
        'license_record':manifest.get('license'), 'credit':'SXS Collaboration; creators and citation keys retained in downloaded record/metadata',
        'clock':'(source coordinate time − first common-horizon sample) / reference total mass; not waveform retarded time',
        'reference_total_mass_in_source_units':mass, 'source_common_horizon_time':common,
        'source_reference_time':reference_time, 'rotation_matrix':ROTATION.tolist(),
        'physical_scale':'65 solar masses is an arbitrary display scale, not a measured event mass; no redshift applied',
        'sampling':'linear per quantity on its own native time grid; no extrapolation; adaptive display simplification',
        'transition':'A/B replaced by C at first available common-horizon sample; common-horizon birth is not waveform amplitude peak',
        'gaps':[{'start':float((binary_end-common)/mass),'end':0,'behavior':'hold last A/B sample; no interpolation across topology'}] if gap>1e-7 else [],
        'archive_note':'Pinned data record; deprecated status from source metadata; no automatic latest-version substitution',
        'deprecated':deprecated, 'superseded_by':meta.get('superseded_by'),
        'pbj_info':meta.get('pbj_info'), 'source_reference_orbital_frequency':meta[prefix+'_orbital_frequency'],
        'source_files':manifest['files'],'citation_keys':meta.get('simulation_bibtex_keys',[])}
    scene = {'schema':'blackhole-scene/1','title':title,
        'units':{'time':'GM/c^3','length':'GM/c^2','reference_mass_msun':65,'mass_frame':'geometric'},
        'provenance':provenance,'frames':frames}
    filename = f'sxs-{id_}-{manifest["level"].lower()}.json'
    write(PUBLIC/filename,scene)
    assert (PUBLIC/filename).stat().st_size <= 5_000_000
    report = {'simulation':manifest['simulation'],'level':manifest['level'],'frames':len(frames),
        'reference_total_mass':mass, 'reference_mass_ratio':meta.get('reference_mass1',meta.get('relaxed_mass1'))/meta.get('reference_mass2',meta.get('relaxed_mass2')),
        'reference_spin1':meta.get('reference_dimensionless_spin1',meta.get('relaxed_dimensionless_spin1')),
        'reference_spin2':meta.get('reference_dimensionless_spin2',meta.get('relaxed_dimensionless_spin2')),
        'common_horizon_display_time':0,'pre_common_data_gap_M':float(gap/mass),'last_mass_fraction':frames[-1]['bodies'][0]['mass'],
        'last_spin_magnitude':float(np.linalg.norm(frames[-1]['bodies'][0]['spin'])), 'compression':stats,
        'scene_sha256':hashlib.sha256((PUBLIC/filename).read_bytes()).hexdigest()}
    return {'id':f'sxs-{id_}-{manifest["level"].lower()}','title':title,'url':'datasets/'+filename},report

def resolution_comparison():
    # Same simulation's source coordinate time, same reference scale, no fitted
    # phase/time alignment. This is a resolution difference, not convergence order.
    low=read(ROOT/'sxs-data/0389/Lev2'); high=read(ROOT/'sxs-data/0389/Lev3')
    lo=max(low[4],high[4]); hi=min(bounds(x[2]['C'])[0] for x in [low,high])-50
    t=np.linspace(lo,hi,2001); mass=high[3]
    vals=[]
    for data in [low,high]:
        a=evaluate(data[2]['A'],t,mass); b=evaluate(data[2]['B'],t,mass)
        vals.append((np.linalg.norm(a[:,:3]-b[:,:3],axis=1),a,b))
    return {'simulation':'SXS:BBH:0389','levels':['Lev2','Lev3'], 'samples':len(t),
        'alignment':'shared original coordinate time; no phase/time fit; both normalized by Lev3 reference total mass',
        'source_time_interval':[lo,hi],
        'separation_abs_difference_max_M':float(np.max(np.abs(vals[0][0]-vals[1][0]))),
        'position_vector_difference_max_M':float(max(np.max(np.linalg.norm(vals[0][j][:,:3]-vals[1][j][:,:3],axis=1)) for j in [1,2])),
        'spin_vector_difference_max':float(max(np.max(np.linalg.norm(vals[0][j][:,4:]-vals[1][j][:,4:],axis=1)) for j in [1,2])),
        'common_horizon_time_difference_source_units':bounds(high[2]['C'])[0]-bounds(low[2]['C'])[0],
        'interpretation':'Two-resolution coordinate diagnostic only; no convergence order, waveform mismatch or observational error inferred'}

if __name__ == '__main__':
    catalog=[]; reports=[]
    for sim,level in DATASETS:
        entry,report=adapt(ROOT/'sxs-data'/sim/level); catalog.append(entry); reports.append(report)
    report={'adapter':'sxs-adapter/2','datasets':reports,'resolution_comparison':resolution_comparison()}
    write(PUBLIC/'catalog.json',catalog)
    write(ROOT/'adapter-results'/'sxs-validation.json',report)
    print(json.dumps(report,indent=2))
