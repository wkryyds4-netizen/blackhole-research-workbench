"""Held-out coordinate alignment diagnostic for pinned SXS:BBH:0389.

Run: simulator/.adapters-venv/Scripts/python.exe research/alignment_study.py
This does not calculate waveform mismatch, convergence order, or observables.
"""
from pathlib import Path
import hashlib
import json
import platform

import numpy as np
import scipy
from scipy.optimize import minimize_scalar
from sxs_adapter import read, bounds

ROOT = Path(__file__).resolve().parent
OUT = ROOT / 'alignment-results'


def rotate_z(v, phi):
    c, s = np.cos(phi), np.sin(phi)
    return v @ np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]]).T


class RelativeTrajectory:
    """Original inertial coordinates, A minus B; one shared mass scale."""
    def __init__(self, data, mass):
        self.data, self.mass = data, mass
        self.lo = max(data[4], *(bounds(data[2][k])[0] for k in 'AB')) / mass
        self.hi = min(*(bounds(data[2][k])[1] for k in 'AB'), bounds(data[2]['C'])[0]) / mass

    def __call__(self, t):
        t = np.asarray(t)
        if t.min() < self.lo or t.max() > self.hi:
            raise ValueError('Coordinate extrapolation forbidden')
        def body(k):
            a = self.data[2][k]['position']
            return np.column_stack([np.interp(t*self.mass, a[:, 0], a[:, j]) for j in (1, 2, 3)]) / self.mass
        return body('A') - body('B')


def fit_fixed_alignment(reference, candidate, train, bound=30., samples=1001):
    t = np.linspace(*train, samples)
    target = reference(t)
    def objective(dt, details=False):
        v = candidate(t + dt)
        # Minimize sum |Rz(phi) v - target|^2 analytically in phi.
        dot = np.sum(v[:, 0]*target[:, 0] + v[:, 1]*target[:, 1])
        cross = np.sum(v[:, 0]*target[:, 1] - v[:, 1]*target[:, 0])
        phi = np.arctan2(cross, dot)
        mse = float(np.mean(np.sum((rotate_z(v, phi)-target)**2, axis=1)))
        return (mse, float(phi)) if details else mse
    coarse = np.linspace(-bound, bound, 601)
    values = np.array([objective(dt) for dt in coarse])
    # Refine every coarse local minimum, then select the global best sampled basin.
    choices = [(float(values[0]), float(coarse[0])), (float(values[-1]), float(coarse[-1]))]
    for i in range(1, len(coarse)-1):
        if values[i] <= values[i-1] and values[i] <= values[i+1]:
            result = minimize_scalar(objective, bounds=(coarse[i-1], coarse[i+1]), method='bounded', options={'xatol':1e-10})
            choices.append((float(result.fun), float(result.x)))
    mse, dt = min(choices)
    _, phi = objective(dt, True)
    return {'delta_t_M':dt, 'phase_rotation_rad':phi, 'training_rms_relative_position_M':float(np.sqrt(mse)),
            'search_bound_M':bound, 'search_boundary_hit':abs(dt) > bound-1e-5, 'fit_samples':samples}


def metrics(reference, candidate, t, dt, phi):
    a, b = reference(t), rotate_z(candidate(t+dt), phi)
    residual = np.linalg.norm(b-a, axis=1)
    separation = np.abs(np.linalg.norm(b, axis=1)-np.linalg.norm(a, axis=1))
    phase = np.arctan2(a[:,0]*b[:,1]-a[:,1]*b[:,0], a[:,0]*b[:,0]+a[:,1]*b[:,1])
    summary = {'rms_relative_position_M':float(np.sqrt(np.mean(residual**2))), 'max_relative_position_M':float(residual.max()),
               'rms_separation_difference_M':float(np.sqrt(np.mean(separation**2))), 'max_separation_difference_M':float(separation.max()),
               'rms_planar_phase_residual_rad':float(np.sqrt(np.mean(phase**2))), 'max_abs_planar_phase_residual_rad':float(np.max(np.abs(phase)))}
    return summary, {'relative_position_difference_M':residual.tolist(), 'separation_difference_M':separation.tolist(), 'planar_phase_residual_rad':phase.tolist()}


def checks():
    queried = {'reference':[], 'candidate':[]}
    def reference(t):
        t = np.asarray(t)
        phase = .07*t + .0008*t*t
        radius = 10 - .015*t
        return np.column_stack([radius*np.cos(phase), radius*np.sin(phase), .03*t])
    dt, phi = 3.25, -.43
    def candidate(t):
        return rotate_z(reference(np.asarray(t)-dt), -phi)
    def tracked_reference(t):
        queried['reference'].extend([float(np.min(t)), float(np.max(t))])
        return reference(t)
    def tracked_candidate(t):
        queried['candidate'].extend([float(np.min(t)), float(np.max(t))])
        return candidate(t)
    fit = fit_fixed_alignment(tracked_reference, tracked_candidate, [20, 80], 10)
    assert min(queried['reference']) == 20 and max(queried['reference']) == 80
    assert min(queried['candidate']) >= 10 and max(queried['candidate']) <= 90
    assert abs(fit['delta_t_M']-dt) < 1e-5
    assert abs(fit['phase_rotation_rad']-phi) < 1e-6
    summary, _ = metrics(reference, candidate, np.linspace(90,120,301), fit['delta_t_M'], fit['phase_rotation_rad'])
    assert summary['max_relative_position_M'] < 1e-5
    np.testing.assert_allclose(rotate_z(np.array([[1.,0,0]]), np.pi/2), [[0,1,0]], atol=1e-14)
    np.testing.assert_allclose(np.linalg.norm(rotate_z(reference(np.arange(5)), .7), axis=1), np.linalg.norm(reference(np.arange(5)), axis=1))
    return {'synthetic_known_delta_t_M':dt, 'synthetic_recovered_delta_t_M':fit['delta_t_M'],
            'synthetic_known_phi_rad':phi, 'synthetic_recovered_phi_rad':fit['phase_rotation_rad'],
            'synthetic_held_out_max_error_M':summary['max_relative_position_M'], 'positive_rotation_and_norm_preservation':True,
            'fit_query_audit':{'reference_range':[min(queried['reference']),max(queried['reference'])],
                               'candidate_range':[min(queried['candidate']),max(queried['candidate'])], 'no_validation_values_requested':True}}


def main():
    validation = checks()
    low = read(ROOT/'sxs-data/0389/Lev2')
    high = read(ROOT/'sxs-data/0389/Lev3')
    mass = high[3]
    reference, candidate = RelativeTrajectory(high, mass), RelativeTrajectory(low, mass)
    bound = 30.
    # The SAME fixed evaluation grid is safe for every allowed shift.
    start = max(reference.lo, candidate.lo+bound)
    end = min(reference.hi, candidate.hi-bound)-50.
    span = end-start
    validation_window = [start+.65*span, end]
    evaluation_t = np.linspace(*validation_window, 2001)
    train = [start, start+.25*span]
    assert train[1]+bound < validation_window[0]-bound
    dt_ch = (bounds(low[2]['C'])[0]-bounds(high[2]['C'])[0])/mass
    fitted = fit_fixed_alignment(reference, candidate, train, bound)
    methods = []
    for id_, label, fit in [
        ('raw', '共同原始坐标时间', {'delta_t_M':0., 'phase_rotation_rad':0.}),
        ('common-horizon', '共同视界时间平移', {'delta_t_M':dt_ch, 'phase_rotation_rad':0.}),
        ('early-fit', '早期窗口固定时间＋相位拟合', fitted)]:
        summary, series = metrics(reference, candidate, evaluation_t, fit['delta_t_M'], fit['phase_rotation_rad'])
        training_summary, _ = metrics(reference, candidate, np.linspace(*train, 1001), fit['delta_t_M'], fit['phase_rotation_rad'])
        methods.append({'id':id_, 'label':label, **fit, 'training':training_summary, 'validation':summary, 'series':series})
    sensitivity = []
    for fraction in (.15, .25, .35):
        window = [start, start+fraction*span]
        assert window[1]+bound < validation_window[0]-bound
        fit = fit_fixed_alignment(reference, candidate, window, bound)
        summary, _ = metrics(reference, candidate, evaluation_t, fit['delta_t_M'], fit['phase_rotation_rad'])
        sensitivity.append({'training_fraction':fraction, 'training_window_M':window, **fit, 'validation':summary})
    dense = fit_fixed_alignment(reference, candidate, train, bound, 2001)
    dense_summary, _ = metrics(reference, candidate, evaluation_t, dense['delta_t_M'], dense['phase_rotation_rad'])
    validation['fit_sample_count_sensitivity'] = {'coarse_samples':1001, 'dense_samples':2001,
        'delta_t_change_M':dense['delta_t_M']-fitted['delta_t_M'], 'phase_change_rad':dense['phase_rotation_rad']-fitted['phase_rotation_rad'],
        'held_out_rms_change_M':dense_summary['rms_relative_position_M']-methods[2]['validation']['rms_relative_position_M']}
    for out_of_bounds in ([candidate.lo-1], [candidate.hi+1]):
        try:
            candidate(out_of_bounds)
        except ValueError:
            pass
        else:
            raise AssertionError('Extrapolation was accepted')
    validation.update({'extrapolation_rejected':True, 'no_training_validation_overlap_even_with_shift_bound':True,
                       'all_fitted_searches_interior':not any(f['search_boundary_hit'] for f in sensitivity)})
    assert validation['all_fitted_searches_interior']
    sources = [{'level':data[1]['level'], 'doi':data[1]['doi'], 'files':data[1]['files']} for data in [low, high]]
    report = {'schema':'blackhole-alignment-study/1', 'simulation':'SXS:BBH:0389', 'reference_level':'Lev3', 'comparison_level':'Lev2',
        'reference_total_mass_source_units':mass, 'time_units':'source coordinate time / Lev3 reference total mass (GM/c^3)',
        'length_units':'GM/c^2 with the same Lev3 reference mass',
        'convention':'r = center(A) - center(B); compare Rz(phi) r_Lev2(t + delta_t) with r_Lev3(t)',
        'rotation_axis':'original inertial +z; maps to workbench -y under Rx(+pi/2)',
        'fit_objective':'mean squared 3D relative-position residual on early window; constant proper rotation and constant time shift only',
        'common_horizon_alignment_note':'uses a late-time landmark and is a retrospective comparison, not an early-only predictive fit',
        'interpolation':'linear on each horizon original native grid, no display compression and no extrapolation',
        'training_window_M':train, 'validation_window_M':validation_window, 'validation_samples':len(evaluation_t),
        'source_common_horizon_difference_high_minus_low':bounds(high[2]['C'])[0]-bounds(low[2]['C'])[0],
        'time_M':evaluation_t.tolist(), 'methods':methods, 'training_window_sensitivity':sensitivity, 'checks':validation,
        'limitations':['Coordinate/gauge dependent; coordinate separation is not proper distance.',
            'Two resolutions do not establish convergence order or a continuum error.',
            'No waveform, detector weighting, noise-weighted mismatch, or formation-channel evidence is evaluated.',
            'Only a near-nonspinning system is studied; constant z rotation is not a general precession alignment.',
            'Linear native-grid interpolation remains in the diagnostic; differences include sampling effects.',
            'Early-window choices and shift bound are protocol choices, not probabilistic uncertainty intervals.'],
        'references':[{'url':'https://arxiv.org/abs/0810.1767','role':'Published waveform comparison demonstrates alignment changes quoted numerical differences; not validation of this coordinate metric.'},
            {'url':'https://sxs.readthedocs.io/en/main/tutorials/03-Horizons/','role':'Official horizon quantities and coordinate trajectories.'}],
        'sources':sources, 'runtime':{'python':platform.python_version(),'numpy':np.__version__,'scipy':scipy.__version__},
        'script_sha256':hashlib.sha256(Path(__file__).read_bytes()).hexdigest()}
    OUT.mkdir(exist_ok=True)
    (OUT/'alignment.json').write_text(json.dumps(report, ensure_ascii=False, separators=(',',':'), allow_nan=False), encoding='utf-8')
    compact = {k:v for k,v in report.items() if k not in ('time_M','methods')}
    compact['methods']=[{k:v for k,v in method.items() if k != 'series'} for method in methods]
    (OUT/'summary.json').write_text(json.dumps(compact, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    print(json.dumps({'methods':compact['methods'],'checks':validation,'training_window_sensitivity':sensitivity},ensure_ascii=False,indent=2))


if __name__ == '__main__':
    main()
