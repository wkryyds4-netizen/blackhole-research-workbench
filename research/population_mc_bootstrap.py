"""Paired empirical finite-sample diagnostics, NOT a new population inference.

Run using root requirements (astropy). Does not modify frozen input products.
The injection bootstrap assumes IID empirical trials, an unaudited assumption.
"""
import os
os.environ.setdefault('OPENBLAS_NUM_THREADS', '1')
from pathlib import Path
import argparse
import hashlib
import json
import sys
import time

import numpy as np
from scipy.integrate import cumulative_trapezoid

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
from population_functionals import aligned_cdf

KEYS = [f'{si}-{ai}-0-1' for ai in [0, 1] for si in range(3)]
SIGMAS = [.3, .6, 1.2]*2
MODES = ['events', 'injections', 'both']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def prepare_weights(data):
    """Recompute from raw samples through the ORIGINAL model functions."""
    import model
    model.ZGRID = np.linspace(0, 2.3, 3001)
    model.ZDV = model.Planck15.differential_comoving_volume(model.ZGRID).value
    model.ZNORM = {lam: np.trapezoid(model.ZDV*(1+model.ZGRID)**(lam-1), model.ZGRID) for lam in model.LAMBDAS}
    names = ['mass_1', 'mass_ratio', 'redshift', 'a_1', 'a_2', 'cos_tilt_1', 'cos_tilt_2', 'prior']
    errors = []
    weights = []
    with np.load(ROOT/'data/posteriors.npz', allow_pickle=False) as archive:
        for ei, event in enumerate(data['events']):
            draws = {name: archive[f'{event["index"]}_{name}'] for name in names}
            ws = []
            for key in KEYS:
                si, ai, mi, zi = map(int, key.split('-'))
                a, b = model.component_weights(draws, ai, mi, model.LAMBDAS[zi], model.SIGMAS[si])
                moments = model.reduce_weights(a, b)
                reference = np.array(data['scenarios'][key]['event_moments'][ei])
                errors.append(float(np.max(abs(moments/reference-1))))
                assert len(a) == data['scenarios'][key]['sample_counts'][ei]
                ws.append(np.column_stack([a, b]))
            weights.append(np.stack(ws, axis=1))  # sample, model, component
    with np.load(ROOT/'data/injections.npz', allow_pickle=False) as archive:
        draws = dict(archive)
    total = int(draws['total_generated'])
    ws = []
    for key in KEYS:
        si, ai, mi, zi = map(int, key.split('-'))
        a, b = model.component_weights(draws, ai, mi, model.LAMBDAS[zi], model.SIGMAS[si])
        moments = np.array([a.sum()/total, b.sum()/total, (a*a).sum()/total**2,
                            (b*b).sum()/total**2, (a*b).sum()/total**2])
        reference = np.array(data['scenarios'][key]['selection'])
        errors.append(float(np.max(abs(moments/reference-1))))
        ws.append(np.column_stack([a, b]))
    assert max(errors) < 1e-11, max(errors)
    return weights, np.stack(ws, axis=1), total, max(errors)


def empirical_counts(rng, retained, total, repeats):
    """Collapse unrecorded zero-weight trials; keep their multinomial covariance."""
    if not (0 < retained <= total):
        raise ValueError('require 0 < retained <= total')
    if total == retained:
        return rng.multinomial(total, np.full(retained, 1/total), size=repeats)
    p = np.r_[np.full(retained, 1/total), 1-retained/total]
    return rng.multinomial(total, p, size=repeats)[:, :retained]


def mean_covariance(weights, total):
    """Exact covariance of an empirical-bootstrap sample mean, including zeros."""
    means = weights.sum(axis=0)/total
    return weights.T @ weights / total**2 - np.outer(means, means)/total


def quadrature_weights(x):
    dx = np.diff(x)
    return np.r_[dx[0]/2, (dx[:-1]+dx[1:])/2, dx[-1]/2]


def r_grid(points):
    r = np.linspace(float(aligned_cdf(0, 1.2)), .5, points)
    fs = np.array([(.5-r)/(.5-float(aligned_cdf(0, s))) for s in SIGMAS])
    return r, fs


def log_curves(event_means, selection, fs):
    """Inputs: event_means [B,E,K,2], selection [B,K,2]; output [B,K,G]."""
    f = fs[None, :, :]
    logs = np.zeros((event_means.shape[0], fs.shape[0], fs.shape[1]))
    for ei in range(event_means.shape[1]):
        a, b = event_means[:, ei, :, 0, None], event_means[:, ei, :, 1, None]
        logs += np.log(f*a+(1-f)*b)
    a, b = selection[:, :, 0, None], selection[:, :, 1, None]
    logs -= event_means.shape[1]*np.log(f*a+(1-f)*b)
    return logs


def curve_stats(r, logs):
    densities = np.exp(logs-logs.max(axis=-1, keepdims=True))
    densities /= np.trapezoid(densities, r, axis=-1)[..., None]
    c = cumulative_trapezoid(densities, r, axis=-1, initial=0)
    c /= c[..., -1, None]
    qs = np.empty(logs.shape[:2]+(3,))
    for bi in range(logs.shape[0]):
        for ki in range(logs.shape[1]):
            qs[bi, ki] = np.interp([.05, .5, .95], c[bi, ki], r)
    return qs, densities


def delta_shape_diagnostics(weights, injection, total, fs, pquad):
    """2 E_p Var(logL) - 2 Var(E_p logL), with same-draw covariance."""
    k, g = fs.shape
    point = {'events': np.zeros((k, g)), 'injections': np.zeros((k, g))}
    shape = {'events': np.zeros(k), 'injections': np.zeros(k)}
    for w, nt, mode, factor in [(w, len(w), 'events', 1) for w in weights]+[(injection, total, 'injections', -len(weights))]:
        for ki in range(k):
            means = w[:, ki].sum(axis=0)/nt
            covariance = mean_covariance(w[:, ki], nt)
            jac = factor*np.column_stack([fs[ki], 1-fs[ki]])/(fs[ki]*means[0]+(1-fs[ki])*means[1])[:, None]
            variance = np.einsum('gi,ij,gj->g', jac, covariance, jac)
            mean_jac = pquad[ki] @ jac
            point[mode][ki] += variance
            shape[mode][ki] += 2*(pquad[ki]@variance-mean_jac@covariance@mean_jac)
    return {KEYS[ki]: {mode: dict(posterior_mean_pointwise_variance=float(pquad[ki]@point[mode][ki]),
                                  posterior_mean_pair_difference_variance=float(shape[mode][ki])) for mode in ['events', 'injections']}
            for ki in range(k)}


def generate_means(weights, injection, total, repeats, seed):
    rng_event, rng_inj = [np.random.default_rng(s) for s in np.random.SeedSequence(seed).spawn(2)]
    k = injection.shape[1]
    event_means = np.empty((repeats, len(weights), k, 2))
    for ei, w in enumerate(weights):
        for start in range(0, repeats, 64):
            stop = min(start+64, repeats)
            counts = empirical_counts(rng_event, len(w), len(w), stop-start)
            event_means[start:stop, ei] = (counts @ w.reshape(len(w), -1)).reshape(stop-start, k, 2)/len(w)
    selection = np.empty((repeats, k, 2))
    for start in range(0, repeats, 64):
        stop = min(start+64, repeats)
        counts = empirical_counts(rng_inj, len(injection), total, stop-start)
        selection[start:stop] = (counts @ injection.reshape(len(injection), -1)).reshape(stop-start, k, 2)/total
    return event_means, selection


def evaluate_replicates(r, fs, event_means, selection, original_events, original_selection, pquad):
    n = len(event_means)
    summaries, shape_stats = {}, {}
    for mode in MODES:
        qs = []
        # Center each random logL curve using the FIXED nominal posterior.
        # This removes arbitrary random offsets before the shape variance.
        sum_center = np.zeros(fs.shape)
        square_center = np.zeros(fs.shape)
        for start in range(0, n, 32):
            stop = min(start+32, n)
            ev = event_means[start:stop] if mode != 'injections' else np.broadcast_to(original_events, (stop-start,)+original_events.shape[1:])
            se = selection[start:stop] if mode != 'events' else np.broadcast_to(original_selection, (stop-start,)+original_selection.shape[1:])
            logs = log_curves(ev, se, fs)
            qs.append(curve_stats(r, logs)[0])
            centered = logs-np.sum(logs*pquad[None, :, :], axis=-1, keepdims=True)
            sum_center += centered.sum(axis=0)
            square_center += (centered*centered).sum(axis=0)
        variance = (square_center-sum_center**2/n)/(n-1)
        shape_stats[mode] = (2*np.sum(pquad*variance, axis=-1)).tolist()
        summaries[mode] = np.concatenate(qs)
    return summaries, shape_stats


def describe(draws):
    return dict(mean=float(np.mean(draws)), sd=float(np.std(draws, ddof=1)),
                quantiles=np.quantile(draws, [.05, .5, .95]).tolist())


def summarize_bootstrap(qs, nominal, nseed):
    rows = {}
    for ki, key in enumerate(KEYS):
        rows[key] = dict(nominal_quantiles=nominal[ki].tolist(), modes={})
        for mode in MODES:
            d = qs[mode][:, ki, 1]
            rows[key]['modes'][mode] = describe(d)
            rows[key]['modes'][mode]['mean_minus_nominal_median'] = float(d.mean()-nominal[ki, 1])
            rows[key]['modes'][mode]['seed_sds'] = [float(np.std(d[i*nseed:(i+1)*nseed], ddof=1)) for i in range(2)]
    pairs = []
    for ai in range(2):
        for si, sj in [(0, 1), (0, 2), (1, 2)]:
            i, j = ai*3+si, ai*3+sj
            pair = dict(keys=[KEYS[i], KEYS[j]], nominal_difference=float(nominal[i, 1]-nominal[j, 1]), modes={})
            for mode in MODES:
                x, y = qs[mode][:, i, 1], qs[mode][:, j, 1]
                pair['modes'][mode] = describe(x-y)
                pair['modes'][mode]['correlation'] = float(np.corrcoef(x, y)[0, 1])
                pair['modes'][mode]['sd_if_covariance_ignored'] = float(np.sqrt(np.var(x, ddof=1)+np.var(y, ddof=1)))
            pairs.append(pair)
    return rows, pairs


def make_figure(result, output):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'svg.hashsalt': 'population-stage2', 'axes.spines.top': False, 'axes.spines.right': False})
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    colors = ['#0072B2', '#D55E00', '#009E73']
    for ai, ax in enumerate(axes):
        pairs = result['pairs'][ai*3:(ai+1)*3]
        for pi, pair in enumerate(pairs):
            ax.scatter([pair['nominal_difference']*100], [pi], marker='D', color='black', s=20, zorder=4)
            for mi, mode in enumerate(MODES):
                lo, median, hi = np.array(pair['modes'][mode]['quantiles'])*100
                y = pi+(mi-1)*.18
                ax.errorbar(median, y, xerr=[[median-lo], [hi-median]], fmt='o', color=colors[mi],
                            markersize=4, capsize=3, label=mode if pi == 0 else None)
        ax.axvline(0, color='gray', linestyle='--', linewidth=1)
        ax.set_yticks([0, 1, 2], ['0.3 minus 0.6', '0.3 minus 1.2', '0.6 minus 1.2'])
        ax.set_ylabel('Tilt-width pair')
        ax.set_xlabel('Difference of posterior r medians (percentage points)')
        ax.set_title(['Low-spin assumption: Beta(1.5, 8.5)', 'Baseline amplitude: Beta(2, 5)'][ai])
        ax.legend(loc='best', fontsize=9)
    fig.suptitle('Paired empirical bootstrap | bars: 5-95% resampling ranges | diamonds: nominal\nFixed 69-event catalog, shared r prior; not astrophysical credible intervals', fontsize=11)
    fig.tight_layout()
    fig.savefig(output/'paired-mc-comparison.svg', metadata={'Date': None})
    fig.savefig(output/'paired-mc-comparison.png', dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--repeats-per-seed', type=int, default=1024)
    parser.add_argument('--points', type=int, default=1001)
    parser.add_argument('--output', type=Path, default=ROOT/'research/population-stage2')
    args = parser.parse_args()
    if args.repeats_per_seed < 2 or args.points < 101:
        parser.error('require >=2 repeats per seed and >=101 grid points')
    start = time.time()
    args.output.mkdir(exist_ok=True, parents=True)
    data = json.loads((ROOT/'results/model.json').read_text(encoding='utf-8'))
    weights, injection, total, error = prepare_weights(data)
    print(f'Raw moments verified: max relative error {error:.3g}', flush=True)
    original_events = np.stack([w.mean(axis=0) for w in weights])[None, ...]
    original_selection = (injection.sum(axis=0)/total)[None, ...]
    r, fs = r_grid(args.points)
    nominal, density = curve_stats(r, log_curves(original_events, original_selection, fs))
    pquad = density[0]*quadrature_weights(r)
    delta = delta_shape_diagnostics(weights, injection, total, fs, pquad)
    ev, se = [], []
    for seed in [20261024, 20261025]:
        print(f'Resampling seed {seed} ({args.repeats_per_seed} repeats)', flush=True)
        e, s = generate_means(weights, injection, total, args.repeats_per_seed, seed)
        ev.append(e)
        se.append(s)
    event_means, selection = np.concatenate(ev), np.concatenate(se)
    np.savez_compressed(args.output/'bootstrap-means.npz', event_means=event_means, selection=selection,
                        original_events=original_events, original_selection=original_selection)
    qs, shape = evaluate_replicates(r, fs, event_means, selection, original_events, original_selection, pquad)
    rows, pairs = summarize_bootstrap(qs, nominal[0], args.repeats_per_seed)
    for ki, key in enumerate(KEYS):
        rows[key]['delta_method'] = delta[key]
        rows[key]['empirical_shape_variance'] = {mode: shape[mode][ki] for mode in MODES}
    # A fixed, evenly-spaced subset tests numerical integration, not another MC run.
    idx = np.linspace(0, len(event_means)-1, min(64, len(event_means))).astype(int)
    r2, fs2 = r_grid(2*args.points-1)
    nom2, den2 = curve_stats(r2, log_curves(original_events, original_selection, fs2))
    q2, _ = evaluate_replicates(r2, fs2, event_means[idx], selection[idx], original_events, original_selection,
                              den2[0]*quadrature_weights(r2))
    max_grid_error = max(float(np.max(abs(q2[m]-qs[m][idx]))) for m in MODES)
    result = dict(schema='paired-empirical-mc/1', date='2026-10-24', keys=KEYS, seeds=[20261024, 20261025],
                  repeats_per_seed=args.repeats_per_seed, total_repeats=len(event_means), points=args.points,
                  common_r_support=[float(r[0]), float(r[-1])], n_events=len(weights),
                  retained_injections=len(injection), total_generated=total,
                  max_raw_moment_relative_error=error, grid_check_repeats=len(idx), max_grid_quantile_error=max_grid_error,
                  scenarios=rows, pairs=pairs, bootstrap_quantiles={m: qs[m].tolist() for m in MODES},
                  assumptions=['IID empirical draws within each event; chain correlations unavailable',
                               'IID empirical injection trials including zero-weight non-retained trials; original stratification unaudited',
                               'Shared counts across all six models and all parameter points',
                               'Fixed catalog and selection; no correction for unseen tails or mismatched provenance',
                               'Resampling intervals describe numerical sensitivity, not astrophysical posterior uncertainty'],
                  input_sha256={str(p.relative_to(ROOT)).replace('\\', '/'): digest(p) for p in
                                [ROOT/'results/model.json', ROOT/'data/posteriors.npz', ROOT/'data/injections.npz', ROOT/'model.py']},
                  script_sha256=digest(Path(__file__)), python=sys.version, numpy=np.__version__, elapsed_seconds=time.time()-start)
    (args.output/'paired-bootstrap-results.json').write_text(json.dumps(result, indent=2, allow_nan=False), encoding='utf-8')
    make_figure(result, args.output)
    print(json.dumps({k: result[k] for k in ['total_repeats', 'max_raw_moment_relative_error', 'max_grid_quantile_error', 'elapsed_seconds']}), flush=True)
    for key, row in rows.items():
        print(key, 'nominal', row['nominal_quantiles'][1], 'both SD', row['modes']['both']['sd'], 'shape', row['empirical_shape_variance'], flush=True)


if __name__ == '__main__':
    main()
