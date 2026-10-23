"""Conditional spin-tilt functional study, with explicit common-support control.

Run with the adapter environment. Only saved likelihood moments are needed.
No MC screen is used to truncate or renormalize any posterior.
"""
from pathlib import Path
import argparse
import hashlib
import json
import sys

import numpy as np
from scipy.integrate import cumulative_trapezoid
from scipy.special import ndtr, xlogy

ROOT = Path(__file__).resolve().parents[1]
SIGMAS = (.3, .6, 1.2)
Q = np.array([.05, .5, .95])


def aligned_cdf(x, sigma):
    return (ndtr((np.asarray(x)-1)/sigma)-ndtr(-2/sigma))/(ndtr(0)-ndtr(-2/sigma))


def functional_coefficients(sigma):
    u = (1-np.cos(np.pi/18))/2
    return {'r': (.5, float(aligned_cdf(0, sigma))-.5),
            'a': (u, float(1-aligned_cdf(np.cos(np.pi/18), sigma))-u)}


def log_likelihood(rec, f, correct=True):
    f = np.asarray(f)
    em = np.asarray(rec['event_moments'])
    mix = f[:, None]*em[:, 0]+(1-f[:, None])*em[:, 1]
    logs = np.log(mix).sum(axis=1)
    if correct:
        a, b = rec['selection'][:2]
        logs -= len(em)*np.log(f*a+(1-f)*b)
    return logs


def normalize(x, logs):
    density = np.exp(logs-np.max(logs))
    return density/np.trapezoid(density, x)


def cdf(x, p):
    c = cumulative_trapezoid(p, x, initial=0)
    return c/c[-1]


def summarize(x, p):
    p = p/np.trapezoid(p, x)
    mean = np.trapezoid(x*p, x)
    sd = np.sqrt(np.trapezoid((x-mean)**2*p, x))
    width = x[-1]-x[0]
    prior_sd = width/np.sqrt(12)
    return dict(quantiles=np.interp(Q, cdf(x, p), x).tolist(), mean=float(mean),
                sd=float(sd), prior_support=[float(x[0]), float(x[-1])],
                kl_to_uniform_prior_nats=float(np.trapezoid(xlogy(p, p*width), x)),
                standardized_shift=float(abs(mean-(x[0]+x[-1])/2)/prior_sd),
                contraction=float(sd/prior_sd))


def transform(x, p, offset, slope):
    y = offset+slope*x
    yp = p/abs(slope)
    return (y[::-1], yp[::-1]) if slope < 0 else (y, yp)


def mc_diagnostics(rec, f):
    """Pointwise delta-method diagnostics, not a bias or covariance estimate."""
    f = np.asarray(f)[:, None]
    em = np.asarray(rec['event_moments'])
    counts = np.asarray(rec['sample_counts'])
    n = len(em)
    mean = f*em[:, 0]+(1-f)*em[:, 1]
    second = f*f*em[:, 2]+(1-f)**2*em[:, 3]+2*f*(1-f)*em[:, 4]
    all_event_ess = counts*mean**2/second
    event_ess = np.min(all_event_ess, axis=1)
    event_var = np.maximum(0, ((second/mean**2-1)/counts).sum(axis=1))
    f = f[:, 0]
    a, b, aa, bb, ab = rec['selection']
    xi = f*a+(1-f)*b
    second_xi = f*f*aa+(1-f)**2*bb+2*f*(1-f)*ab
    inj_ess = xi*xi/second_xi
    sd = np.sqrt(event_var+n*n*np.maximum(0, second_xi/xi**2-1/rec['total_generated']))
    return dict(passed=(event_ess >= 50) & (inj_ess >= 4*n) & (sd <= 1),
                min_event_ess=event_ess, min_event_index=np.argmin(all_event_ess, axis=1),
                injection_ess=inj_ess, log_likelihood_mc_sd=sd)


def mc_screen(rec, f):
    return mc_diagnostics(rec, f)['passed']


def screen_masses(x, p, mask):
    return dict(prior=float(np.trapezoid(mask.astype(float), x)/(x[-1]-x[0])),
                posterior=float(np.trapezoid(p*mask, x)))


def run_study(data, points=4001):
    f = np.linspace(0, 1, points)
    lower = max(float(aligned_cdf(0, s)) for s in SIGMAS)
    r = np.linspace(lower, .5, points)
    records, common_densities, common_masks = {}, {}, {}
    for key, rec in data['scenarios'].items():
        sigma = SIGMAS[int(key.split('-')[0])]
        coefficients = functional_coefficients(sigma)
        p = normalize(f, log_likelihood(rec, f))
        native = {'f': summarize(f, p)}
        for quantity, (offset, slope) in coefficients.items():
            native[quantity] = summarize(*transform(f, p, offset, slope))
        fp = (.5-r)/(.5-float(aligned_cdf(0, sigma)))
        # Ascending r implies descending f: evaluate L directly on ascending r.
        pc = normalize(r, log_likelihood(rec, fp))
        fmax = float(fp[0])
        common = summarize(r, pc)
        common['retained_original_prior_mass'] = fmax
        common['retained_original_posterior_mass'] = float(np.interp(fmax, f, cdf(f, p)))
        diagnostics = mc_diagnostics(rec, fp)
        common_masks[key] = diagnostics['passed']
        indices, counts = np.unique(diagnostics['min_event_index'], return_counts=True)
        common['mc_diagnostic_ranges'] = {name: [float(diagnostics[name].min()), float(diagnostics[name].max())]
                                          for name in ['min_event_ess', 'injection_ess', 'log_likelihood_mc_sd']}
        common['minimum_ess_event_grid_counts'] = {data['events'][i]['name']: int(count)
                                                  for i, count in zip(indices, counts)}
        common['mc_passing_mass'] = screen_masses(r, pc, common_masks[key])
        native['mc_passing_mass'] = screen_masses(f, p, mc_screen(rec, f))
        records[key] = dict(sigma=sigma, native=native, common_r=common)
        common_densities[key] = pc
    groups = []
    for ai in range(3):
        for mi in range(3):
            for zi in range(3):
                keys = [f'{si}-{ai}-{mi}-{zi}' for si in range(3)]
                common_pass = np.logical_and.reduce([common_masks[k] for k in keys])
                original = [records[k]['native']['r']['quantiles'][1] for k in keys]
                matched = [records[k]['common_r']['quantiles'][1] for k in keys]
                pairwise = {f'{i}:{j}': float(.5*np.trapezoid(abs(common_densities[keys[i]]-
                             common_densities[keys[j]]), r)) for i, j in [(0, 1), (0, 2), (1, 2)]}
                groups.append(dict(keys=keys, native_r_median_span=float(np.ptp(original)),
                                   common_r_median_span=float(np.ptp(matched)), common_r_tv=pairwise,
                                   common_mc_intersection_prior_mass=float(np.trapezoid(common_pass.astype(float), r)/(r[-1]-r[0])),
                                   common_mc_intersection_posterior_mass={k: float(np.trapezoid(common_pass*common_densities[k], r)) for k in keys}))
    default = data['default']
    p0 = normalize(f, log_likelihood(data['scenarios'][default], f, correct=False))
    base_sigma = records[default]['sigma']
    baseline_uncorrected = {'f': summarize(f, p0), 'r': summarize(*transform(f, p0, *functional_coefficients(base_sigma)['r']))}
    base_keys = [f'{si}-'+default.split('-', 1)[1] for si in range(3)]
    curves = dict(common_r_grid=r.tolist(), baseline_triplet={k: common_densities[k].tolist() for k in base_keys})
    invariant_error = max(abs(row['native'][x][metric]-row['native']['f'][metric])
                          for row in records.values() for x in ['r', 'a']
                          for metric in ['kl_to_uniform_prior_nats', 'standardized_shift', 'contraction'])
    summary = dict(default=default, scenarios=len(records), matched_groups=len(groups),
                   native_f_median_range=[min(v['native']['f']['quantiles'][1] for v in records.values()), max(v['native']['f']['quantiles'][1] for v in records.values())],
                   native_r_median_range=[min(v['native']['r']['quantiles'][1] for v in records.values()), max(v['native']['r']['quantiles'][1] for v in records.values())],
                   baseline=records[default], baseline_without_selection=baseline_uncorrected,
                   max_affine_invariance_error=invariant_error,
                   groups_with_smaller_median_span=sum(g['common_r_median_span'] < g['native_r_median_span'] for g in groups),
                   common_tv_max_range=[min(max(g['common_r_tv'].values()) for g in groups), max(max(g['common_r_tv'].values()) for g in groups)],
                   min_retained_original_posterior_mass=min(v['common_r']['retained_original_posterior_mass'] for v in records.values()),
                   groups_with_empty_mc_intersection=sum(g['common_mc_intersection_prior_mass'] == 0 for g in groups))
    return dict(schema='population-functional-study/1', points=points, common_r_support=[lower, .5],
                summary=summary, scenarios=records, groups=groups, curves=curves)


def make_figures(data, result, dest):
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    plt.rcParams.update({'font.size': 10, 'axes.spines.top': False, 'axes.spines.right': False, 'svg.hashsalt': 'population-functional-study-v1'})
    colors = ['#0072B2', '#D55E00', '#009E73']
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    f = np.linspace(0, 1, result['points'])
    baseline = result['summary']['default']
    for si, sigma in enumerate(SIGMAS):
        key = f'{si}-'+baseline.split('-', 1)[1]
        p = normalize(f, log_likelihood(data['scenarios'][key], f))
        x, y = transform(f, p, *functional_coefficients(sigma)['r'])
        axes[0].plot(x, y, color=colors[si], label=fr'$\sigma={sigma}$')
        axes[0].hlines(1/(x[-1]-x[0]), x[0], x[-1], color=colors[si], linestyle='--', alpha=.6)
        axes[1].plot(result['curves']['common_r_grid'], result['curves']['baseline_triplet'][key], color=colors[si])
    lo, hi = result['common_r_support']
    axes[1].hlines(1/(hi-lo), lo, hi, color='gray', linestyle='--', label='Common prior')
    for ax in axes[:2]:
        ax.set_xlabel(r'Intrinsic primary-BH fraction $r=P(\cos\theta_1<0)$')
        ax.set_ylabel('Density')
    axes[0].legend()
    axes[1].legend()
    axes[0].set_title('Same f prior; different r supports\nDashed: induced priors')
    axes[1].set_title('Uniform r on shared support\nA restricted-prior diagnostic')
    groups = result['groups']
    xx = np.array([g['native_r_median_span']*100 for g in groups])
    yy = np.array([g['common_r_median_span']*100 for g in groups])
    empty = np.array([g['common_mc_intersection_prior_mass'] == 0 for g in groups])
    axes[2].scatter(xx[~empty], yy[~empty], c='#55479a', s=28, label='Nonempty shared MC-pass domain')
    axes[2].scatter(xx[empty], yy[empty], edgecolors='#b75532', facecolors='none', s=36,
                    label='Empty shared MC-pass domain')
    extent = max(max(xx), max(yy))*1.08
    axes[2].plot([0, extent], [0, extent], '--', color='gray')
    axes[2].set(xlim=(0, extent), ylim=(0, extent), xlabel='Original-prior median span (percentage points)',
                ylabel='Common-prior median span (percentage points)', title='27 matched width triplets\nFixed spin-magnitude, mass and redshift settings')
    axes[2].legend(fontsize=8, loc='upper left')
    fig.suptitle('69-event tutorial baseline | conditional models | no MC-based truncation', fontsize=12)
    fig.tight_layout()
    fig.savefig(dest/'prior-support-control.svg', metadata={'Date': None})
    fig.savefig(dest/'prior-support-control.png', dpi=160)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('--points', type=int, default=4001)
    parser.add_argument('--output', type=Path, default=ROOT/'research/population-stage1')
    args = parser.parse_args()
    if args.points < 101:
        parser.error('--points must be >= 101')
    src = ROOT/'results/model.json'
    data = json.loads(src.read_text(encoding='utf-8'))
    result = run_study(data, args.points)
    result['provenance'] = dict(input='results/model.json', sha256=hashlib.sha256(src.read_bytes()).hexdigest(),
                                script_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                                date='2026-10-24', python=sys.version, numpy=np.__version__)
    result['scope'] = 'Conditional phenomenological tilt mixture. Common-r prior restricts support; no model averaging; no astrophysical formation-channel fractions. MC screen is reported only.'
    args.output.mkdir(parents=True, exist_ok=True)
    (args.output/'functional-results.json').write_text(json.dumps(result, ensure_ascii=False, indent=2, allow_nan=False), encoding='utf-8')
    make_figures(data, result, args.output)
    print(json.dumps(result['summary'], indent=2))


if __name__ == '__main__':
    main()
