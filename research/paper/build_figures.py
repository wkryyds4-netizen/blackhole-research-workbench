"""Create paper-sized original figures from frozen, published numerical outputs."""
from pathlib import Path
import json
import sys

import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
sys.path.insert(0, str(ROOT/'research'))
from population_functionals import functional_coefficients, log_likelihood, normalize, transform


def save(fig, stem):
    fig.tight_layout()
    fig.savefig(HERE/'figures'/f'{stem}.png', dpi=240)
    fig.savefig(HERE/'figures'/f'{stem}.svg', metadata={'Date': None})
    plt.close(fig)


def main():
    (HERE/'figures').mkdir(exist_ok=True)
    plt.rcParams.update({'font.size': 9, 'svg.hashsalt': 'blackhole-paper-v02',
                         'axes.spines.top': False, 'axes.spines.right': False})
    model = json.loads((ROOT/'results/model.json').read_text(encoding='utf-8'))
    stage1 = json.loads((ROOT/'research/population-stage1/functional-results.json').read_text())
    stage2 = json.loads((ROOT/'research/population-stage2/paired-bootstrap-results.json').read_text())
    colors = ['#0072B2', '#D55E00', '#009E73']
    f = np.linspace(0, 1, 4001)
    fig, axes = plt.subplots(1, 2, figsize=(7.1, 3.1))
    for i, sigma in enumerate([.3, .6, 1.2]):
        key = f'{i}-1-0-1'
        p = normalize(f, log_likelihood(model['scenarios'][key], f))
        x, y = transform(f, p, *functional_coefficients(sigma)['r'])
        axes[0].plot(x, y, color=colors[i], label=fr'$\sigma={sigma}$')
        axes[0].hlines(1/np.ptp(x), x[0], x[-1], color=colors[i], linestyle='--', alpha=.65)
        axes[1].plot(stage1['curves']['common_r_grid'], stage1['curves']['baseline_triplet'][key], color=colors[i])
    lo, hi = stage1['common_r_support']
    axes[1].hlines(1/(hi-lo), lo, hi, color='gray', linestyle='--', label='Common prior')
    for ax in axes:
        ax.set_xlabel(r'$r=P(\cos\theta_1<0)$')
        ax.set_ylabel('Density')
        ax.legend(fontsize=8)
    axes[0].set_title('Original induced priors')
    axes[1].set_title('Shared-support prior')
    save(fig, 'prior-support')
    fig, ax = plt.subplots(figsize=(6.4, 3.3))
    groups = stage1['groups']
    x = np.array([g['native_r_median_span']*100 for g in groups])
    y = np.array([g['common_r_median_span']*100 for g in groups])
    empty = np.array([g['common_mc_intersection_prior_mass'] == 0 for g in groups])
    ax.scatter(x[~empty], y[~empty], color=colors[0], s=22, label='Nonempty shared MC-pass domain')
    ax.scatter(x[empty], y[empty], facecolors='none', edgecolors=colors[1], s=30, label='Empty shared MC-pass domain')
    ax.plot([0, 13], [0, 13], '--', color='gray', linewidth=1)
    ax.set(xlim=(0, 13), ylim=(0, 13), xlabel='Original-prior median span (percentage points)',
           ylabel='Common-prior median span (percentage points)')
    ax.legend(fontsize=8, loc='upper left')
    save(fig, 'common-prior-groups')
    fig, axes = plt.subplots(2, 1, figsize=(7, 5.8))
    modes = ['events', 'injections', 'both']
    for ai, ax in enumerate(axes):
        for pi, pair in enumerate(stage2['pairs'][ai*3:(ai+1)*3]):
            ax.scatter(pair['nominal_difference']*100, pi, marker='D', s=20, color='black', zorder=4)
            for mi, mode in enumerate(modes):
                lo, med, hi = np.array(pair['modes'][mode]['quantiles'])*100
                ax.errorbar(med, pi+(mi-1)*.18, xerr=[[med-lo], [hi-med]], fmt='o', markersize=3,
                            capsize=3, color=colors[mi], label=mode if pi == 0 else None)
        ax.axvline(0, color='gray', linestyle='--', linewidth=1)
        ax.set_yticks([0, 1, 2], ['0.3 minus 0.6', '0.3 minus 1.2', '0.6 minus 1.2'])
        ax.set_xlabel('Difference of posterior r medians (percentage points)')
        ax.set_title(['Low-spin amplitude: Beta(1.5, 8.5)', 'Baseline amplitude: Beta(2, 5)'][ai])
        ax.legend(fontsize=8, loc='best', ncol=3)
    save(fig, 'paired-bootstrap')


if __name__ == '__main__':
    main()
