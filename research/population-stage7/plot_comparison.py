"""Four-event source representation comparison and joint conditional shifts."""
from pathlib import Path
import argparse
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parent


def main():
    parser=argparse.ArgumentParser();parser.add_argument('--prefix',default='');args=parser.parse_args()
    result=json.loads((BASE/(args.prefix+'comparison-results.json')).read_text())
    fig,axes=plt.subplots(2,2,figsize=(12,8.5),layout='constrained')
    names=[e['event'] for e in result['diagnostics']];x=np.arange(len(names));width=.22
    colors=['#6e7981','#d58b26','#277b8e']
    for offset,(label,title) in enumerate([('original','Tutorial'),('cosmo','Full cosmo'),('nocosmo','Full original-prior')]):
        rows=[e['representations'][label][0] for e in result['diagnostics']]
        axes[0,0].bar(x+(offset-1)*width,[e['ess'] for e in rows],width,color=colors[offset],label=title)
        axes[0,1].bar(x+(offset-1)*width,[100*e['maximum_contribution'] for e in rows],width,color=colors[offset],label=title)
    for ax in axes[0]:ax.set_xticks(x,[n.split('_')[0] for n in names]);ax.grid(axis='y',alpha=.2)
    axes[0,0].set_yscale('log');axes[0,0].set_ylabel('Population-integral contribution ESS')
    axes[0,0].set_title('(a) Low-spin sigma=0.3, original median');axes[0,0].legend(fontsize=8)
    axes[0,1].set_ylabel('Largest single-row contribution (%)');axes[0,1].set_title('(b) Same model and evaluation point')
    k=np.arange(6)
    for label,title,color in [('all_new_cosmo','New events: cosmo','#d58b26'),('all_new_nocosmo','New events: original-prior','#277b8e'),('new_plus_38_nocosmo','New events + GW190805: original-prior','#b5526d')]:
        y=result['scenarios'][label]['median_shift_pp']
        axes[1,0].scatter(k,y,label=title,color=color,s=35)
    axes[1,0].axhline(0,color='black',lw=.7);axes[1,0].set_ylabel('Conditional median shift from tutorial (pp)')
    axes[1,0].set_xticks(k,['low\n0.3','low\n0.6','low\n1.2','base\n0.3','base\n0.6','base\n1.2'])
    axes[1,0].set_title('(c) Joint replacements; selection unchanged');axes[1,0].legend(fontsize=8);axes[1,0].grid(alpha=.2)
    scopes=['original','all_new_cosmo','all_new_nocosmo','new_plus_38_nocosmo'];y=np.arange(4)
    for j,label in enumerate(['Low-spin magnitudes','Baseline-spin magnitudes']):
        values=[result['scenarios'][s]['sigma03_minus12_gap_pp'][j] for s in scopes]
        axes[1,1].plot(values,y,marker='o',label=label)
    axes[1,1].set_yticks(y,['Tutorial baseline','New events: cosmo','New events: original-prior','New + GW190805: original-prior'])
    axes[1,1].invert_yaxis();axes[1,1].set_xlabel('Median difference: sigma=0.3 minus 1.2 (pp)')
    axes[1,1].set_title('(d) Conditional model-width contrast');axes[1,1].legend(fontsize=8);axes[1,1].grid(alpha=.2)
    fig.suptitle('Targeted multi-event posterior-representation audit',fontsize=15)
    fig.supxlabel('Same 69-event conditional model; native spin references retained.\nESS is not a convergence proof; finite-array contrasts are not calibrated uncertainty intervals.',fontsize=9)
    fig.savefig(BASE/(args.prefix+'representation-comparison.svg'))
    fig.savefig(BASE/(args.prefix+'representation-comparison.png'),dpi=150)
    plt.close(fig)


if __name__=='__main__':main()
