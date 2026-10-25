"""Scientific figure for full-O3 selection and representation controls."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parent


def main():
    d=json.loads((BASE/'selection-results.json').read_text())
    fig,axes=plt.subplots(2,2,figsize=(12,8.5),layout='constrained');x=np.arange(6)
    labels=['low\n0.3','low\n0.6','low\n1.2','base\n0.3','base\n0.6','base\n1.2']
    for name,title,color in [('five_frozen','Five event replacements','#b5526d'),('five_full_o3','Also use full O3 injections','#277b8e')]:
        axes[0,0].scatter(x,d['scenarios'][name]['median_shift_from_frozen_pp'],label=title,color=color)
    axes[0,0].set_title('(a) Joint shifts from frozen tutorial');axes[0,0].set_ylabel('Conditional r median shift (pp)');axes[0,0].legend(fontsize=8)
    axes[0,1].scatter(x,d['five_event_full_o3_shift_pp'],color='#277b8e')
    axes[0,1].set_title('(b) Full O3 change alone, with five replacements');axes[0,1].set_ylabel('Conditional r median shift (pp)')
    for ax in axes[0]:ax.set_xticks(x,labels);ax.axhline(0,color='black',lw=.7);ax.grid(alpha=.2)
    scopes=['tutorial_frozen','five_frozen','five_full_o3','five_full_o3_zcap'];y=np.arange(4)
    for i,title in enumerate(['Low-spin magnitudes','Baseline-spin magnitudes']):
        axes[1,0].plot([d['scenarios'][s]['sigma03_minus12_gap_pp'][i] for s in scopes],y,'o-',label=title)
    axes[1,0].set_yticks(y,['Frozen tutorial','Five event replacements','Also full O3','Also joint z <= 1.9 cap'])
    axes[1,0].invert_yaxis();axes[1,0].set_xlabel('Median difference: sigma=0.3 minus 1.2 (pp)')
    axes[1,0].set_title('(c) Model-width contrast');axes[1,0].legend(fontsize=8);axes[1,0].grid(alpha=.2)
    for off,key,title,color in [(-.18,'o3_combined_ess','Combined-file O3','#6e7981'),(.18,'o3_full_ess','Full O3','#277b8e')]:
        axes[1,1].bar(x+off,[v[key] for v in d['selection_diagnostics']],width=.36,label=title,color=color)
    axes[1,1].set_yscale('log');axes[1,1].set_xticks(x,labels);axes[1,1].set_title('(d) O3 contribution ESS at prior-stage medians')
    axes[1,1].set_ylabel('Contribution ESS');axes[1,1].legend(fontsize=8);axes[1,1].grid(axis='y',alpha=.2)
    fig.suptitle('Selection representation matters for the small model contrast',fontsize=15)
    fig.supxlabel('Same 69 events, six fixed models, native spin references. Related finite injection arrays.\nA small z-cap shift does not bound missing high-redshift detection efficiency.',fontsize=9)
    for ext in ['svg','png']:fig.savefig(BASE/f'selection-comparison.{ext}',dpi=150)
    plt.close(fig)
    thinning=BASE/'thinning-results.json'
    if thinning.exists():
        t=json.loads(thinning.read_text())
        with np.load(BASE/'thinning-checkpoints.npz',allow_pickle=False) as cp:q=cp['quantiles']
        gaps=100*np.stack([q[:,0,1]-q[:,2,1],q[:,3,1]-q[:,5,1]],axis=1)
        fig,axes=plt.subplots(1,2,figsize=(11,4.5),layout='constrained')
        for k,title in enumerate(['Low-spin magnitudes','Baseline-spin magnitudes']):
            ax=axes[k];ax.hist(gaps[:,k],bins=32,color='#a4ccd0',edgecolor='white')
            ax.axvline(t['full_bank_gaps_pp'][k],color='#277b8e',lw=2,label='Full saved bank')
            ax.axvline(t['actual_combined_subset_gaps_pp'][k],color='#b5526d',lw=2,label='Actual combined subset')
            ax.axvline(0,color='black',lw=.7,ls=':');ax.set_title(title);ax.grid(axis='y',alpha=.2)
            ax.set_xlabel('Median difference: sigma=0.3 minus 1.2 (pp)');ax.set_ylabel('Artificial thinning replicates')
            ax.legend(fontsize=8)
        fig.suptitle('512 paired random thinnings of the fixed O3 bank',fontsize=14)
        fig.supxlabel('Uniform sampling without replacement; 128,454 of 284,672 accepted injections.\nThese distributions are conditional numerical benchmarks, not astrophysical uncertainty intervals.',fontsize=9)
        for ext in ['svg','png']:fig.savefig(BASE/f'thinning-comparison.{ext}',dpi=150)
        plt.close(fig)


if __name__=='__main__':main()
