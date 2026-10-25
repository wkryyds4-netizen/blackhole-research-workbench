"""Plot common-random-row policy contrasts, with conditional scope visible."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parent


def main():
    result=json.loads((BASE/'paired-policy-results.json').read_text())
    with np.load(BASE/'paired-policy-checkpoints.npz',allow_pickle=False) as bank:gaps=bank['gaps']
    with np.load(BASE/'policy-influence-checkpoints.npz',allow_pickle=False) as bank:
        presence=bank['presence'];top_rows=bank['top_rows']
    full=np.array(result['full_gaps_pp'])
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(11.5,8),layout='constrained')
    for k,label in enumerate(['Low spin: Beta(1.5, 8.5)','Baseline spin: Beta(2, 5)']):
        ax=axes[k,0];x=gaps[:,0,k];y=gaps[:,1,k]
        ax.scatter(x,y,s=13,c='#176b87',alpha=.35,label='512 shared-row thinnings')
        lim=[min(x.min(),y.min())-.15,max(x.max(),y.max())+.15]
        ax.plot(lim,lim,color='0.4',ls='--',lw=1,label='Equal gap')
        ax.scatter(full[0,k],full[1,k],marker='*',s=180,c='#c65035',edgecolor='white',label='Full bank',zorder=3)
        ax.set(xlim=lim,ylim=lim,xlabel='Width-model gap, union with cWB (pp)',ylabel='Width-model gap, matched filters (pp)',title=label)
        ax.text(.04,.95,f'Paired correlation = {result["diagnostics"][k]["union_matched_gap_correlation"]:.3f}',transform=ax.transAxes,va='top')
        ax.grid(alpha=.15)
        if k==0:ax.legend(frameon=False,fontsize=8,loc='lower right')
        ax=axes[k,1];delta=y-x
        yes=presence[:,k]
        ax.hist([delta[~yes],delta[yes]],bins=30,stacked=True,color=['#176b87','#d49c42'],alpha=.8,edgecolor='white',linewidth=.4,
            label=[f'Row {top_rows[k]} excluded',f'Row {top_rows[k]} included'])
        ax.axvline(full[1,k]-full[0,k],color='#c65035',lw=2,label='Full-bank policy change')
        ax.axvline(0,color='0.4',ls='--',lw=1)
        d=result['diagnostics'][k]
        ax.set(xlabel='Matched-filter minus union gap (pp)',ylabel='Artificial thinnings',title='Search-policy contrast using the same rows')
        ax.text(.04,.95,f'Paired SD: {d["policy_change_std_pp"]:.3f} pp\nIgnoring covariance: {d["unpaired_variance_sum_std_pp"]:.3f} pp',transform=ax.transAxes,va='top',fontsize=9)
        ax.legend(frameon=False,fontsize=8,loc='upper right')
    fig.suptitle('Search policy and finite-bank sampling must be compared jointly',fontsize=15)
    fig.supxlabel('Conditional artificial SRS without replacement; not astrophysical intervals or historical-campaign coverage',fontsize=9)
    fig.savefig(BASE/'paired-policy-thinning.svg');fig.savefig(BASE/'paired-policy-thinning.png',dpi=170)
    plt.close(fig)


if __name__=='__main__':main()
