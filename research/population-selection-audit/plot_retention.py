"""Show evidence for an inferred retention model and its limited consequences."""
from pathlib import Path
import json
import numpy as np
from scipy.stats import ncx2
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parent


def main():
    result=json.loads((BASE/'cumulative-retention-results.json').read_text())
    with np.load(BASE/'cumulative-bank.npz',allow_pickle=False) as b:
        e=b['early'];obs=b['semianalytic_observed_phase_maximized_snr_net'][e];opt=b['estimated_optimal_snr_net'][e]
        dof=2*sum((b['estimated_optimal_snr_'+i][e]>0).astype(int) for i in 'HLV')
    raw=ncx2.cdf(obs**2,dof,opt**2);conditional=1-ncx2.sf(obs**2,dof,opt**2)/ncx2.sf(8.5**2,dof,opt**2)
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(11.8,8),layout='constrained')
    ax=axes[0,0];bins=np.linspace(7.5,14,66)
    ax.hist(obs,bins=bins,color='#277e9e',edgecolor='white',linewidth=.4)
    ax.axvline(8.5,color='#cc6b36',ls='--',label='Hypothesized retention at 8.5')
    ax.set(xlim=(7.5,14),xlabel='Saved O1/O2 observed network SNR',ylabel='Rows per bin',title='A sharp empirical floor in the released rows')
    ax.text(.97,.95,f'Minimum = {obs.min():.9f}\n593,005 early rows',ha='right',va='top',transform=ax.transAxes);ax.legend(frameon=False,fontsize=9)
    ax=axes[0,1];grid=np.linspace(0,1,1001)
    for x,c,label in [(raw,'#cc6b36','Unconditioned Gaussian model'),(conditional,'#277e9e','Conditioned on observed SNR > 8.5')]:
        ax.plot(grid[1:],np.cumsum(np.histogram(x,grid)[0])/len(x)-grid[1:],color=c,label=label)
    ax.axhline(0,color='0.4',ls='--');ax.set(xlabel='Probability integral transform u',ylabel='Empirical CDF(u) minus u',title='Retention explains most of the distribution mismatch')
    ax.legend(frameon=False,fontsize=8,loc='center',bbox_to_anchor=(.5,.68))
    labels=[f'{g}, width {s}' for g in ['Low','Baseline'] for s in ['0.3','0.6','1.2']];y=np.arange(6)
    ax=axes[1,0];delta=np.array(result['corrected_optimal_minus_retained_optimal_median_pp'])
    ax.barh(y,delta,color=['#277e9e']*3+['#cc6b36']*3);ax.axvline(0,color='0.4',lw=1)
    ax.set(yticks=y,yticklabels=labels,xlabel='Corrected minus retained-optimal median r (pp)',title='Inverse-retention correction under the hypothesis');ax.invert_yaxis();ax.ticklabel_format(axis='x',style='sci',scilimits=(0,0))
    ax=axes[1,1];q=np.array(result['conditional_noise_simulation']['median_shift_5_50_95_pp']);actual=result['released_observed_minus_expected_median_pp']
    ax.hlines(y,q[:,0],q[:,2],color='#277e9e',lw=3,label='5%-95% of conditional simulations')
    ax.scatter(q[:,1],y,c='#277e9e',s=25,label='Simulation median')
    ax.scatter(actual,y,c='#cc6b36',marker='x',s=50,label='Released noise realization')
    ax.axvline(0,color='0.4',ls='--');ax.set(yticks=y,yticklabels=labels,xlabel='Median r minus conditional-expectation result (pp)',title='Conditional noise redraws (256 replicates)',ylim=(7.2,-.5))
    ax.legend(frameon=False,fontsize=8,loc='lower left')
    fig.suptitle('Observed-SNR retention: empirical evidence and a conditional correction',fontsize=14)
    fig.supxlabel('Generation configuration not recovered. Gaussian noise + retention at 8.5 is a hypothesis; these are not astrophysical intervals.',fontsize=9)
    fig.savefig(BASE/'cumulative-retention.svg');fig.savefig(BASE/'cumulative-retention.png',dpi=170);plt.close(fig)


if __name__=='__main__':main()
