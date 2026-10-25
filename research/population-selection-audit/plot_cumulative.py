"""Distinguish bundled release changes from a paired SNR-threshold control."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parent


def main():
    result=json.loads((BASE/'cumulative-results.json').read_text())
    with np.load(BASE/'cumulative-checkpoints.npz',allow_pickle=False) as cp:r=cp['r'];fs=cp['fs'];selection=cp['selections']
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(11.8,8),layout='constrained')
    names=list(result['scenarios']);colors=['#727e87','#277e9e','#cc6b36']
    labels=['Historical + full O3\noptimal SNR','New bank retained rows\noptimal SNR mask','Same new bank\nobserved SNR']
    for k,title in enumerate(['Low spin: Beta(1.5, 8.5)','Baseline spin: Beta(2, 5)']):
        ax=axes[0,k];v=[result['scenarios'][n]['width_gaps_pp'][k] for n in names]
        ax.plot([0,1],v[:2],color='0.6',ls=':',lw=2,label='Bundled representation change')
        ax.plot([1,2],v[1:],color='#277e9e',lw=2,label='Retained-row mask change only')
        for x,y,c in zip(range(3),v,colors):
            ax.scatter(x,y,s=65,c=c,zorder=3);ax.annotate(f'{y:.5f}',(x,y),xytext=(0,12),textcoords='offset points',ha='center')
        span=max(v)-min(v)
        ax.set(xlim=(-.4,2.4),ylim=(min(v)-span*.35,max(v)+span*.45),xticks=[0,1,2],xticklabels=labels,
            ylabel=r'Median $r$ gap: $\sigma=0.3$ minus $1.2$ (pp)',title=title)
        ax.grid(axis='y',alpha=.15)
        if k==0:ax.legend(loc='lower right',frameon=False,fontsize=8)
    ax=axes[1,0];delta=result['paired_observed_minus_optimal_median_pp']
    y=np.arange(6);labels=[f'{group}, width {s}' for group in ['Low','Baseline'] for s in ['0.3','0.6','1.2']]
    ax.axvline(0,color='0.4',lw=1)
    ax.barh(y,delta,color=['#277e9e']*3+['#cc6b36']*3,height=.58)
    for yi,v in enumerate(delta):ax.text(v+(.0005 if v>=0 else -.0005),yi,f'{v:+.4f}',ha='left' if v>=0 else 'right',va='center',fontsize=9)
    ax.set(yticks=y,yticklabels=labels,xlim=(-.026,.008),xlabel='Observed minus retained-optimal median r (pp)',title='Mask change; missing rows are not restored');ax.invert_yaxis()
    ax=axes[1,1]
    for k in range(6):
        opt=fs[k]*selection[1,k,0]+(1-fs[k])*selection[1,k,1]
        obs=fs[k]*selection[2,k,0]+(1-fs[k])*selection[2,k,1]
        loglr=-69*np.log(obs/opt);loglr-=loglr[-1]
        ax.plot(r,loglr,color=['#277e9e','#cc6b36'][k//3],ls=['-', '--', ':'][k%3],label=labels[k])
    ax.axhline(0,color='0.4',lw=.7);ax.set(xlabel=r'$r=P(\cos\theta_1<0)$',ylabel='Change in log likelihood relative to r = 0.5',title='Most of the efficiency change cancels in shape')
    ax.legend(frameon=False,fontsize=8,ncol=2);ax.grid(alpha=.15)
    fig.suptitle('O1/O2 threshold comparisons must account for which rows were retained',fontsize=14)
    fig.supxlabel('All saved early rows have observed SNR > 8.500004. Optimal-mask values are not a complete counterfactual. Fixed 69 events; no O4.',fontsize=9)
    fig.savefig(BASE/'cumulative-control.svg');fig.savefig(BASE/'cumulative-control.png',dpi=170);plt.close(fig)


if __name__=='__main__':main()
