"""Visual diagnostics for dependent representations of one PE analysis."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parent

def main():
    j=json.loads((BASE/'tail-results.json').read_text())
    fig,axes=plt.subplots(2,2,figsize=(12,9),layout='constrained')
    colors={'original':'#777777','cosmo':'#bc761b','nocosmo':'#236b9a'}
    labels={'original':'Tutorial subset (6,913)','cosmo':'Cosmo release (11,786)','nocosmo':'Original-prior release (169,040)'}
    ax=axes[0,0]
    with np.load(BASE/'tail-checkpoints.npz',allow_pickle=False) as cp:
        for name in colors:
            a=cp[f'sorted_fraction_{name}'];x=np.arange(1,len(a)+1)
            take=np.unique(np.r_[np.arange(min(100,len(a))),np.geomspace(1,len(a),500).astype(int)-1])
            ax.plot(x[take],100*np.cumsum(a)[take],label=labels[name],color=colors[name],lw=1.8)
    ax.set_xscale('log');ax.set_ylim(0,102);ax.axhline(50,color='#bbb',ls=':',lw=1)
    ax.set_xlabel('Number of largest-weight draws (log scale)')
    ax.set_ylabel('Cumulative event-integral contribution (%)')
    ax.set_title('Weight concentration at the original r median')
    ax.legend(fontsize=8,loc='lower right',frameon=False)
    ax=axes[0,1];y=np.arange(6)
    for name,offset,color,label in [('cosmo_lal',-.08,colors['cosmo'],'Cosmo release'),('nocosmo_lal',.08,colors['nocosmo'],'Original-prior release')]:
        ax.scatter(j['scenarios'][name]['median_shift_from_original_pp'],y+offset,color=color,s=35,label=label)
    ax.axvline(0,color='#aaa',lw=1)
    ax.set_yticks(y,[f'Low spin, sigma={s}' for s in [.3,.6,1.2]]+[f'Baseline spin, sigma={s}' for s in [.3,.6,1.2]])
    ax.invert_yaxis();ax.set_xlabel('r median shift from tutorial baseline (percentage points)')
    ax.set_title('Only GW190805 replaced; other 68 events fixed')
    ax.set_ylim(6.1,-.5);ax.legend(fontsize=8,loc='lower right',frameon=False)
    ax=axes[1,0]
    with np.load(BASE/'representations.npz',allow_pickle=False) as d:
        a,b=d['nocosmo_a_1'],d['nocosmo_a_2'];low=(a<.2)&(b<.2)
        kept=np.zeros(len(a),bool);kept[d['cosmo_to_nocosmo_rows']]=True
        ax.scatter(a[low&~kept],b[low&~kept],s=12,alpha=.5,color=colors['nocosmo'],label='Not in cosmo subset')
        ax.scatter(a[low&kept],b[low&kept],s=25,alpha=.8,color=colors['cosmo'],label='Retained in cosmo subset')
        i=j['tracked_nocosmo_row']
        ax.scatter([a[i]],[b[i]],marker='*',s=190,color='#ad3154',edgecolor='black',lw=.4,label='Previously influential draw',zorder=5)
    ax.set_xlim(0,.2);ax.set_ylim(0,.2)
    ax.set_xlabel('Primary spin magnitude a1');ax.set_ylabel('Secondary spin magnitude a2')
    ax.set_title('Prespecified low-spin region: 547 original / 51 cosmo draws')
    ax.legend(fontsize=8,loc='upper left',frameon=True)
    ax=axes[1,1];values=[]
    for label in ['cosmo','nocosmo']:
        keep=j['scenarios'][label+'_lal']['median_shift_from_original_pp'][0]
        drop=j['scenarios'][label+'_drop_tracked']['median_shift_from_original_pp'][0]
        values.append(drop-keep)
    bars=ax.bar(['Cosmo release','Original-prior release'],values,color=[colors['cosmo'],colors['nocosmo']],width=.5)
    for bar,value in zip(bars,values):ax.text(bar.get_x()+bar.get_width()/2,value+.003,f'{value:.5f}',ha='center',fontsize=10)
    ax.set_ylim(0,.15);ax.set_ylabel('Change in r median after removing one draw\n(percentage points)')
    ax.set_title('Same draw, reduced single-point influence')
    ax.text(.5,.93,'Removal is a diagnostic only; all rows retained in main results.',ha='center',va='top',transform=ax.transAxes,fontsize=8)
    fig.suptitle('GW190805: preserve weighted original draws to avoid resampling loss\nSame waveform and likelihood; conditional diagnostics, not independent observations',fontsize=13)
    fig.savefig(BASE/'prior-representation-tail.svg')
    fig.savefig(BASE/'prior-representation-tail.png',dpi=145)
    plt.close(fig)

if __name__=='__main__':main()
