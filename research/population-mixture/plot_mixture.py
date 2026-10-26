"""Scientific summary of fixed-input O3 mixture comparisons; no error bars."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parent


def main():
    d=json.loads((BASE/'results.json').read_text())
    v=json.loads((BASE/'validation.json').read_text())
    assert v['status']=='pass' and v['result_sha256']==hashlib.sha256((BASE/'results.json').read_bytes()).hexdigest()
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,
                         'svg.fonttype':'none','pdf.fonttype':42})
    fig,axes=plt.subplots(2,2,figsize=(10,7.2),layout='constrained')
    colors={'matched_filter':'#007f82','with_cwb':'#b16a20'}
    labels={'matched_filter':'Matched-filter union','with_cwb':'Union with cWB'}
    x=np.arange(3)
    ax=axes[0,0]
    for policy,color in colors.items():
        cases=d['metrics']['g']['scenarios'][policy]
        for case,ls,marker in [('pure_bbh_baseline','--','o'),('mixture_full_union','-','s')]:
            y=cases[case]['model_median_contrasts_pp']
            ax.plot(x,y,ls=ls,marker=marker,color=color,label=labels[policy]+' / '+('BBH' if case=='pure_bbh_baseline' else 'mixture'),ms=4)
    ax.set(xticks=x,xticklabels=['0','4','16'],xlabel=r'Even-shape parameter $\lambda$',
           ylabel='Same-shape amplitude contrast (pp)',title='(a) Negative-effective-spin fraction g')
    ax.legend(fontsize=7,loc='upper left',frameon=False)
    ax.set_ylim(11.5,16.2);ax.grid(alpha=.2)
    ax=axes[0,1]
    for i,(policy,color) in enumerate(colors.items()):
        cases=d['metrics']['r']['scenarios'][policy]
        y=[cases[case]['model_median_contrasts_pp'][0] for case in ['pure_bbh_baseline','mixture_full_union']]
        ax.plot([0,1],y,marker='o',color=color,label=labels[policy])
        for pos,value in enumerate(y):ax.annotate(f'{value:+.4f}',(pos,value),xytext=(5,6 if i else -13),textcoords='offset points',color=color,fontsize=8)
    ax.axhline(0,color='.45',lw=.8)
    ax.set(xticks=[0,1],xticklabels=['Pure BBH bank','Four-population mixture'],xlim=(-.25,1.35),ylim=(-.11,.11),
           ylabel='Native-tilt width contrast (pp)',title='(b) Low-spin r diagnostic: tiny ordering')
    ax.text(.02,.96,r'Beta(1.5,8.5): $\sigma=0.3$ minus $\sigma=1.2$',transform=ax.transAxes,va='top',fontsize=8)
    ax.legend(fontsize=7,loc='lower right',frameon=False);ax.grid(alpha=.2)
    ax=axes[1,0]
    names=[f'Beta({aa},{bb}), $\\lambda={lam}$' for aa,bb in [('1.5','8.5'),('2','5')] for lam in [0,4,16]]
    for policy,color in colors.items():
        shift=d['metrics']['g']['scenarios'][policy]['mixture_full_union']['median_shift_from_pure_bbh_pp']
        ax.plot(shift,np.arange(6),marker='o',color=color,label=labels[policy],ls='none')
    ax.set(yticks=np.arange(6),yticklabels=names,xlabel='Shift in g posterior median (pp)',
           title='(c) Change of saved injection representation',xlim=(-.01,.16))
    ax.axvline(0,color='.5',lw=.8);ax.invert_yaxis();ax.grid(axis='x',alpha=.2)
    ax.legend(fontsize=7,loc='lower right',frameon=False)
    ax=axes[1,1]
    s=d['support']['magnitude_families'][0]
    values=[100*s['old_bbh_covered_target_mass'],100*s['mixture_added_target_mass'],100*s['still_uncovered_target_mass']]
    left=0
    for value,color,label in zip(values,['#6da8b0','#e9b15f','#d8dce0'],['Previously covered','Added by IMBH proposal','Still outside union']):
        ax.barh([0],value,left=left,color=color,height=.45,label=label);left+=value
    ax.set(xlim=(0,100),ylim=(-.65,.8),yticks=[],xlabel='Target probability mass (%)',
           title='(d) Proposal support, before hopeless cuts')
    ax.text(.02,.94,f"Added: {values[1]:.4f} percentage points\nStill uncovered: {values[2]:.4f}%",transform=ax.transAxes,va='top',fontsize=9)
    ax.text(.02,.06,'178 added-region rows; zero saved detections.\nNeither count bounds actual omitted sensitivity.',transform=ax.transAxes,fontsize=8,va='bottom')
    ax.legend(fontsize=7,loc='lower center',bbox_to_anchor=(.5,-.36),frameon=False)
    fig.suptitle('Official O3 mixture: conditional all-69-source comparison',fontsize=13)
    for suffix in ['svg','pdf','png']:
        fig.savefig(BASE/f'mixture-comparison.{suffix}',dpi=170)
    plt.close(fig)
    print('Saved mixture-comparison.svg/pdf and local PNG preview.')


if __name__=='__main__':main()
