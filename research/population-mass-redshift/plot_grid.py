"""Vector scientific summary; model envelopes are explicitly not intervals."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parent


def main():
    d=json.loads((BASE/'results.json').read_text());v=json.loads((BASE/'validation.json').read_text())
    assert v['status']=='pass' and v['result_sha256']==hashlib.sha256((BASE/'results.json').read_bytes()).hexdigest()
    name='official_mixture:matched_filter';case=d['cases'][name]
    q=np.array(case['quantiles']).reshape((3,3,2,3,3))
    plt.rcParams.update({'font.size':9,'axes.spines.top':False,'axes.spines.right':False,
                         'pdf.fonttype':42,'svg.fonttype':'none'})
    fig,axes=plt.subplots(2,2,figsize=(11.5,7.9),layout='constrained')
    ax=axes[0,0]
    values=(100*(q[:,:,1,:,1]-q[:,:,0,:,1])).reshape((9,3)).T
    im=ax.imshow(values,aspect='auto',vmin=10,vmax=17.5,cmap='viridis',origin='lower')
    for row in range(3):
        for col in range(9):ax.text(col,row,f'{values[row,col]:.2f}',ha='center',va='center',color='white' if values[row,col]<14.6 else '#14212b',fontsize=8)
    ax.set(xticks=range(9),xticklabels=[f'M{mi}\n{k}' for mi in range(3) for k in [0,2,4]],
           yticks=range(3),yticklabels=['0','4','16'],xlabel=r'Mass model / redshift exponent $\kappa$',
           ylabel=r'Even-spin shape $\lambda$',title='(a) Same-shape amplitude contrasts (pp)')
    fig.colorbar(im,ax=ax,shrink=.85,pad=.025)
    ax=axes[0,1]
    y=0
    for si,mag in enumerate(['Beta(1.5,8.5)','Beta(2,5)']):
        for li,shape in enumerate([0,4,16]):
            centers=100*q[:,:,si,li,1]
            ax.plot([centers.min(),centers.max()],[y,y],color='#007f82',lw=3,solid_capstyle='round')
            ax.scatter([centers[0,1]],[y],marker='o',color='#b16a20',s=22,zorder=3)
            y+=1
    ax.set(yticks=range(6),yticklabels=[f'{mag}, $\\lambda={shape}$' for mag in ['Beta(1.5,8.5)','Beta(2,5)'] for shape in [0,4,16]],
           xlabel='g posterior median (%)',title='(b) Median ranges over nine model choices',xlim=(10,34))
    ax.invert_yaxis();ax.grid(axis='x',alpha=.2)
    ax.set_ylim(6.8,-.9)
    ax.plot([],[],color='#007f82',lw=3,label='Mass/redshift model range')
    ax.scatter([],[],color='#b16a20',s=22,label=r'M0, $\kappa=2$ reference')
    ax.legend(frameon=False,fontsize=7,loc='lower right')
    ax.text(.02,.98,'Model range; not a credible interval',transform=ax.transAxes,va='top',fontsize=8)
    ax=axes[1,0]
    s=case['median_sensitivity_pp']
    amounts=[s['mass_only_maximum_span'],s['redshift_only_maximum_span'],s['joint_mass_redshift_maximum_span'],d['bank_comparisons']['matched_filter']['maximum_median_shift_pp']]
    ax.barh(range(4),amounts,color=['#6c9fab','#6c9fab','#007f82','#b16a20'],height=.58)
    for i,value in enumerate(amounts):ax.text(value+.04,i,f'{value:.3f}',va='center',fontsize=9)
    ax.set(yticks=range(4),yticklabels=['Mass prescription','Redshift exponent','Joint model span','O3 injection bank'],
           xlabel='Largest median change (pp)',title='(c) Sizes of declared model perturbations',xlim=(0,max(amounts)*1.24))
    ax.invert_yaxis();ax.grid(axis='x',alpha=.2)
    ax.set_ylim(4.0,-.55)
    ax.text(.99,.02,'Grid maxima; not uncertainty budgets',ha='right',transform=ax.transAxes,fontsize=8)
    ax=axes[1,1]
    support=np.array([d['support']['models'][mi*18+zi*6]['banks']['official_mixture']['uncovered_target_probability']*100 for mi in range(3) for zi in range(3)]).reshape((3,3))
    im=ax.imshow(support,vmin=0,vmax=50,cmap='Oranges',aspect='auto',origin='lower')
    for row in range(3):
        for col in range(3):ax.text(col,row,f'{support[row,col]:.3f}%',ha='center',va='center',color='white' if support[row,col]>32 else '#14212b',fontsize=10)
    ax.set(xticks=range(3),xticklabels=['0','2','4'],yticks=range(3),yticklabels=['M0','M1','M2'],
           xlabel=r'Redshift exponent $\kappa$',ylabel='Mass prescription',title='(d) Target probability outside O3 proposal')
    fig.colorbar(im,ax=ax,shrink=.85,pad=.025,label='Population probability (%)')
    fig.suptitle('Full69-source spin controls: 54 mass/redshift models\nOfficial O3 mixture / matched-filter union; fixed Uniform[0,1] prior on g',fontsize=12)
    fig.supxlabel(r'M0: $\alpha=2$ + 10% peak; M1: $\alpha=3$ + 10% peak; M2: $\alpha=2$, no peak. Missing population mass is not detection error.',fontsize=8)
    for suffix in ['svg','pdf','png']:fig.savefig(BASE/f'mass-redshift-comparison.{suffix}',dpi=170)
    plt.close(fig)
    print('Saved mass-redshift-comparison.svg/pdf and local PNG preview.')


if __name__=='__main__':main()
