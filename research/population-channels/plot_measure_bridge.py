"""Scientific interface diagnostics; no observed channel posterior is plotted."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from fast_measure import effective_spin_density

BASE=Path(__file__).resolve().parent


def main():
    norm=json.loads((BASE/'channel-log-normalization.json').read_text())
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(11.8,4.6))
    x=np.linspace(-.989,.989,1201)
    for family,color,label in [('uniform_magnitude','#285ba6','Uniform magnitude'),('beta2_5_magnitude','#c05a21','Beta(2,5) magnitude')]:
        for q,style in [(.1,'-'),(1.,'--')]:
            axes[0].plot(x,effective_spin_density(x,q,.99,family),style,color=color,lw=1.7,label=f'{label}, q={q:g}')
    axes[0].set(xlabel=r'Effective spin $\chi_{\rm eff}$',ylabel=r'Conditional density $f_Q(\chi_{\rm eff}\mid q)$',title='A  Two explicit hidden-spin assumptions')
    axes[0].legend(frameon=False,fontsize=8.5);axes[0].grid(alpha=.16)
    values=[r['fraction'] for r in norm['channels']];errors=[1.96*r['bernoulli_mc_standard_error'] for r in norm['channels']]
    positions=np.arange(5)
    axes[1].errorbar(positions,values,yerr=errors,fmt='o',color='#285ba6',capsize=4,ms=6)
    axes[1].set(xticks=positions,xticklabels=[r['channel'] for r in norm['channels']],xlim=(-.5,4.5),ylim=(0,1),ylabel=r'Emulator mass inside declared domain $Z_c$',title='B  Restricted-domain channel normalization')
    axes[1].grid(axis='y',alpha=.16)
    for i,v in enumerate(values):axes[1].text(i,v+.045,f'{v:.4f}',ha='center',fontsize=9)
    fig.suptitle('Physical-channel interface checks: fixed birth spin 0.1, CE efficiency 2',fontsize=12,y=.97)
    fig.text(.5,.025,'No observed branching fit. Panel A shows auxiliary priors, not channel spin distributions.\nPanel B: 262,144 emulator draws/channel; bars are 1.96 MC standard errors.\nDomain: source m1 5-100, m2 >=3 solar masses; q <=1, z <=1.9, |effective spin| <0.99.',ha='center',va='bottom',fontsize=8.5)
    fig.subplots_adjust(left=.07,right=.98,top=.83,bottom=.24,wspace=.3)
    fig.savefig(BASE/'measure-bridge-diagnostics.svg')
    fig.savefig(BASE/'measure-bridge-diagnostics.png',dpi=160)
    plt.close(fig)


if __name__=='__main__':main()
