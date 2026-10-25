"""Scientific figure from saved summaries; no raw-data recomputation."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parent


def main():
    study=json.loads((BASE/'reference-results.json').read_text())
    init=json.loads((BASE/'initialization-results.json').read_text())
    audit=json.loads((BASE/'reference-audit.json').read_text())
    primary=json.loads((BASE/'reference-validation.json').read_text())
    plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'svg.fonttype':'none'})
    fig,axes=plt.subplots(2,2,figsize=(12.5,8.5),layout='constrained')
    ax=axes[0,0];proof=study['counterexample'];phi=np.array([r['phi12'] for r in proof])/np.pi
    for j,color in [(1,'#176b87'),(2,'#b85042')]:
        ax.plot(phi,[r[f'cos_tilt_{j}_infinity'] for r in proof],'.-',color=color,label=f'BH {j}')
    ax.axhline(0,color='0.65',lw=.8);ax.set(xlabel=r'Relative spin azimuth $\phi_{12}/\pi$',ylabel=r'$\cos\theta_{i,\infty}$',
        title='A  Missing azimuth prevents unique transport',ylim=(-1.05,1.05));ax.legend(frameon=False)
    ax.text(.03,.04,'Same masses, spin magnitudes and native tilts\nDifferent azimuth; original precession initialization',transform=ax.transAxes,fontsize=8,bbox=dict(facecolor='white',edgecolor='none',alpha=.85))

    ax=axes[0,1];labels=['Original r(f), Newtonian L','Newtonian L(f)','1PN L, no spin terms','2PN L, no spin terms','2PN L, spin magnitude','2PN L, spin vector']
    stats=[study['validation']['released_implementation_abs_difference_quantiles']]+[r['abs_difference_quantiles'] for r in init['summaries']]
    for n,row in enumerate(stats):
        ax.plot([row[0],row[3]],[n,n],color='#176b87',lw=2)
        ax.scatter(row[0],n,c='#176b87',s=25);ax.scatter(row[3],n,c='#b85042',s=25,marker='x')
    ax.set_xscale('log');ax.set_yticks(range(6),labels);ax.invert_yaxis()
    ax.set(xlabel='Absolute cosine difference: median (dot) to max (cross)',title='B  Match the initial angular momentum (80 rows)')
    ax.grid(axis='x',alpha=.2)

    ax=axes[1,0];x=np.arange(5);width=.35
    for j,color in [(1,'#176b87'),(2,'#b85042')]:
        values=[100*e['spins'][j-1]['sign_change_fraction'] for e in study['events']]
        ax.bar(x+(j-1.5)*width,values,width,color=color,label=f'BH {j}')
    ax.set_xticks(x,[e['event'].split('_')[0] for e in study['events']],rotation=15)
    ax.set(ylabel='Posterior samples changing tilt-cosine sign (%)',title='C  Native vs released infinity: 903,026 rows');ax.legend(frameon=False)
    ax.text(.03,.95,'Unweighted event posteriors; not population fractions',transform=ax.transAxes,va='top',fontsize=8)
    ax.set_ylim(0,18)

    ax=axes[1,1];q=np.array([r['inputs']['q'] for r in primary['records']]);err=np.array([r['tolerance_difference'] for r in primary['records']])
    ax.scatter(1-q,np.maximum(err,1e-16),s=16,alpha=.6,c='#176b87',label='Initial tolerance check')
    lookup={(r['index'],r['source_row']):1-r['inputs']['q'] for r in primary['records']}
    refined=[r for r in audit['rows'] if r['refined']]
    ax.scatter([lookup[r['index'],r['source_row']] for r in refined],[r['tolerance_difference'] for r in refined],marker='x',s=50,c='#b85042',label='5 refined checks')
    ax.axhline(1e-6,c='0.45',ls='--',lw=1);ax.set_xscale('log');ax.set_yscale('log')
    ax.set(xlabel=r'Mass asymmetry $1-q$',ylabel='Cosine change between integration tolerances',title='D  Near-equal masses require tighter integration')
    ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.15)
    fig.suptitle('Spin-reference audit: conventions, missing information and numerical limits',fontsize=15)
    fig.savefig(BASE/'spin-reference-audit.svg');fig.savefig(BASE/'spin-reference-audit.png',dpi=170)
    plt.close(fig)


if __name__=='__main__':main()
