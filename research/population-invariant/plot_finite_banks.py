"""Visualize conditional finite-bank sensitivity without physical error bars."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parent


def main():
    path=BASE/'paired-pe-results.json';report=json.loads(path.read_text())
    validation=json.loads((BASE/'paired-pe-validation.json').read_text())
    assert validation['result_sha256']==hashlib.sha256(path.read_bytes()).hexdigest()
    rows={r['name']:r for r in report['contrasts']}
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(13,9),layout='constrained')
    cases=[('event_only_5_50_95_pp',-.18,'#57916b','Event thinning'),('selection_only_5_50_95_pp',0,'#277e9e','O3 thinning'),('joint_5_50_95_pp',.18,'#cc6b36','Joint thinning')]
    def ranges(ax,names,labels,title,xlabel):
        y=np.arange(len(names))
        for field,offset,color,label in cases:
            q=np.array([rows[name][field] for name in names])
            ax.errorbar(q[:,1],y+offset,xerr=np.array([q[:,1]-q[:,0],q[:,2]-q[:,1]]),fmt='o',markersize=4,capsize=2,color=color,label=label)
        ax.scatter([rows[name]['full_bank_value_pp'] for name in names],y,c='black',s=25,marker='D',label='Full banks',zorder=4)
        ax.axvline(0,c='.5',lw=.8);ax.set(yticks=y,yticklabels=labels,title=title,xlabel=xlabel)
        ax.invert_yaxis();ax.grid(axis='x',alpha=.15)
    ranges(axes[0,0],['matched:native_width:low_spin','matched:native_width:baseline_spin'],['Beta(1.5,8.5)','Beta(2,5)'],
        'Native-tilt width contrasts:\nthe small ordering is unstable','Median difference in r, sigma 0.3 minus 1.2 (pp)')
    axes[0,0].set_ylim(1.8,-.6);axes[0,0].legend(frameon=False,fontsize=8,loc='lower left',ncol=2)
    ranges(axes[0,1],[f'matched:effective_amplitude:lambda{l}' for l in [0,4,16]],['lambda0','lambda4','lambda16'],
        'Effective-spin amplitude contrasts:\nordering persists in the specified thinnings','Median difference in g, Beta(2,5) minus Beta(1.5,8.5) (pp)')
    axes[0,1].set_xlim(0,18)
    ranges(axes[1,0],[f'effective_cwb_minus_matched:model{k}' for k in range(6)],
        [f'{a}, lambda{l}' for a in ['low spin','baseline spin'] for l in [0,4,16]],
        'Search-policy contrasts:\nmost central ranges cross zero','Including cWB minus matched-filter g median (pp)')
    ax=axes[1,1]
    names=['matched:native_width:low_spin','matched:native_width:baseline_spin','matched:effective_amplitude:lambda0','matched:effective_amplitude:lambda16']
    labels=['Native width, low spin','Native width, baseline spin','Effective amplitude, lambda0','Effective amplitude, lambda16']
    values=[rows[n]['interaction_rms_pp'] for n in names]
    ax.barh(np.arange(4),values,color='#966b9d');ax.set(yticks=range(4),yticklabels=labels,xlabel='RMS joint-minus-event-minus-O3 interaction (pp)',title='Nonlinear propagation:\ninteraction of event and O3 perturbations')
    ax.invert_yaxis();ax.grid(axis='x',alpha=.15)
    fig.suptitle('Fixed 69 events: paired posterior-bank and injection-bank sensitivity',fontsize=15)
    fig.supxlabel('Lines: 5th–95th percentiles of 512 specified artificial thinnings; points: medians. O1/O2 fixed.\nEach event: 6,913 retained rows. O3: 128,454 of 284,672 accepted rows. These are not astrophysical intervals or real-campaign errors.',fontsize=9)
    fig.savefig(BASE/'finite-bank-comparison.svg');fig.savefig(BASE/'finite-bank-comparison.png',dpi=170);plt.close(fig)


if __name__=='__main__':main()
