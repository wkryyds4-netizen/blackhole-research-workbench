"""Scientific figure for a completed and independently validated control."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parent


def main():
    path=BASE/'control-results.json';report=json.loads(path.read_text())
    validation=json.loads((BASE/'control-validation.json').read_text())
    assert validation['result_sha256']==hashlib.sha256(path.read_bytes()).hexdigest()
    sources=json.loads((BASE/'event-integrals.json').read_text())['events']
    with np.load(BASE/'control-checkpoints.npz',allow_pickle=False) as cp:sel=cp['selection_matched']
    labels=[f'Beta({a},{b}), lambda={lam}' for a,b in [(1.5,8.5),(2,5)] for lam in [0,4,16]]
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12.7,8.5),layout='constrained');y=np.arange(6)
    ax=axes[0,0]
    for name,dy,color,label in [('matched_filter',-.16,'#277e9e','Matched-filter selection'),('with_cwb',0,'#cc6b36','Selection including cWB'),('without_selection_diagnostic',.16,'#92908c','Uncorrected diagnostic')]:
        q=np.asarray(report['scenarios'][name]['quantiles'])*100
        ax.errorbar(q[:,1],y+dy,xerr=np.array([q[:,1]-q[:,0],q[:,2]-q[:,1]]),fmt='o',markersize=3.5,capsize=2,color=color,label=label)
    ax.set(yticks=y,yticklabels=labels,xlabel='Binary population fraction g = P(chi_eff < 0), percent',xlim=(0,100),title='Same 69 full sources; conditional 90% intervals')
    ax.invert_yaxis();ax.grid(axis='x',alpha=.15);ax.legend(frameon=False,fontsize=8)
    ax=axes[0,1]
    ax.barh(y,report['no_selection_minus_matched_median_pp'],color='#92908c',label='Uncorrected minus matched filter')
    ax.scatter(report['cwb_minus_matched_median_pp'],y,c='#cc6b36',s=30,label='Including cWB minus matched filter',zorder=3)
    ax.axvline(0,c='.4',lw=.8);ax.set(yticks=y,yticklabels=[str(l) for l in [0,4,16]*2],ylabel='Even-shape lambda',xlabel='Median displacement (percentage points)',title='Selection correction and policy changes differ')
    ax.invert_yaxis();ax.legend(frameon=False,fontsize=8);ax.grid(axis='x',alpha=.15)
    ax=axes[1,0];ax.barh(y,sel[:,1]/sel[:,0],color=['#277e9e']*3+['#57916b']*3)
    ax.axvline(1,c='.4',lw=.8,ls='--');ax.set(yticks=y,yticklabels=[str(l) for l in [0,4,16]*2],ylabel='Even-shape lambda',xlabel='Positive / negative component detection efficiency',title='Selection asymmetry within each fixed family')
    ax.invert_yaxis();ax.grid(axis='x',alpha=.15)
    ax=axes[1,1]
    rows=sorted(report['event_concentration'],key=lambda r:min(r['ess_at_medians']))[:15]
    ax.bar(np.arange(len(rows)),[min(r['ess_at_medians']) for r in rows],color='#57916b')
    ax.set(yscale='log',xticks=np.arange(len(rows)),xticklabels=[str(r['index']) for r in rows],xlabel='Source index, ordered by concentration',ylabel='Minimum contribution ESS at the six medians',title='15 most concentrated of 69 full posteriors')
    ax.grid(axis='y',alpha=.15)
    fig.suptitle('Negative effective spin: an explicitly 2PN reference-invariant control',fontsize=15)
    fig.supxlabel('Six restricted sign-shape families; fixed mass/redshift models; uniform g prior. Uncorrected curves are diagnostics.\n2PN conservation and isotropic-measure preservation are assumptions; this is not a formation-channel fraction or full-GR transport.',fontsize=9)
    fig.savefig(BASE/'invariant-control.svg');fig.savefig(BASE/'invariant-control.png',dpi=170);plt.close(fig)
    transport=dict(events=len(sources),rows=sum(r['rows'] for r in sources),
        infinity_nonfinite_rows=sum(r['released_infinity_nonfinite_rows'] for r in sources),
        infinity_sign_disagreement_rows=sum(r['released_infinity_sign_disagreement_rows'] for r in sources),
        maximum_effective_spin_difference=max(r['maximum_released_effective_spin_difference'] for r in sources if r['maximum_released_effective_spin_difference'] is not None),
        maximum_full_ratio_difference=max(x for r in sources for row in r['maximum_full_ratio_transport_difference_at_g01_05_09'] for x in row if x is not None),
        limits='Released-column invariance diagnostic; assumes an invariant isotropic PN measure and does not establish the injection generation reference.')
    (BASE/'reference-invariance-summary.json').write_text(json.dumps(transport,indent=2),encoding='utf-8');print(json.dumps(transport,indent=2))


if __name__=='__main__':main()
