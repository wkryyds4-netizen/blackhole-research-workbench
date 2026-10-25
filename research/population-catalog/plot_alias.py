"""Visualize label-derived differences and retain the failed numerical gate."""
from pathlib import Path
import argparse
import json
import sys
import numpy as np

BASE=Path(__file__).resolve().parent
ROOT=BASE.parents[1]


def prepare():
    import h5py
    sys.path.insert(0,str(ROOT/'research/population-stage7'))
    from extract_events import join,SOURCE
    rec=next(r for r in json.loads((BASE/'download-results.json').read_text()) if r['index']==14)
    with h5py.File(BASE/'source-data'/rec['source']['file'],'r') as h:
        a=h['C01:IMRPhenomXPHM/posterior_samples'][:];b=h['C01:Mixed/posterior_samples'][:]
        rows,info=join([a[k] for k in SOURCE],[b[k] for k in SOURCE]);assert info['matches']==len(a) and info['ambiguous']==0
    np.savez_compressed(BASE/'alias-plot-cache.npz',canonical=np.array([a[f'cos_tilt_{j}_infinity_only_prec_avg'] for j in [1,2]]),alias=np.array([b[f'cos_tilt_{j}_infinity_only_prec_avg'][rows] for j in [1,2]]))


def plot():
    import matplotlib
    matplotlib.use('Agg')
    import matplotlib.pyplot as plt
    with np.load(BASE/'alias-plot-cache.npz',allow_pickle=False) as cp:a=cp['canonical'];b=cp['alias']
    diagnosis=json.loads((BASE/'alias-evolution-diagnosis.json').read_text());validation=json.loads((BASE/'alias-heldout-validation.json').read_text())
    historical=json.loads((BASE/'historical-defaults-validation.json').read_text())
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(11.8,8),layout='constrained')
    delta=b-a
    ax=axes[0,0]
    for j,c in enumerate(['#277e9e','#cc6b36']):ax.hist(delta[j],bins=np.linspace(-.14,.14,91),histtype='step',lw=1.6,color=c,label=f'Component {j+1}')
    ax.axvline(0,color='0.5',lw=.8);ax.set(yscale='log',xlabel='Mixed minus IMRPhenomXPHM infinity cosine',ylabel='Posterior rows per bin',title='147,856 matching native rows; different derived tilts')
    ax.legend(frameon=False);ax.grid(axis='y',alpha=.15)
    ax=axes[0,1];hb=ax.hexbin(a[0],delta[0],gridsize=60,mincnt=1,bins='log',cmap='Blues')
    ax.axhline(0,color='0.4',lw=.8);ax.set(xlabel='IMRPhenomXPHM infinity cosine, component 1',ylabel='Mixed minus IMRPhenomXPHM',title='Largest component-1 difference = 0.12980')
    fig.colorbar(hb,ax=ax,label='Matching rows per hexagon')
    ax=axes[1,0];ds=diagnosis['summaries'];freq=[v['frequency'] for v in ds]
    ax.plot(freq,[v['max_difference_from_canonical'] for v in ds],marker='o',c='#277e9e',label='Compared with IMRPhenomXPHM')
    ax.plot(freq,[v['max_difference_from_alias'] for v in ds],marker='s',c='#cc6b36',label='Compared with Mixed')
    ax.set(yscale='log',xticks=freq,xlabel='Reference frequency used in our transport (Hz)',ylabel='Maximum infinity-cosine residual',title='20 diagnostic rows: the two conventions separate')
    ax.legend(frameon=False,fontsize=8);ax.grid(alpha=.15)
    ax=axes[1,1];v=validation['summaries'];x=np.arange(2);width=.24
    ax.bar(x-width,[r['max_target_cosine_difference'] for r in v],width,color='#cc6b36',label='Independent solver vs release')
    ax.bar(x,[r['max_final_tolerance_difference'] for r in v],width,color='#277e9e',label='Final tolerance-pair change')
    ax.bar(x+width,[r['maximum_released_cosine_residual'] for r in historical['summaries']],width,color='#57916b',label='Historical defaults vs release')
    ax.axhline(1e-6,c='0.3',ls='--',label='Original 1e-6 gate')
    ax.set(xticks=x,xticklabels=['16 Hz vs Phenom','20 Hz vs Mixed'],ylabel='Maximum absolute cosine difference',ylim=(0,1.85e-6),title='70 fixed rows: historical defaults explain the residual')
    ax.ticklabel_format(axis='y',style='sci',scilimits=(0,0));ax.legend(frameon=False,fontsize=8,loc='upper left')
    fig.suptitle('The same native posterior can carry different reference-dependent derived columns',fontsize=14)
    fig.supxlabel('GW190421_213856. Original 1e-6 gate remains failed. Historical-default replay is a follow-up on the same rows, not a full historical runtime.',fontsize=9)
    fig.savefig(BASE/'alias-reference.svg');fig.savefig(BASE/'alias-reference.png',dpi=170);plt.close(fig)


if __name__=='__main__':
    parser=argparse.ArgumentParser();parser.add_argument('--prepare',action='store_true');args=parser.parse_args()
    if args.prepare:prepare()
    else:plot()
