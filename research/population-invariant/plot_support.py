"""Plot declared support sensitivities without implying a measured budget."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.special import expit

BASE=Path(__file__).resolve().parent


def main():
    path=BASE/'support-results.json';r=json.loads(path.read_text())
    v=json.loads((BASE/'support-validation.json').read_text())
    assert v['result_sha256']==hashlib.sha256(path.read_bytes()).hexdigest()
    with np.load(BASE/'support-checkpoints.npz',allow_pickle=False) as cp:
        g=cp['g'];cdf=cp['cdf'][0]
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(1,2,figsize=(12,5.6),layout='constrained')
    colors=['#277e9e','#57916b','#cc6b36']
    B=np.linspace(0,1.5,301)
    for k,(lam,color) in enumerate(zip([0,4,16],colors)):
        low=np.array([np.interp(expit(-b),cdf[k+3],g)-np.interp(expit(b),cdf[k],g) for b in B])*100
        high=np.array([np.interp(expit(b),cdf[k+3],g)-np.interp(expit(-b),cdf[k],g) for b in B])*100
        axes[0].fill_between(B,low,high,color=color,alpha=.09)
        axes[0].plot(B,low,color=color,label=rf'$\lambda={lam}$ lower bound')
        threshold=next(x['critical_log_correction_oscillation'] for x in r['ordering_certificate_thresholds'] if x['policy']==0 and x['lambda_shape']==lam)
        axes[0].scatter([threshold],[0],s=25,color=color)
    axes[0].axhline(0,color='.4',lw=.8)
    axes[0].set(xlabel='Assumed log-correction oscillation B',ylabel='Amplitude-model median difference in g (pp)',
        title='Outer envelopes under a declared error budget\nB has not been measured')
    axes[0].legend(fontsize=8,frameon=False)
    finite=r['event_completions'][:-1];etas=np.array([x['eta'] for x in finite])
    shifts=np.array([x['median_shift_pp'][0] for x in finite])
    for k,(lam,color) in enumerate(zip([0,4,16],colors)):
        axes[1].plot(etas,shifts[:,k],color=color,marker='o',ms=3,ls='--',label=rf'Low spin, $\lambda={lam}$')
        axes[1].plot(etas,shifts[:,k+3],color=color,marker='o',ms=3,label=rf'Baseline spin, $\lambda={lam}$')
    axes[1].set_xscale('symlog',linthresh=1e-4,linscale=.35)
    axes[1].set(xlabel=r'Assumed missing integral ratio $\eta=\mu_-/\mu_+$',ylabel='Shift of the g posterior median (pp)',
        title='Complete one event at integral level\nGW191204_171526; no physical likelihood is inferred')
    axes[1].legend(fontsize=8,frameon=False,ncol=2)
    for ax in axes:ax.grid(alpha=.15)
    fig.suptitle('Missing support: conditional bounds and a declared stress family',fontsize=14)
    fig.supxlabel('Matched-filter policy; all 69 events retained. Shading is a sensitivity envelope, not a credible interval.\nEach stress value is an assumption. Neither B nor the event completion is estimated by finite-bank resampling.',fontsize=9)
    fig.savefig(BASE/'support-sensitivity.svg');fig.savefig(BASE/'support-sensitivity.png',dpi=170);plt.close(fig)


if __name__=='__main__':main()
