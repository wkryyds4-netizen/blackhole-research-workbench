"""Deterministic shifts, conditional thinning intervals and missing high weights."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BASE=Path(__file__).resolve().parent

def main():
    j=json.loads((BASE/'sensitivity-results.json').read_text())
    s=j['scenarios'];p=np.array(s['phenom_both']['median_shift_from_original_pp'])
    mix=np.array(s['mixed_both']['median_shift_from_original_pp'])-p
    ranges=np.array(j['finite_release_thinning']['median_minus_full_phenom_pp_quantiles'])
    fig,(ax,tail)=plt.subplots(1,2,figsize=(12,5.7),gridspec_kw={'width_ratios':[1.65,1]},layout='constrained')
    y=np.arange(6)
    ax.hlines(y,ranges[:,0],ranges[:,2],color='#bbbfc6',lw=7,label='5-95% finite-release thinning range')
    ax.scatter(-p,y-.10,color='#b96b16',marker='o',s=43,label='Original tutorial subsets',zorder=3)
    ax.scatter(mix,y+.10,color='#24669b',marker='D',s=40,label='Released Mixed, assumed common prior',zorder=3)
    ax.axvline(0,color='#555',lw=.9)
    ax.set_yticks(y,[f'Beta(1.5,8.5), sigma={v}' for v in [.3,.6,1.2]]+[f'Beta(2,5), sigma={v}' for v in [.3,.6,1.2]])
    ax.invert_yaxis();ax.set_ylim(6.7,-.6)
    ax.set_xlabel('Change in r median relative to full IMRPhenomXPHM\n(percentage points; two events replaced)')
    ax.set_title('Conditional population-output sensitivity',fontsize=12)
    ax.legend(loc='lower right',fontsize=8,frameon=False)
    t=j['weight_tail_audit'];top=t['top_five']
    tail.barh(np.arange(5),[100*x['fraction_of_integral'] for x in top],
              color=['#24669b' if x['included_in_original_subset'] else '#b96b16' for x in top])
    tail.set_yticks(np.arange(5),[f'Row {x["source_row"]}' for x in top]);tail.invert_yaxis()
    tail.set_xlabel('Contribution to this event integral (%)')
    tail.set_title(f'GW190805: five largest weights\nLow-spin sigma=0.3, r={t["r"]:.4f}',fontsize=11)
    for i,x in enumerate(top):
        tail.text(100*x['fraction_of_integral']+.35,i,'kept' if x['included_in_original_subset'] else 'absent',va='center',fontsize=9)
    tail.set_xlim(0,25)
    tail.text(.02,.025,'Absent = missing from original 6,913 draws.\nLargest draw is present in full 11,786 release.',
              transform=tail.transAxes,va='bottom',fontsize=8)
    fig.suptitle('Prior recovery and two-event stress test\nThinning ranges are not astrophysical credible intervals or waveform error bars',fontsize=12)
    fig.savefig(BASE/'prior-and-pe-sensitivity.svg')
    fig.savefig(BASE/'prior-and-pe-sensitivity.png',dpi=150)
    plt.close(fig)

if __name__=='__main__':main()
