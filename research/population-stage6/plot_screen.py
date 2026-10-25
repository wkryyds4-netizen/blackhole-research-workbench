"""Static scientific figure: concentration and two distinct influence diagnostics."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parent


def main():
    data=json.loads((BASE/'catalog-screen.json').read_text())
    byid={e['index']:e for e in data['events']}
    ids=data['rankings']['lowest_median_ess'][:8]
    labels=[byid[i]['event'] for i in ids]
    fig,axes=plt.subplots(2,2,figsize=(12.5,9),layout='constrained')
    fig.suptitle('69-event audit: concentration and influence are different diagnostics',fontsize=15)
    ax=axes[0,0]
    y=np.arange(len(ids))
    for k in range(6):
        values=[byid[i]['settings'][k]['ess_at_median'] for i in ids]
        ax.scatter(values,y+(k-2.5)*.09,s=22,label=f'{"low spin" if k<3 else "baseline spin"}, sigma={[.3,.6,1.2][k%3]}')
    ax.set_yticks(y,labels);ax.set_ylim(10.4,-.5);ax.set_xscale('log')
    ax.set_xlabel('Contribution ESS at each original conditional median')
    ax.set_title('(a) Eight smallest cross-setting median ESS')
    ax.grid(axis='x',alpha=.2);ax.legend(fontsize=7,loc='lower right')
    ax=axes[0,1]
    for e in data['events']:
        color='#d45d32' if e['index'] in data['followup_indices'] else '#48749c'
        ax.scatter(e['min_median_ess'],e['max_abs_leave_event_shift_pp'],c=color,s=24,alpha=.8)
    offsets={19:(5,5),20:(8,12),25:(-5,-22),50:(8,0)}
    for ix in data['followup_indices']:
        e=byid[ix];ax.annotate(e['event'],(e['min_median_ess'],e['max_abs_leave_event_shift_pp']),
                              xytext=offsets[ix],textcoords='offset points',fontsize=7)
    ax.set_xscale('log');ax.set_xlabel('Minimum median ESS across six settings')
    ax.set_ylabel('Largest absolute leave-event-out shift (pp)')
    ax.set_title('(b) Each dot is one stored event');ax.grid(alpha=.2)
    ax=axes[1,0]
    chosen=data['followup_indices'];yy=np.arange(len(chosen))
    ax.barh(yy-.17,[byid[i]['max_abs_drop_row_shift_pp'] for i in chosen],height=.32,label='Remove one influential row',color='#d49b32')
    ax.barh(yy+.17,[byid[i]['max_abs_leave_event_shift_pp'] for i in chosen],height=.32,label='Remove event; N and selection exponent = 68',color='#48749c')
    ax.set_yticks(yy,[byid[i]['event'] for i in chosen]);ax.set_ylim(4.2,-.5)
    ax.set_xlabel('Largest absolute conditional median shift (pp)')
    ax.set_title('(c) Follow-up priority; no exclusions in main analysis')
    ax.legend(fontsize=7,loc='lower right');ax.grid(axis='x',alpha=.2)
    ax=axes[1,1]
    for ix in chosen:
        e=byid[ix]
        values=[100*x['max_row_fraction_at_median'] for x in e['settings']]
        line,=ax.plot(np.arange(3),values[:3],marker='o',label=e['event'])
        ax.plot(np.arange(3,6),values[3:],marker='o',color=line.get_color())
    ax.set_xticks(np.arange(6),['low\n0.3','low\n0.6','low\n1.2','base\n0.3','base\n0.6','base\n1.2'])
    ax.set_xlabel('Spin-magnitude assumption and tilt width sigma')
    ax.set_ylabel('Largest single-row integral contribution (%)')
    ax.set_title('(d) Same 6,913 tutorial rows per event')
    ax.legend(fontsize=7);ax.grid(alpha=.2)
    fig.supxlabel('Fixed six-model family and original tutorial representations; provisional event identities.\nESS ignores chain correlation and missing tails; diagnostic shifts are not uncertainty intervals.',fontsize=9)
    fig.savefig(BASE/'catalog-priority.svg')
    fig.savefig(BASE/'catalog-priority.png',dpi=150)
    plt.close(fig)


if __name__=='__main__':main()
