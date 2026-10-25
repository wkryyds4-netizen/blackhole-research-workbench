"""Full-catalog candidate metadata overview, with verification scope visible."""
from pathlib import Path
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap

BASE=Path(__file__).resolve().parent


def main():
    result=json.loads((BASE/'catalog-summary.json').read_text());events=result['events'];n=len(events)
    frequencies=np.array([e['native_reference_hz'] for e in events])
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(11.8,6.9),gridspec_kw={'height_ratios':[1.7,1]},layout='constrained')
    ax=axes[0,0];values=sorted(set(frequencies));colors=['#277e9e','#cc6b36','#57916b','#966b9d']
    for f,c in zip(values,colors):
        mask=frequencies==f;ax.scatter(np.flatnonzero(mask)+1,frequencies[mask],c=c,s=35,label=f'{f:g} Hz ({int(mask.sum())})')
    ax.set(xlabel='Candidate event, chronological order',ylabel='Native spin reference frequency (Hz)',title='References must be read event by event',yticks=values,xlim=(0,n+1))
    ax.legend(frameon=False,fontsize=9);ax.grid(axis='y',alpha=.15)
    ax=axes[0,1]
    ax.bar([str(v).rstrip('0').rstrip('.') for v in values],[(frequencies==v).sum() for v in values],color=colors[:len(values)])
    ax.set(xlabel='Native reference frequency (Hz)',ylabel='Candidate released analyses',title='One common default would overwrite real metadata')
    for j,v in enumerate(values):ax.text(j,int((frequencies==v).sum())+.5,str(int((frequencies==v).sum())),ha='center')
    ax=axes[1,0];checks=['native_spin_fields_complete','infinity_tilts_present','full_core_prior_strings_present']
    matrix=np.array([[e[k] for e in events] for k in checks])
    ax.imshow(matrix,aspect='auto',interpolation='nearest',cmap=ListedColormap(['#ddd8cf','#277e9e']),vmin=0,vmax=1)
    ax.set(yticks=range(3),yticklabels=['Native spins + relative azimuth','Released infinity tilts','Core analytic prior strings'],xlabel='Candidate event, chronological order',xticks=[0,9,29,49,68],xticklabels=['1','10','30','50','69'],title='Blue: metadata present; gray: declaration absent')
    ax=axes[1,1];ax.axis('off')
    verified=result['counts']['native_posterior_source_identity_verified']
    text=('69 fixed official candidate files\n'
          'C01:IMRPhenomXPHM inventories\n\n'
          'Metadata availability is now testable.\n'
          f'{verified} native sources verified; {n-verified} pending.\n'
          'Injection spin references remain unresolved.')
    ax.text(.02,.95,text,va='top',fontsize=11,linespacing=1.5)
    fig.suptitle('Full-catalog spin-reference and prior metadata inventory',fontsize=15)
    fig.supxlabel('Metadata and analysis-lineage inventory; common-reference tilt inference is not claimed.',fontsize=9)
    fig.savefig(BASE/'catalog-metadata.svg');fig.savefig(BASE/'catalog-metadata.png',dpi=170);plt.close(fig)


if __name__=='__main__':main()
