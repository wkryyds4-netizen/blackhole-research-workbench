"""Visualize the current validated native-source replacement snapshot."""
from pathlib import Path
import hashlib
import json
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

BASE=Path(__file__).resolve().parent


def main():
    path=BASE/'native-catalog-results.json';result=json.loads(path.read_text())
    validation=json.loads((BASE/'native-catalog-validation.json').read_text())
    assert validation['result_sha256']==hashlib.sha256(path.read_bytes()).hexdigest()
    integrals=json.loads((BASE/'native-event-integrals.json').read_text())
    rows=[r for r in integrals['events'] if r['index'] in result['newly_recomputed_source_indices']]
    labels=[f'{group}, sigma={sigma}' for group in ['Beta(1.5,8.5)','Beta(2,5)'] for sigma in [.3,.6,1.2]]
    old=np.asarray(result['baseline_quantiles'])*100;new=np.asarray(result['quantiles'])*100
    n=result['completed_full_source_events'];total=result['fixed_total_events']
    plt.rcParams.update({'font.size':10,'svg.fonttype':'none','axes.spines.top':False,'axes.spines.right':False})
    fig,axes=plt.subplots(2,2,figsize=(12.5,8.4),layout='constrained')
    y=np.arange(6)
    ax=axes[0,0]
    for q,dy,c,label in [(old,-.13,'#92908c','Previous 5-source baseline'),(new,.13,'#277e9e',f'{n} full sources, {total-n} tutorial arrays')]:
        ax.errorbar(q[:,1],y+dy,xerr=np.array([q[:,1]-q[:,0],q[:,2]-q[:,1]]),fmt='o',markersize=4,c=c,capsize=2,label=label)
    ax.set(yticks=y,yticklabels=labels,xlabel='r = P(cos tilt 1 < 0), percent',title='Conditional medians and 90% credible intervals')
    ax.invert_yaxis();ax.legend(frameon=False,fontsize=8);ax.grid(axis='x',alpha=.15)
    ax=axes[0,1]
    ax.barh(y,result['median_shift_pp'],color=['#277e9e']*3+['#cc6b36']*3)
    ax.axvline(0,c='.4',lw=.7);ax.set(yticks=y,yticklabels=[str(s) for s in [.3,.6,1.2]*2],xlabel='Median shift from 5-source baseline (pp)',ylabel='Tilt width sigma',title='Small displacements, shown on their own scale')
    ax.invert_yaxis();ax.grid(axis='x',alpha=.15)
    ax=axes[1,0];ordered=sorted(rows,key=lambda r:r['min_component_ess'])[:20]
    ax.bar(np.arange(len(ordered)),[r['min_component_ess'] for r in ordered],color='#57916b')
    ax.set(yscale='log',xticks=np.arange(len(ordered)),xticklabels=[str(r['index']) for r in ordered],xlabel='Restored source index, ordered by concentration',ylabel='Minimum component contribution ESS',title=f'{len(ordered)} most concentrated of {len(rows)} recomputed sources\nMinimum over all 12 model components')
    ax.grid(axis='y',alpha=.15)
    ax=axes[1,1];ax.axis('off')
    gaps=result['width_gaps_pp'];errors=[r['mean_relative_error'] for r in validation['checks']]
    note=(f'Original-prior full sources: {n}/{total}\n'
          f'Recomputed rows: {result["source_rows"]:,}\n'
          'Previously restored sources retained in baseline.\n\n'
          f'Median width gaps: {gaps[0]:.5f} / {gaps[1]:.5f} pp\n'
          f'Independent event-mean discrepancy: {max(errors):.2e}\n'
          f'Independent quantile discrepancy: {validation["full_product_quantile_error"]:.2e}\n\n'
          'ESS describes contribution concentration;\n'
          'it excludes PE-chain autocorrelation.\n'
          'No source row or event has been discarded.')
    ax.text(.02,.96,note,va='top',linespacing=1.5)
    fig.suptitle(f'Original-prior restoration: {n} full sources in a fixed {total}-event likelihood',fontsize=15)
    fig.supxlabel('Native spin references; fixed matched-filter selection and six conditional models. No common-reference or formation-channel inference.',fontsize=9)
    fig.savefig(BASE/'native-catalog-comparison.svg');fig.savefig(BASE/'native-catalog-comparison.png',dpi=170);plt.close(fig)


if __name__=='__main__':main()
