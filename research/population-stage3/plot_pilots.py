"""Plot paired, unchanged event-posterior draws; not population inference."""
from pathlib import Path
import json
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
import numpy as np

BASE = Path(__file__).resolve().parent

def main():
    report = json.loads((BASE / 'pe-pilot-match.json').read_text())
    fig, axes = plt.subplots(2, 2, figsize=(10, 9), sharex=True, sharey=True,
                             layout='constrained')
    with np.load(BASE / 'pilot-tilt-diagnostics.npz', allow_pickle=False) as data:
        for row, event in enumerate(report['events']):
            matched = next(x for x in event['labels'] if x['full_match'])
            for column, component in enumerate([1, 2]):
                ax = axes[row, column]
                k = f'cos_tilt_{component}'
                x = data[f'{event["original_index"]}_{k}']
                y = data[f'{event["original_index"]}_{k}_infinity']
                counts, _, _ = np.histogram2d(x, y, bins=40, range=[[-1,1],[-1,1]])
                shown = ax.imshow(np.ma.masked_equal(counts.T, 0), origin='lower',
                    extent=[-1,1,-1,1], norm=LogNorm(vmin=1,vmax=300), cmap='viridis')
                ax.plot([-1,1],[-1,1],color='#ad5050',lw=1)
                ax.axhline(0,color='#777777',lw=.6)
                ax.axvline(0,color='#777777',lw=.6)
                d = matched['stored_infinity_tilt_diagnostic'][k]
                ax.set_title(f'{event["event"]} | BH {component}', fontsize=11)
                ax.set_xlabel(f'cos(tilt) at stored reference ({matched["f_ref_hz"][0]:g} Hz)')
                ax.set_ylabel('Stored cos(tilt) at infinity\nonly_prec_avg column')
                ax.text(.03,.97,f'Sign crossings: {d["sign_crossing_count"]}/{d["n"]}\n'
                    f'({100*d["sign_crossing_fraction"]:.2f}%)',va='top',transform=ax.transAxes,
                    fontsize=9,bbox=dict(facecolor='white',alpha=.9,edgecolor='none'))
    fig.colorbar(shown,ax=axes.ravel().tolist(),shrink=.75,ticks=[1,10,100,300],
                 label='Paired posterior samples per bin (log scale)')
    fig.suptitle('Two provenance pilots: reference-dependent spin tilts\n'
                 'Event posterior diagnostic; not population fractions or 10-to-20 Hz evolution',fontsize=13)
    fig.savefig(BASE / 'paired-tilt-reference.svg')
    fig.savefig(BASE / 'paired-tilt-reference.png', dpi=145)
    plt.close(fig)

if __name__ == '__main__':
    main()
