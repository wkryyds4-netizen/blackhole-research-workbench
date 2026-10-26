"""Fixed54-model extension; no existing scientific cache is rewritten."""
from pathlib import Path
import sys
import numpy as np
from scipy.special import ndtr

ROOT=Path(__file__).resolve().parents[2]
sys.path.insert(0,str(ROOT/'research/population-invariant'))
from run_control import weights as spin_weights
import model

MASS_IDS=['alpha2_peak10','alpha3_peak10','alpha2_nopeak']
KAPPAS=[0,2,4]
MODELS=[dict(index=mi*18+zi*6+si*3+li,mass_id=MASS_IDS[mi],mass_index=mi,
             alpha=alpha,peak_fraction=peak,kappa=kappa,magnitude=[aa,bb],shape=shape,
             key=f'{MASS_IDS[mi]}_kappa{kappa}_beta{aa:g}-{bb:g}_shape{shape}')
        for mi,(alpha,peak) in enumerate([(2,.1),(3,.1),(2,0)])
        for zi,kappa in enumerate(KAPPAS) for si,(aa,bb) in enumerate([(1.5,8.5),(2.,5.)])
        for li,shape in enumerate([0,4,16])]


def mass_pdf(m,alpha,peak):
    normal=(100**(1-alpha)-3**(1-alpha))/(1-alpha)
    gaussian=np.exp(-.5*((m-35)/5)**2)/(5*np.sqrt(2*np.pi)*(ndtr(13)-ndtr(-6.4)))
    return ((1-peak)*m**-alpha/normal+peak*gaussian)*((m>=3)&(m<=100))


def weights(d):
    # Reweight the existing six normalized spin controls by exact target-density
    # ratios in mass and redshift. No PE posterior resampling or truncation.
    base=spin_weights(d)
    reference=mass_pdf(d['mass_1'],2,.1)
    masses=np.column_stack([np.divide(mass_pdf(d['mass_1'],alpha,peak),reference,
                                     out=np.zeros_like(reference),where=reference>0)
                            for alpha,peak in [(2,.1),(3,.1),(2,0)]])
    redshifts=np.column_stack([(1+d['redshift'])**(k-2)*model.ZNORM[2]/model.ZNORM[k] for k in KAPPAS])
    ratios=masses[:,:,None]*redshifts[:,None,:]
    out=(base[:,None,None,:,:]*ratios[:,:,:,None,None]).reshape((len(reference),54,2))
    assert np.isfinite(out).all() and np.all(out>=0)
    return out
