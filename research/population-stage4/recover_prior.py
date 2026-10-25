"""Recover the tutorial density formula; never execute fetched source code."""
from pathlib import Path
import hashlib
import json
import numpy as np
from scipy.interpolate import interp1d
from astropy.cosmology import Planck15, FlatLambdaCDM

BASE = Path(__file__).resolve().parent
ROOT = BASE.parents[1]
LAL = FlatLambdaCDM(H0=67.9, Om0=.3065, name='Planck15_lal')

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def redshift_grid(archive):
    maximum = max(float(archive[f'{i}_redshift'].max()) for i in range(73))
    return np.linspace(0, maximum * 1.01, 1000)

def density_function(grid, cosmology=Planck15, euclidean=False):
    if euclidean:
        dl = cosmology.luminosity_distance(grid).value
        unnorm = dl**2 * (dl/(1+grid) + (1+grid)*cosmology.hubble_distance.value/cosmology.efunc(grid))
    else:
        unnorm = cosmology.differential_comoving_volume(grid).value * 4*np.pi/(1+grid)
    norm = np.trapezoid(unnorm, grid)
    return interp1d(grid, unnorm/norm, bounds_error=True), norm

def prior(draws, interp):
    z = draws['redshift']
    return interp(z) * (1+z)**2 * draws['mass_1']/4

def shape_residual(stored, proposed):
    ratio = stored/proposed
    scale = float(np.median(ratio))
    return dict(median_scale=scale, maximum_relative_shape_residual=float(np.max(abs(ratio/scale-1))))

def main():
    source = BASE/'sources/pipe-april2024.py.txt'
    records = []
    with np.load(ROOT/'data/posteriors.npz',allow_pickle=False) as a:
        grid = redshift_grid(a)
        fun, norm = density_function(grid)
        euc, _ = density_function(grid, euclidean=True)
        for i in range(73):
            d = {k:a[f'{i}_{k}'] for k in ['mass_1','redshift','prior']}
            p = prior(d,fun)
            err = float(np.max(abs(p/d['prior']-1)))
            assert err < 1e-12
            records.append(dict(index=i,n=len(p),maximum_relative_error=err,
                negative_controls={
                    'omit_mass_jacobian':shape_residual(d['prior'],p/d['mass_1']),
                    'omit_redshift_mass_jacobian':shape_residual(d['prior'],p/(1+d['redshift'])**2),
                    'euclidean_distance':shape_residual(d['prior'],prior(d,euc))}))
    result = dict(schema='blackhole-population-prior-recovery/1',
        source_commit='5dbcb49d602393f23d91e4d2b9a4cbae52d843a8', source_sha256=sha(source),
        posterior_sha256=sha(ROOT/'data/posteriors.npz'), grid=dict(n=1000,min=0,max=float(grid[-1]),
        rule='maximum stored redshift over all 73 arrays, multiplied by 1.01'),
        redshift_normalization_mpc3=norm,
        formula='interp[pz_Planck15](z) * (1+z)^2 * mass_1 / 4; pz proportional to 4pi*dVc/dz/(1+z)',
        measure='d mass_1_source d mass_ratio d redshift d a_1 d a_2 d cos_tilt_1 d cos_tilt_2',
        total_samples=sum(x['n'] for x in records), maximum_relative_error=max(x['maximum_relative_error'] for x in records),
        events=records,
        scope='Reproduces stored values on observed rows; does not prove historical command/version, all event priors, or unobserved support.',
        cosmology_caution='Recovered density uses Astropy Planck15; pilot PE files label Planck15_lal. These are not identical cosmologies.')
    (BASE/'prior-recovery.json').write_text(json.dumps(result,indent=2),encoding='utf-8')
    print(json.dumps({k:v for k,v in result.items() if k!='events'},indent=2))

if __name__=='__main__':
    main()
