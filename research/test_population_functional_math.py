"""Independent mathematical checks; no survey/formation-channel validity claim.

Run with the scientific environment:
  simulator/.adapters-venv/Scripts/python.exe research/test_population_functional_math.py

Numerical quadrature deliberately constructs the density independently of model.py.
"""
import json
import math
from pathlib import Path
import unittest

import numpy as np
from scipy.integrate import quad
from scipy.special import ndtr

ROOT = Path(__file__).resolve().parents[1]


def component_cdf(x, sigma):
    return (ndtr((x - 1) / sigma) - ndtr(-2 / sigma)) / (0.5 - ndtr(-2 / sigma))


def affine(sigma, quantity):
    if quantity == "negative":
        return 0.5, component_cdf(0, sigma) - 0.5
    u = (1 - math.cos(math.pi / 18)) / 2
    return u, 1 - component_cdf(math.cos(math.pi / 18), sigma) - u


def physical_probability(f, sigma, quantity):
    b, s = affine(sigma, quantity)
    return b + s * f


class FunctionalMathTests(unittest.TestCase):
    def test_pdf_integrals_and_joint_semantics(self):
        for sigma in [0.3, 0.6, 1.2]:
            z = quad(lambda c: math.exp(-0.5 * ((c - 1) / sigma) ** 2), -1, 1,
                     epsabs=1e-13)[0]
            for f in [0, 0.17, 0.51, 1]:
                pdf = lambda c: f * math.exp(-0.5 * ((c - 1) / sigma) ** 2) / z + (1-f)/2
                self.assertAlmostEqual(quad(pdf, -1, 1)[0], 1, places=12)
                self.assertAlmostEqual(quad(pdf, -1, 0)[0], physical_probability(f, sigma, "negative"), places=12)
                self.assertAlmostEqual(quad(pdf, math.cos(math.pi/18), 1)[0], physical_probability(f, sigma, "aligned10"), places=12)
                b = component_cdf(0, sigma)
                joint_both = f*b*b + (1-f)/4
                joint_either = f*(2*b-b*b) + (1-f)*3/4
                self.assertAlmostEqual(joint_either + joint_both, 2*physical_probability(f,sigma,"negative"), places=12)

    def test_kl_and_standardized_change_invariant(self):
        f = (np.arange(10001) + 0.5) / 10001
        prior = 6*f*(1-f)
        prior /= prior.mean()
        posterior = prior * np.exp(-0.5*((f-0.21)/0.12)**2)
        posterior /= posterior.mean()
        df = 1/len(f)
        kl_f = np.sum(posterior*np.log(posterior/prior))*df
        mean = lambda x,p,dx: np.sum(x*p)*dx
        var = lambda x,p,dx: np.sum((x-mean(x,p,dx))**2*p)*dx
        effect_f = abs(mean(f,posterior,df)-mean(f,prior,df))/math.sqrt(var(f,prior,df))
        shrink_f = math.sqrt(var(f,posterior,df)/var(f,prior,df))
        for sigma in [.3,.6,1.2]:
            for quantity in ["negative","aligned10"]:
                b,s = affine(sigma,quantity)
                x = b+s*f
                dx = abs(s)*df
                pg, qg = posterior/abs(s), prior/abs(s)
                kl_g = np.sum(pg*np.log(pg/qg))*dx
                self.assertAlmostEqual(kl_g, kl_f, places=12)
                self.assertAlmostEqual(abs(mean(x,pg,dx)-mean(x,qg,dx))/math.sqrt(var(x,qg,dx)),effect_f,places=11)
                self.assertAlmostEqual(math.sqrt(var(x,pg,dx)/var(x,qg,dx)),shrink_f,places=11)

    def test_common_prior_constant_likelihood(self):
        for quantity in ["negative", "aligned10"]:
            supports = [sorted([physical_probability(0,s,quantity),physical_probability(1,s,quantity)]) for s in [.3,.6,1.2]]
            lo, hi = max(x[0] for x in supports), min(x[1] for x in supports)
            self.assertGreater(hi,lo)
            x = lo+(np.arange(2001)+.5)/2001*(hi-lo)
            for sigma in [.3,.6,1.2]:
                b,s = affine(sigma,quantity)
                f = (x-b)/s
                self.assertTrue(np.all((f>=0)&(f<=1)))
                p_f = np.full_like(f,abs(s)/(hi-lo))
                p_x = p_f/abs(s)
                self.assertTrue(np.allclose(p_x,1/(hi-lo),rtol=1e-14))

    def test_existing_scenario_likelihood_transform(self):
        source = ROOT/'results/model.json'
        if not source.exists():
            self.skipTest('No existing scenario moment artifact')
        data = json.loads(source.read_text(encoding='utf-8'))
        f = np.array(data['f_grid'])
        for key, record in data['scenarios'].items():
            m=np.asarray(record['event_moments']); a,b=record['selection'][:2]
            ll=np.log(f[:,None]*m[:,0]+(1-f[:,None])*m[:,1]).sum(axis=1)-len(m)*np.log(f*a+(1-f)*b)
            p=np.exp(ll-np.max(ll)); p/=p.sum()
            prior=np.full(len(f),1/len(f))
            kl=np.sum(p[p>0]*np.log(p[p>0]/prior[p>0]))
            sigma=[.3,.6,1.2][int(key.split('-')[0])]
            for quantity in ['negative','aligned10']:
                intercept,slope=affine(sigma,quantity)
                x=intercept+slope*f
                dx=abs(slope)/len(f)
                pd=p/dx; qd=prior/dx
                kl_transformed=np.sum(pd[p>0]*np.log(pd[p>0]/qd[p>0]))*dx
                self.assertAlmostEqual(kl_transformed,kl,places=11)
                self.assertGreater(np.ptp(x),0)


if __name__ == '__main__':
    unittest.main(verbosity=2)
