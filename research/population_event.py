"""Weighted finite posterior representations; no resampling or support inference.

For posterior weights s and its matching PE density pi, event integrals are
sum(s * population_density / pi) / sum(s). Omitted PE normalization constants
cancel only in conditional hyperposterior shapes, not absolute evidence/rates.
"""
from dataclasses import dataclass, replace
from pathlib import Path
import copy
import json
import numpy as np
from scipy.special import logsumexp


def _log_vector(value, n, name, finite=False):
    a = np.array(value, dtype=float, copy=True)
    if a.shape != (n,) or np.any(np.isnan(a)) or np.any(np.isposinf(a)):
        raise ValueError(f'{name}: expected {n} finite or negative-infinite log values')
    if finite and not np.all(np.isfinite(a)):
        raise ValueError(f'{name}: zero/unknown PE density on retained rows is unsupported')
    if not np.any(np.isfinite(a)):
        raise ValueError(f'{name}: no positive mass')
    a.flags.writeable = False
    return a


@dataclass(frozen=True)
class WeightedEvent:
    samples: dict
    row_ids: np.ndarray
    log_pe_prior: np.ndarray
    log_sample_weight: np.ndarray
    metadata: dict

    def __post_init__(self):
        ids = np.array(self.row_ids, copy=True)
        if ids.ndim != 1 or not len(ids) or ids.dtype.kind not in 'iu' or np.any(ids < 0):
            raise ValueError('row_ids must be a nonempty vector of nonnegative source row indices')
        if len(np.unique(ids)) != len(ids):
            raise ValueError('repeated source rows must be aggregated as weights, not independent rows')
        ids.flags.writeable = False
        n = len(ids)
        arrays = {}
        for key, value in self.samples.items():
            a = np.array(value, dtype=float, copy=True)
            if a.shape != (n,) or not np.all(np.isfinite(a)):
                raise ValueError(f'invalid coordinate array: {key}')
            a.flags.writeable = False
            arrays[key] = a
        if not arrays:
            raise ValueError('coordinates required')
        meta = copy.deepcopy(self.metadata)
        for key in ('event', 'analysis_label', 'source', 'row_basis', 'measure',
                    'reference', 'cosmology', 'prior', 'support_status', 'identity_status'):
            if key not in meta:
                raise ValueError(f'missing metadata: {key}')
        if set(meta['measure']['coordinates']) != set(arrays):
            raise ValueError('measure coordinates must match supplied arrays')
        source = meta['source']
        if not source.get('file') or not source.get('url') or len(source.get('sha256', '')) != 64:
            raise ValueError('source file, URL and SHA256 required')
        try:
            int(source['sha256'], 16)
        except ValueError as error:
            raise ValueError('invalid SHA256') from error
        json.dumps(meta, allow_nan=False)
        object.__setattr__(self, 'row_ids', ids)
        object.__setattr__(self, 'samples', arrays)
        object.__setattr__(self, 'metadata', meta)
        object.__setattr__(self, 'log_pe_prior', _log_vector(self.log_pe_prior, n, 'PE prior', True))
        object.__setattr__(self, 'log_sample_weight', _log_vector(self.log_sample_weight, n, 'sample weight'))

    def log_integral(self, log_population_density):
        """First axis is source row; trailing axes may be population settings."""
        p = np.asarray(log_population_density, dtype=float)
        if p.ndim < 1 or p.shape[0] != len(self.row_ids) or np.any(np.isnan(p)) or np.any(np.isposinf(p)):
            raise ValueError('invalid population log density')
        shape = (len(self.row_ids),) + (1,) * (p.ndim - 1)
        terms = p + (self.log_sample_weight - self.log_pe_prior).reshape(shape)
        return logsumexp(terms, axis=0) - logsumexp(self.log_sample_weight)

    def diagnostics(self, log_population_density):
        p = np.asarray(log_population_density, dtype=float)
        if p.shape != self.row_ids.shape or np.any(np.isnan(p)) or np.any(np.isposinf(p)):
            raise ValueError('diagnostics require one population density per row')
        terms = p + self.log_sample_weight - self.log_pe_prior
        total = logsumexp(terms)
        if not np.isfinite(total):
            raise ValueError('zero event integral has no normalized contributions')
        fraction = np.exp(terms - total)
        sample = np.exp(self.log_sample_weight - logsumexp(self.log_sample_weight))
        return dict(n_rows=len(fraction), contribution_ess=float(1 / np.sum(fraction**2)),
                    representation_weight_ess=float(1 / np.sum(sample**2)),
                    max_contribution=float(fraction.max()),
                    largest_source_row=int(self.row_ids[fraction.argmax()]),
                    n_for_half_integral=int(np.searchsorted(np.cumsum(np.sort(fraction)[::-1]), .5) + 1))

    def reexpress_prior(self, log_new_prior, prior_description, support_status):
        """Same coordinates/rows/likelihood, with continuous posterior reweighting.

        Cannot recover missing support, change spin reference, or compute the
        normalization of a physical PE prior. Positive density is required on
        every retained row; the caller must separately assess target support.
        """
        new = _log_vector(log_new_prior, len(self.row_ids), 'new PE prior', True)
        meta = copy.deepcopy(self.metadata)
        meta.setdefault('prior_history', []).append(meta['prior'])
        meta['prior'] = copy.deepcopy(prior_description)
        meta['support_status'] = support_status
        return replace(self, log_pe_prior=new,
                       log_sample_weight=self.log_sample_weight + new - self.log_pe_prior,
                       metadata=meta)

    def save(self, path):
        """Local cache: JSON metadata and numerical arrays, never Python pickle."""
        np.savez_compressed(path, row_ids=self.row_ids, log_pe_prior=self.log_pe_prior,
                            log_sample_weight=self.log_sample_weight,
                            metadata=np.array(json.dumps(self.metadata, allow_nan=False)),
                            **{'sample_' + k: v for k, v in self.samples.items()})

    @classmethod
    def load(cls, path):
        with np.load(Path(path), allow_pickle=False) as a:
            return cls({k[7:]: a[k] for k in a.files if k.startswith('sample_')},
                       a['row_ids'], a['log_pe_prior'], a['log_sample_weight'],
                       json.loads(str(a['metadata'])))
