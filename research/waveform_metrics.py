"""Finite-record waveform comparisons; no parameter estimation or detector network.

See waveform-metrics-method.md for conventions and interpretation.
"""
from __future__ import annotations

import numpy as np
from scipy.fft import next_fast_len
from scipy.signal.windows import tukey
from scipy.optimize import minimize_scalar

MTSUN_S = 4.9254909476412675e-6
PSD_SOURCE = "https://lscsoft.docs.ligo.org/lalsuite/lal/_l_a_l_adv_l_i_g_o_psd_8c_source.html"


def aligo_historical_psd(f_hz):
    """Historical LALAdvLIGOPsd analytic shape times 1e-49 /Hz.

    For controlled sensitivity experiments only, not measured detector noise.
    Our deliberately limited policy domain is 10--2048 Hz.
    """
    f = np.asarray(f_hz, dtype=float)
    if np.any(~np.isfinite(f)) or np.any((f < 10) | (f > 2048)):
        raise ValueError("Historical PSD is restricted here to 10 <= f <= 2048 Hz")
    x = f / 215.0
    return 1e-49 * (x**(-4.14) - 5/x**2 + 111*(1-x**2+x**4/2)/(1+x**2/2))


def common_grid(t_a, a, t_b, b, dt, *, interval=None):
    """Linear interpolation on overlap only; units inherited from supplied t.

    No extrapolation, frame rotation, peak alignment or anti-alias filtering.
    The caller must establish a common frame/time convention and convergence
    with sampling before interpreting the returned waveforms.
    """
    arrays = []
    for t, y in [(t_a, a), (t_b, b)]:
        t, y = np.asarray(t, float), np.asarray(y)
        if t.ndim != 1 or y.shape != t.shape or len(t) < 4:
            raise ValueError("Each waveform needs matching one-dimensional time/sample arrays")
        if not np.all(np.isfinite(t)) or not np.all(np.isfinite(y)) or np.any(np.diff(t) <= 0):
            raise ValueError("Finite samples and strictly increasing times required")
        arrays.append((t, y))
    if not np.isfinite(dt) or dt <= 0:
        raise ValueError("dt must be positive and finite")
    lo, hi = max(arrays[0][0][0], arrays[1][0][0]), min(arrays[0][0][-1], arrays[1][0][-1])
    if interval is not None:
        if len(interval) != 2 or not np.all(np.isfinite(interval)) or interval[0] >= interval[1]:
            raise ValueError("Invalid interval")
        if interval[0] < lo or interval[1] > hi:
            raise ValueError("Requested interval exceeds common support")
        lo, hi = interval
    if hi <= lo or (hi-lo)/dt < 3:
        raise ValueError("Insufficient common support")
    grid = lo + np.arange(int(np.floor((hi-lo)/dt))+1)*dt
    return grid, *(np.interp(grid, t, y) for t, y in arrays)


def compare_waveforms(a, b, dt_s, *, fmin_hz, fmax_hz,
                      signal_kind="complex_mode", psd="white",
                      taper_fraction=0.1, max_time_shift_s=0.05,
                      oversample=4, padding_factor=4, refine=True):
    """Normalized time/global-phase optimized overlap of same-grid records.

    For complex_mode: weighted two-sided complex correlation using PSD(|f|).
    For real_strain: positive-frequency correlation and Fourier-quadrature phase.
    Returned time_shift_s and phase_shift_rad apply to b(t-time_shift_s).
    The phase multiplies the complex mode by exp(i phase_shift_rad); for real
    strain it rotates the positive-frequency quadrature, NOT orbital phase.
    """
    a, b = np.asarray(a), np.asarray(b)
    if a.ndim != 1 or a.shape != b.shape or len(a) < 8:
        raise ValueError("Waveforms must be same-length one-dimensional arrays, at least 8 samples")
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)):
        raise ValueError("Waveforms contain nonfinite values")
    if not np.isfinite(dt_s) or dt_s <= 0:
        raise ValueError("dt_s must be positive and finite")
    nyquist = .5/dt_s
    if not (np.isfinite(fmin_hz) and np.isfinite(fmax_hz) and 0 < fmin_hz < fmax_hz < nyquist):
        raise ValueError("Require 0 < fmin_hz < fmax_hz < Nyquist; DC/Nyquist bins excluded")
    if not np.isfinite(taper_fraction) or not 0 <= taper_fraction <= 1:
        raise ValueError("taper_fraction must be in [0,1]")
    if not np.isfinite(max_time_shift_s) or not 0 <= max_time_shift_s < (len(a)-1)*dt_s/2:
        raise ValueError("Time search must be nonnegative and less than half finite record duration")
    if int(oversample) != oversample or not 1 <= oversample <= 32:
        raise ValueError("oversample must be an integer in [1,32]")
    if not np.isfinite(padding_factor) or not 2 <= padding_factor <= 32:
        raise ValueError("padding_factor must be in [2,32]")
    if signal_kind not in ["complex_mode", "real_strain"]:
        raise ValueError("Unknown signal_kind")
    if signal_kind == "real_strain" and (np.iscomplexobj(a) or np.iscomplexobj(b)):
        raise ValueError("real_strain requires real inputs; project explicitly first")
    n = len(a)
    nfft = next_fast_len(int(np.ceil(padding_factor*n)))
    nlag = nfft*int(oversample)
    if nlag > 2**24:
        raise ValueError("Requested FFT exceeds memory safety limit")
    window = tukey(n, alpha=taper_fraction, sym=True)
    def removed_fraction(y):
        power = np.sum(np.abs(y)**2)
        return float(1-np.sum(np.abs(y*window)**2)/power) if power > 0 else 1.0
    fa = np.fft.fft(a*window, nfft)*dt_s
    fb = np.fft.fft(b*window, nfft)*dt_s
    f = np.fft.fftfreq(nfft, dt_s)
    mask = (np.abs(f) >= fmin_hz) & (np.abs(f) <= fmax_hz)
    if signal_kind == "real_strain":
        mask &= f > 0
    if np.count_nonzero(mask) < 3:
        raise ValueError("Frequency band contains fewer than three FFT bins")
    ff = f[mask]
    if isinstance(psd, str):
        if psd == "white":
            sn, psd_name = np.ones(len(ff)), "white-unit-PSD"
        elif psd == "aligo_historical":
            sn, psd_name = aligo_historical_psd(np.abs(ff)), "LALAdvLIGOPsd-historical-analytic"
        else:
            raise ValueError("Unknown PSD name")
    elif callable(psd):
        sn, psd_name = np.asarray(psd(np.abs(ff)), float), getattr(psd, "__name__", "custom-callable")
        if sn.ndim == 0:
            sn = np.full(len(ff), sn)
    else:
        raise ValueError("PSD must be a supported name or callable on positive frequencies")
    if sn.shape != ff.shape or not np.all(np.isfinite(sn)) or np.any(sn <= 0):
        raise ValueError("PSD must return finite strictly positive values on every selected bin")
    # Common scale cancels from normalized matches and avoids extreme weights.
    weights = np.median(sn)/sn
    aa, bb = fa[mask], fb[mask]
    na, nb = np.sum(np.abs(aa)**2*weights), np.sum(np.abs(bb)**2*weights)
    if not np.isfinite(na+nb) or min(na, nb) <= np.finfo(float).tiny:
        raise ValueError("Zero or numerically invalid in-band norm")
    cross = aa*np.conj(bb)*weights/np.sqrt(na*nb)
    embedded = np.zeros(nlag, complex)
    # Frequency integer indices retain their signs when embedding for finer lags.
    bins = np.rint(ff*nfft*dt_s).astype(np.int64) % nlag
    embedded[bins] = cross
    corr = np.fft.ifft(embedded)*nlag
    lag = np.fft.fftfreq(nlag, d=1/(nfft*dt_s))
    allowed = np.abs(lag) <= max_time_shift_s+1e-14
    idx = np.flatnonzero(allowed)
    best = int(idx[np.argmax(np.abs(corr[idx]))])
    tau, z = float(lag[best]), complex(corr[best])
    step = dt_s/oversample

    def correlation_at(shift):
        return np.sum(cross*np.exp(2j*np.pi*ff*shift))

    # Refine the best grid peak only, retaining boundaries and the sampled maximum.
    if refine and max_time_shift_s > 0:
        left, right = max(-max_time_shift_s, tau-step), min(max_time_shift_s, tau+step)
        opt = minimize_scalar(lambda x: -abs(correlation_at(x)), bounds=(left, right),
                              method="bounded", options={"xatol": max(1e-15, dt_s*1e-9)})
        for candidate in [float(opt.x), -max_time_shift_s, max_time_shift_s]:
            value = correlation_at(candidate)
            if abs(value) > abs(z)+5e-15:
                tau, z = candidate, complex(value)
    raw_match = float(abs(z))
    if raw_match > 1+1e-10:
        raise ArithmeticError("Normalized overlap violates Cauchy-Schwarz")
    match = float(np.clip(raw_match, 0, 1))
    zero = complex(np.sum(cross))
    return {
        "schema": "blackhole-waveform-metric/1", "signal_kind": signal_kind,
        "match": match, "mismatch": 1-match,
        "time_shift_s": tau, "phase_shift_rad": float(np.angle(z)),
        "amplitude_scale_b_to_a": float(np.sqrt(na/nb)*abs(z)),
        "fixed_time_phase_match": float(zero.real),
        "phase_only_match": float(min(1, abs(zero))),
        "at_time_boundary": bool(max_time_shift_s > 0 and abs(tau) >= max_time_shift_s-step/2),
        "config": {"dt_s": float(dt_s), "samples": n, "record_duration_s": (n-1)*dt_s,
                   "fmin_hz": float(fmin_hz), "fmax_hz": float(fmax_hz),
                   "actual_abs_frequency_range_hz": [float(min(abs(ff))), float(max(abs(ff)))],
                   "selected_frequency_bins": len(ff), "nyquist_hz": nyquist,
                   "psd": psd_name, "psd_source": PSD_SOURCE if psd_name.startswith("LAL") else None,
                   "taper": "symmetric Tukey", "taper_fraction": float(taper_fraction),
                   "unweighted_power_removed_by_taper": {"a": removed_fraction(a), "b": removed_fraction(b)},
                   "detrend": "none", "nfft": nfft, "padding_factor": float(padding_factor),
                   "frequency_step_hz": 1/(nfft*dt_s), "lag_grid_step_s": step,
                   "max_time_shift_s": float(max_time_shift_s), "refine_best_grid_peak": bool(refine),
                   "finite_support": "Same finite record; taper then zero-pad; fixed full-record norms for all lags",
                   "phase_semantics": "global complex mode phase" if signal_kind == "complex_mode" else "positive-frequency Fourier quadrature phase; not general orbital phase",
                   "interpretation": "conditional finite-record shape comparison; not detector-network or population inference"},
    }
