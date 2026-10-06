"""
gnlse_dimensionless.py
======================

Python port of the dimensionless GNLSE solver and entropy machinery from the
SOTHE data and code package.

This is a *reference implementation* meant to validate the MATLAB code in
``../src/`` and to make the formalism reproducible without a MATLAB licence.
It uses only ``numpy`` and ``scipy``.

The solved equation is Eq. (9) of the paper:

    d psi / d xi = i * ( (1/2) d^2 psi / d tau^2
                       + |psi|^2 psi
                       + i * delta3 * d^3 psi / d tau^3 )

with the *integrability-consistent anomalous-dispersion sign* (the +i delta3
TOD term, NOT the -i delta3 of the original 2025 submission).

With numpy's ``fft`` convention exp(-i k tau), partial_tau <-> +i k, so
partial_tau^2 <-> -k^2 and partial_tau^3 <-> -i k^3. The spectral linear
operator is therefore

    D(k) = i * [ (1/2)(-k^2) + i delta3 (-i k^3) ]
         = i * [ -(1/2) k^2 + delta3 k^3 ]
         = -(i/2) k^2 + i delta3 k^3 .

Reference
---------
H. Oguz, "Generalized Thermodynamics of Solitonic Event Horizons in
Dispersive Field Theories," Classical and Quantum Gravity (2026).

Part of the SOTHE data and code package. Author: H. Oguz, 2026. MIT licence.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import Callable


# -----------------------------------------------------------------------------
# Core library: spectral entropy and eta_GSL
# -----------------------------------------------------------------------------

def spectral_entropy(k: np.ndarray, P: np.ndarray, mask: np.ndarray) -> float:
    """Differential Shannon entropy of a masked, normalized spectral PDF.

    S = - integral_{k in support(mask)} p(k) ln p(k) dk,
    p(k) = (P(k) * mask(k)) / integral( P * mask, dk ).
    Returns S in nats.
    """
    P_sub = np.maximum(P * mask, 0.0)
    Z = np.trapezoid(P_sub, k)
    if Z < 1e-20 * np.max(P + 1e-30):
        return 0.0
    p = P_sub / Z
    valid = p > 1e-20
    return float(-np.trapezoid(p[valid] * np.log(p[valid]), k[valid]))


def eta_gsl(xi: np.ndarray, S_tot: np.ndarray) -> np.ndarray:
    """Running entropy-production efficiency, Eq. (8) of the paper.

    eta(xi) = integral_0^xi max(0, dS/dxi') dxi' / integral_0^xi |dS/dxi'| dxi'.

    Returns the same-length array; eta[0] = 0 by convention.
    """
    xi = np.asarray(xi, dtype=float).ravel()
    S_tot = np.asarray(S_tot, dtype=float).ravel()
    if xi.shape != S_tot.shape:
        raise ValueError("xi and S_tot must have the same shape.")
    if xi.size < 3:
        raise ValueError("Need at least 3 samples to compute eta_GSL.")

    # Centered finite differences with one-sided endpoints
    dS = np.empty_like(S_tot)
    dS[1:-1] = (S_tot[2:] - S_tot[:-2]) / (xi[2:] - xi[:-2])
    dS[0]    = (S_tot[1] - S_tot[0])     / (xi[1] - xi[0])
    dS[-1]   = (S_tot[-1] - S_tot[-2])   / (xi[-1] - xi[-2])

    pos    = np.maximum(dS, 0.0)
    absd   = np.abs(dS)
    # cumtrapz-equivalent
    num    = _cumtrapz(pos,  xi)
    denom  = _cumtrapz(absd, xi)
    eta    = np.zeros_like(xi)
    nz     = denom > 0
    eta[nz] = num[nz] / denom[nz]
    return eta


def _cumtrapz(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    """Cumulative trapezoidal integral with leading zero (MATLAB cumtrapz)."""
    out = np.zeros_like(y, dtype=float)
    out[1:] = np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))
    return out


def soliton_mask(k: np.ndarray, P: np.ndarray,
                 width: float, order: int = 10):
    """Dynamically centred super-Gaussian system/bath spectral mask.

    Returns (mask_S, mask_R, k_c) where k_c is the first moment of P on k.
    """
    Ptot = np.trapezoid(P, k)
    k_c = 0.0 if Ptot < 1e-20 else float(np.trapezoid(k * P, k) / Ptot)
    mS = np.exp(-((k - k_c) / width) ** order)
    return mS, 1.0 - mS, k_c


def find_rr_lobe(k: np.ndarray, P: np.ndarray, k_c: float,
                 mask_width: float, side: str = "positive",
                 smoothing_sigma: float = 3.0,
                 search_radius_frac: float = 0.15) -> float:
    """Robust locator for the resonant-radiation (Cherenkov) lobe centre.

    Implements the "automatically located lobe fit on the log-spectral envelope"
    extraction procedure reported in Sec. IV.D of the paper.

    Strategy:
      1. Restrict to the dispersive-wave side, well outside the soliton mask
         (k > k_c + 3 * mask_width  for the +delta3 anomalous-dispersion
         convention).
      2. Gaussian-smooth the log power spectrum to suppress grid-scale jitter.
      3. Find the argmax k_0 of the smoothed log-spectrum on the side region.
      4. Return the first moment (centroid) of the surrounding lobe within
         a fractional radius search_radius_frac * |k_0|.

    This is one of three independent extraction procedures used in the paper
    (the others are a fixed-window tail fit and a comoving-frame fit). The
    paper reports c_RR = 1.01, R^2 = 0.93, and invariant k_RR * delta3 = 1.02
    on delta3 in [0.06, 0.10].
    """
    try:
        from scipy.ndimage import gaussian_filter1d
    except ImportError as e:
        raise ImportError("find_rr_lobe needs scipy") from e

    if side == "positive":
        sel = (k > k_c + 3.0 * mask_width)
    elif side == "negative":
        sel = (k < k_c - 3.0 * mask_width)
    else:
        raise ValueError("side must be 'positive' or 'negative'")
    sel &= (P > P.max() * 1e-6)
    if not np.any(sel):
        return float("nan")

    logP = np.log(P + 1e-30)
    logP_smooth = gaussian_filter1d(logP, sigma=smoothing_sigma)
    P_search = np.where(sel, logP_smooth, -np.inf)
    j0 = int(np.argmax(P_search))
    k0 = k[j0]

    r = search_radius_frac * abs(k0)
    lobe = (np.abs(k - k0) < r) & sel
    if not np.any(lobe):
        return float(k0)
    return float(np.trapezoid(k[lobe] * P[lobe], k[lobe])
                 / np.trapezoid(P[lobe], k[lobe]))


# -----------------------------------------------------------------------------
# Solver
# -----------------------------------------------------------------------------

@dataclass
class GNLSEResult:
    xi: np.ndarray
    tau: np.ndarray
    k: np.ndarray
    psi_hist: np.ndarray
    spec_hist: np.ndarray
    S_hor: np.ndarray
    S_rad: np.ndarray
    S_tot: np.ndarray
    eta_gsl: np.ndarray
    photon_number: np.ndarray
    k_c: np.ndarray
    Delta_S_tot: float = field(init=False)
    eta_GSL_final: float = field(init=False)
    params: dict = field(default_factory=dict)

    def __post_init__(self):
        self.Delta_S_tot   = float(self.S_tot[-1] - self.S_tot[0])
        self.eta_GSL_final = float(self.eta_gsl[-1])


def _initial_pulse(tau: np.ndarray, shape: str, N_sol: float) -> np.ndarray:
    if shape == "sech":
        return N_sol * (1.0 / np.cosh(tau))
    if shape == "gaussian":
        return N_sol * np.exp(-tau ** 2 / 2)
    if shape == "super_gaussian":
        return N_sol * np.exp(-((tau ** 2 / 2) ** 2))
    if shape == "chirped":
        chirp = 0.5
        return N_sol * (1.0 / np.cosh(tau)) * np.exp(-1j * chirp * tau ** 2 / 2)
    raise ValueError(f"Unknown pulse shape {shape!r}")


def gnlse_dimensionless(N_sol: float = 3.5,
                        delta3: float = 0.02,
                        xi_max: float = 12.0,
                        n_steps: int = 6000,
                        Nt: int = 2 ** 14,
                        tau_window: float = 20.0,
                        n_save: int = 200,
                        pulse_shape: str = "sech",
                        mask_width: float = 3.0,
                        mask_order: int = 10,
                        verbose: bool = True) -> GNLSEResult:
    """Symmetric SSFM integration of the dimensionless GNLSE [Eq. (9)]."""

    # Grid
    dtau = tau_window / Nt
    tau  = (np.arange(Nt) - Nt // 2) * dtau
    k_natural  = (np.arange(Nt) - Nt // 2) * (2 * np.pi / tau_window)
    k_natural  = np.fft.fftshift(k_natural)              # math-natural order
    k_disp     = np.fft.fftshift(k_natural)              # monotone for diagnostics

    psi   = _initial_pulse(tau, pulse_shape, N_sol)
    psi_f = np.fft.fft(psi)

    # Integrability-consistent dispersion operator
    D = -1j * 0.5 * k_natural ** 2 + 1j * delta3 * k_natural ** 3

    dxi       = xi_max / n_steps
    half_step = np.exp(D * dxi / 2)

    save_id = np.unique(np.round(np.linspace(1, n_steps, n_save)).astype(int))
    n_save  = len(save_id)
    xi_vec  = np.linspace(0, xi_max, n_save)

    psi_hist  = np.zeros((Nt, n_save), dtype=complex)
    spec_hist = np.zeros((Nt, n_save))
    S_hor     = np.zeros(n_save)
    S_rad     = np.zeros(n_save)
    photon_N  = np.zeros(n_save)
    k_c_hist  = np.zeros(n_save)

    if verbose:
        print(f"gnlse_dimensionless (py): N={N_sol}, delta3={delta3}, "
              f"xi_max={xi_max}, Nt={Nt}")

    c = 0
    next_prog = 0.1
    for n in range(1, n_steps + 1):
        # SSFM step
        psi_f = psi_f * half_step
        psi_t = np.fft.ifft(psi_f)
        psi_t = psi_t * np.exp(1j * np.abs(psi_t) ** 2 * dxi)
        psi_f = np.fft.fft(psi_t)
        psi_f = psi_f * half_step

        if c < n_save and n == save_id[c]:
            psi_t_save = np.fft.ifft(psi_f)
            psi_hist[:, c]  = psi_t_save
            P = np.abs(np.fft.fftshift(psi_f)) ** 2
            spec_hist[:, c] = P

            mS, mR, kc = soliton_mask(k_disp, P, mask_width, mask_order)
            S_hor[c] = spectral_entropy(k_disp, P, mS)
            S_rad[c] = spectral_entropy(k_disp, P, mR)
            photon_N[c] = float(np.trapezoid(np.abs(psi_t_save) ** 2, tau))
            k_c_hist[c] = kc
            c += 1

        if verbose and n / n_steps >= next_prog:
            print(f"  ... {round(100 * n / n_steps):3d}% done")
            next_prog += 0.1

    S_tot = S_hor + S_rad
    eta   = eta_gsl(xi_vec, S_tot)

    res = GNLSEResult(
        xi=xi_vec, tau=tau, k=k_disp,
        psi_hist=psi_hist, spec_hist=spec_hist,
        S_hor=S_hor, S_rad=S_rad, S_tot=S_tot,
        eta_gsl=eta, photon_number=photon_N, k_c=k_c_hist,
        params=dict(N_sol=N_sol, delta3=delta3, xi_max=xi_max,
                    n_steps=n_steps, Nt=Nt, tau_window=tau_window,
                    n_save=n_save, pulse_shape=pulse_shape,
                    mask_width=mask_width, mask_order=mask_order),
    )
    if verbose:
        print(f"Done. Delta S_tot = {res.Delta_S_tot:+.3f} nats, "
              f"eta_GSL(xi_f) = {res.eta_GSL_final:.3f}")
        print(f"Photon-number drift: {(photon_N[-1]/photon_N[0]-1):.2e}")
    return res


def rr_phase_match_fit(delta3_vec: np.ndarray, k_RR_vec: np.ndarray) -> dict:
    """Through-origin Akhmediev-Karlsson fit, k_RR = c_RR / delta3."""
    d3 = np.asarray(delta3_vec, dtype=float).ravel()
    kR = np.asarray(k_RR_vec, dtype=float).ravel()
    inv_d3 = 1.0 / d3
    c_RR   = float(np.sum(kR * inv_d3) / np.sum(inv_d3 ** 2))
    y_hat  = c_RR * inv_d3
    R2     = float(1.0 - np.sum((kR - y_hat) ** 2) / np.sum((kR - kR.mean()) ** 2))
    inv    = kR * d3
    return dict(c_RR=c_RR, R2=R2,
                invariant=inv,
                invariant_mean=float(inv.mean()),
                invariant_spread=float((inv.max() - inv.min()) / inv.mean()),
                residuals=kR - y_hat)


if __name__ == "__main__":
    # Smoke test
    r = gnlse_dimensionless(N_sol=3.5, delta3=0.02, xi_max=8.0,
                            Nt=2 ** 12, n_steps=2000, n_save=100, verbose=True)
    print(f"\nSmoke test result: Delta S_tot = {r.Delta_S_tot:.3f}, "
          f"eta_GSL = {r.eta_GSL_final:.3f}")
