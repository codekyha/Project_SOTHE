"""Dispersion relation, resonance condition, local-soliton scaling and first-order (Born) emission rate.

Equation (the convention of the article):

    d_xi psi = i ( 1/2 psi_tt + |psi|^2 psi + i delta3 psi_ttt ),

with linear dispersion beta0(k) = k^2/2 - delta3 k^3 (a plane wave exp(i k tau - i beta0 xi)) and group
velocity v_g0(k) = k - 3 delta3 k^2.

Local soliton.  A soliton of mean wavenumber k_s is, about its carrier, a soliton of the same equation with
second-order coefficient D2 = 1 - 6 delta3 k_s; rescaling tau by sqrt(D2) and by its amplitude N maps it onto the
zero-carrier soliton of amplitude 1 of the equation with

    delta_tilde = N delta3 / D2^(3/2),

and its norm is Q = 2 N sqrt(D2).  Its velocity (centre of mass, exact for the linear part) is
V = v_g0(k_s) - delta3 N^2 / D2, and the phase-matching condition of the resonant (Cherenkov) radiation, in the
scaled wavenumber q measured from the carrier (k = k_s + q N / sqrt(D2)), reads

    drift-corrected:  delta_tilde q^3 - q^2/2 - delta_tilde q - 1/2 = 0,
    no drift term:    delta_tilde q^3 - q^2/2 - 1/2 = 0.

First-order (Born) rate.  The radiation is driven by the third-order term acting on the sech soliton; its
Fourier amplitude at the phase-matched wavenumber, delta_tilde q^3 * pi sech(pi q / 2), radiates norm at the
rate |S(q_r)|^2 / |v_rel| (golden rule; v_rel the group velocity of the radiation relative to the soliton).
Relative to the soliton norm 2 (scaled units):

    R_Born(delta_tilde) = pi^2 delta_tilde^2 q_r^6 sech^2(pi q_r / 2) / (2 |q_r - 3 delta_tilde q_r^2 + delta_tilde|).

In physical units the relative rate is R = N^2 R(delta_tilde) per unit xi.  The exponent pi q_r is fixed by the
pole of the sech at tau = i pi / (2N); the prefactor of the first-order formula is not, and the measured rates
(see analysis.py) determine it.
"""
import math

import numpy as np


def beta0(k, d3):
    return 0.5 * np.asarray(k) ** 2 - d3 * np.asarray(k) ** 3


def vg0(k, d3):
    return np.asarray(k) - 3.0 * d3 * np.asarray(k) ** 2


def D2(ks, d3):
    return 1.0 - 6.0 * d3 * np.asarray(ks)


def delta_tilde(N, ks, d3):
    """Scaled third-order coefficient of a local soliton of amplitude N and mean wavenumber ks."""
    return N * d3 / D2(ks, d3) ** 1.5


def amplitude_from_norm(Q, ks, d3):
    """Amplitude N of the local sech soliton of norm Q = 2 N sqrt(D2)."""
    return Q / (2.0 * np.sqrt(D2(ks, d3)))


def soliton_velocity(N, ks, d3):
    """Centre-of-mass velocity of the local soliton: v_g0(ks) - d3 N^2 / D2 (laboratory frame)."""
    return vg0(ks, d3) - d3 * N ** 2 / D2(ks, d3)


def _largest_real_root(coeffs):
    r = np.roots(coeffs)
    r = np.real(r[np.abs(np.imag(r)) < 1e-9 * max(1.0, np.max(np.abs(r)))])
    if r.size == 0:
        raise ValueError("no real root")
    return float(np.max(r))


def q_res(dt, variant="drift"):
    """Scaled resonant wavenumber q_r(delta_tilde): largest real root of the phase-matching cubic.
    variant 'drift': dt q^3 - q^2/2 - dt q - 1/2 = 0 (soliton velocity -dt included);
    variant 'nodrift': dt q^3 - q^2/2 - 1/2 = 0."""
    if dt <= 0:
        raise ValueError("delta_tilde must be positive")
    if variant == "drift":
        return _largest_real_root([dt, -0.5, -dt, -0.5])
    if variant == "nodrift":
        return _largest_real_root([dt, -0.5, 0.0, -0.5])
    raise ValueError(variant)


def q_res_vec(dts, variant="drift"):
    return np.array([q_res(float(d), variant) for d in np.atleast_1d(dts)])


def v_rel_scaled(q, dt, variant="drift"):
    """Group velocity of the radiation relative to the soliton, scaled units."""
    return q - 3.0 * dt * q ** 2 + (dt if variant == "drift" else 0.0)


def rate_born_scaled(dt, variant="drift"):
    """First-order (Born) relative emission rate in scaled units (soliton of amplitude 1, norm 2)."""
    q = q_res(dt, variant)
    s = dt * q ** 3 * math.pi / math.cosh(0.5 * math.pi * q)
    return s * s / (2.0 * abs(v_rel_scaled(q, dt, variant)))


def rate_born(N, ks, d3, variant="drift"):
    """First-order relative rate in physical units: N^2 R_Born(delta_tilde)."""
    return float(N) ** 2 * rate_born_scaled(float(delta_tilde(N, ks, d3)), variant)


def k_res_lab(N, ks, d3, variant="drift"):
    """Resonant wavenumber in the laboratory frame: ks + q_r N / sqrt(D2)."""
    dt = float(delta_tilde(N, ks, d3))
    return float(ks) + q_res(dt, variant) * float(N) / math.sqrt(float(D2(ks, d3)))


def k_res_measured(phase_rate, V, d3):
    """Resonant wavenumber from two measured quantities: the phase rotation rate at the soliton centre,
    phi_c = d(arg psi(tau_c))/dxi, and the soliton velocity V.  The radiation is phase-matched where
    beta0(k) + phi_c - k V = 0, i.e. d3 k^3 - k^2/2 + V k - phi_c = 0 (largest real root)."""
    return _largest_real_root([d3, -0.5, V, -phase_rate])


def exponent_scaled(dt, variant="drift"):
    """pi q_r(delta_tilde): the exponent of the rate, R ~ exp(-pi q_r)."""
    return math.pi * q_res(dt, variant)


def universal_rate(dt, C, p, variant="drift"):
    """Two-parameter form of the scaled rate: C dt^p exp(-pi q_r(dt))."""
    return C * dt ** p * math.exp(-exponent_scaled(dt, variant))


def dlnrate_dlndt(fun, dt, h=1e-4):
    """d ln R / d ln dt of a scaled rate function, by a centred difference in ln dt."""
    a, b = math.log(dt) - h, math.log(dt) + h
    return (math.log(fun(math.exp(b))) - math.log(fun(math.exp(a)))) / (2 * h)
