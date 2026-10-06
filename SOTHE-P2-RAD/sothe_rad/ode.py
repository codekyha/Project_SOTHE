"""Recoil dynamics of a resonantly radiating soliton: adiabatic two-variable model.

The soliton is followed through its norm Q and its mean wavenumber k_s.  It radiates norm at the relative rate R
into the resonant wavenumber k_r and conserves momentum,

    dQ/dxi = -R Q,        d(k_s Q)/dxi = k_r dQ/dxi   =>   dk_s/dxi = -(k_r - k_s) R,

and R depends on the state through the scaled third-order coefficient of the local soliton (model.py):

    D2 = 1 - 6 delta3 k_s,  N = Q / (2 sqrt(D2)),  delta_tilde = N delta3 / D2^(3/2),
    k_r - k_s = q_r(delta_tilde) N / sqrt(D2),      R = N^2 Rs(delta_tilde),

Rs the scaled rate function (measured on prepared solitons, analysis.py).  This is the 'norm' amplitude law.  The
'deplete' law instead lets the amplitude follow the norm, dN/dxi = -N R, and keeps D2 out of it (a simpler
variant, also evaluated).  Asymptotically R decays as 1/(K xi) with K = d(1/R)/dxi, so -d ln Q/dxi -> c/xi with
c = 1/K; K is evaluated along the solution.
"""
import math

import numpy as np

from . import model

try:
    from scipy.integrate import solve_ivp
except ImportError:      # pragma: no cover
    solve_ivp = None


def state_quantities(Q, ks, d3, root="drift", amp_law="norm", N_dep=None):
    D = 1.0 - 6.0 * d3 * ks
    N = Q / (2.0 * math.sqrt(D)) if amp_law == "norm" else N_dep
    dt = N * d3 / D ** 1.5
    q = model.q_res(dt, root)
    return D, N, dt, q


def rhs(xi, y, d3, rate_scaled, root="drift", amp_law="norm", scale=1.0):
    if amp_law == "norm":
        lnQ, ks = y
        Q = math.exp(lnQ)
        D, N, dt, q = state_quantities(Q, ks, d3, root, "norm")
        R = scale * N ** 2 * rate_scaled(dt)
        return [-R, -(N / math.sqrt(D)) * q * R]
    lnQ, ks, lnN = y
    N = math.exp(lnN)
    D = 1.0 - 6.0 * d3 * ks
    dt = N * d3 / D ** 1.5
    q = model.q_res(dt, root)
    R = scale * N ** 2 * rate_scaled(dt)
    return [-R, -(N / math.sqrt(D)) * q * R, -R]


def rate_now(Q, ks, d3, rate_scaled, root="drift", amp_law="norm", N=None, scale=1.0):
    D = 1.0 - 6.0 * d3 * ks
    if amp_law == "norm":
        N = Q / (2.0 * math.sqrt(D))
    dt = N * d3 / D ** 1.5
    return scale * N ** 2 * rate_scaled(dt)


def integrate(d3, Q0, ks0, xi0, xi_end, rate_scaled, root="drift", amp_law="norm", scale=1.0, N0=None, xi_eval=None):
    """Integrate the model from (Q0, ks0) at xi0 to xi_end.  Returns a dict of arrays on xi_eval (default: 0.25 grid)."""
    if solve_ivp is None:
        raise RuntimeError("SciPy is required for the recoil model")
    if xi_eval is None:
        xi_eval = np.arange(xi0, xi_end + 1e-9, 0.25)
    if amp_law == "norm":
        y0 = [math.log(Q0), ks0]
    else:
        if N0 is None:
            N0 = Q0 / (2.0 * math.sqrt(1.0 - 6.0 * d3 * ks0))
        y0 = [math.log(Q0), ks0, math.log(N0)]
    sol = solve_ivp(rhs, (xi0, xi_end), y0, method="DOP853", t_eval=xi_eval, rtol=1e-11, atol=1e-13,
                    args=(d3, rate_scaled, root, amp_law, scale))
    if not sol.success:
        raise RuntimeError(sol.message)
    Q = np.exp(sol.y[0])
    ks = sol.y[1]
    D = 1.0 - 6.0 * d3 * ks
    N = Q / (2.0 * np.sqrt(D)) if amp_law == "norm" else np.exp(sol.y[2])
    dt = N * d3 / D ** 1.5
    R = np.array([scale * Ni ** 2 * rate_scaled(d) for Ni, d in zip(N, dt)])
    V = model.soliton_velocity(N, ks, d3)
    invR = 1.0 / R
    K = np.gradient(invR, sol.t)
    with np.errstate(divide="ignore"):
        c_inst = 1.0 / K          # infinite where 1/R does not change at the precision of the solution
    return {"xi": sol.t, "Q": Q, "ks": ks, "N": N, "dt": dt, "R": R, "V": V, "K": K, "c_inst": c_inst}


def calibrate_scale(R_meas, Q, ks, d3, rate_scaled, root="drift", amp_law="norm", N=None):
    """Multiplier that makes the model rate equal to a measured rate in the state (Q, ks)."""
    return R_meas / rate_now(Q, ks, d3, rate_scaled, root, amp_law, N, 1.0)


def window_exponents(xi, Q, edges):
    """c over windows [a, b]: -(ln Q(b) - ln Q(a)) / (ln b - ln a)."""
    out = []
    for a, b in zip(edges[:-1], edges[1:]):
        qa, qb = np.interp(a, xi, Q), np.interp(b, xi, Q)
        out.append(-(math.log(qb) - math.log(qa)) / (math.log(b) - math.log(a)))
    return out
