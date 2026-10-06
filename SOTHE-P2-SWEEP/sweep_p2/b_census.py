"""B4 (comoving channel census with measured drift epochs; Hawking window) and B2 (entanglement layer, closed form).

census_rows() calls kit.channel_census of SOTHE-P2 1.0.0 (311-point band [0.05, 1.6], Eq. (11')) in the laboratory
frame and in the frame comoving with each measured drift epoch.  At delta3 = 0 the closed-form edge of kit is 0/0,
so the dispersionless edge V^2/2 - mu is used (the exact delta3 -> 0 limit of Eq. (11')).

Hawking window (REVIEW_G5 S1-1): upstream escape needs omega > mu = N^2/2 and the downstream pair channel needs
omega < omega_max_pair, so emission inside the band needs [max(mu, w_lo), min(omega_max_pair, w_hi)] to be non-empty.
At the threshold omega = mu, Eq. (11) with kappa = N S(k_op) gives |beta_mu|^2 = 1/(exp(pi N / S) - 1)."""
import math

import numpy as np

from sothe_p2 import kit, suite

W_LO, W_HI = 0.05, 1.6


def kappa_closed(N, k_op, d3, kg):
    return N * float(suite.slow_light(k_op, d3, kg))


def beta2_thermal(w, kappa):
    return 1.0 / math.expm1(2 * math.pi * w / kappa)


def EN_thermal(w, kappa):
    return 2.0 * math.asinh(math.sqrt(beta2_thermal(w, kappa)))


def window(N, k_op, d3, kg, wmax):
    mu = 0.5 * N ** 2
    kap = kappa_closed(N, k_op, d3, kg)
    lo, hi = max(mu, W_LO), min(wmax, W_HI)
    ok = (hi > lo) and math.isfinite(wmax)
    n_mu = beta2_thermal(mu, kap) if kap > 0 else float("nan")
    out = {"mu": mu, "kappa": kap, "window_lo": lo if ok else None, "window_hi": hi if ok else None,
           "window_frac": (hi - lo) / (W_HI - W_LO) if ok else 0.0, "piN_over_S": math.pi * N / (kap / N) if kap > 0 else float("inf"),
           "n_at_mu": n_mu, "EN_at_mu": 2 * math.asinh(math.sqrt(n_mu)) if kap > 0 else float("nan")}
    if ok:
        out["EN_window_lo"] = EN_thermal(lo, kap)
        out["EN_window_hi"] = EN_thermal(hi, kap)
    return out


def census_one(N, k_op, d3, Vd, kg, nw):
    c = kit.channel_census(N, k_op, d3, Vd, kg, band=(W_LO, W_HI), nw=nw)
    wmax = c["omega_max_pair"] if d3 > 0 else c["omega_max_pair_no_tod"]
    c["omega_max_pair_used"] = wmax
    c.update(window(N, k_op, d3, kg, wmax))
    return c


def census_rows(N, d3, KOP, KG, epochs, nw):
    """epochs: list of (label, Vd, window) with the laboratory frame first ('lab', 0.0, None)."""
    rows = []
    for kop in KOP:
        for kg in KG:
            for lab, Vd, win in epochs:
                c = census_one(N, kop, d3, Vd, kg, nw)
                c["frame"] = lab
                c["epoch_window"] = win
                rows.append(c)
    return rows


def entanglement_layer(kappa, omega, eta=(1.0, 0.9, 0.5, 0.1)):
    """Closed forms of Sec. III D at kappa over the band: E_N, nu_-, nbar_max, and the loss-robust E_N under symmetric
    pure loss eta, nu_-' = 1 - eta (1 - exp(-2 r))."""
    w = np.asarray(omega, dtype=float)
    b2 = 1.0 / np.expm1(2 * np.pi * w / kappa)
    EN = 2 * np.arcsinh(np.sqrt(b2))
    nu = np.exp(-EN)
    out = {"omega": w.tolist(), "EN": EN.tolist(), "nu": nu.tolist(), "nbar_max": (0.5 * (np.exp(EN) - 1)).tolist(),
           "EN_max": float(EN[0]), "nu_max": float(nu[-1]), "IR_asymptote": math.log(2 * kappa / (math.pi * w[0])),
           "EN_at_3kappa": 2 * math.asinh(math.sqrt(1.0 / math.expm1(6 * math.pi)))}
    for e in eta:
        nup = 1 - e * (1 - nu)
        out["EN_loss_eta%.1f" % e] = (-np.log(nup)).tolist()
    return out
