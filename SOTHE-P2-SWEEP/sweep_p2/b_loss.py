"""B3: radiative loss and drift of the launched soliton psi0 = N sech(N tau) under Eq. (1).

One propagation per task (kit.splitstep_loss of SOTHE-P2 1.0.0, sub-grid peak).  reduce_point() is the per-point
reduction of sothe_p2/g4c.py run() (window rates and drifts, c = mean(rate * xi_mid) over the last three windows,
loss over [0, 40], loss over one Hawking period after xi = 10, recoil check), unchanged, so that the points of
job 530047 are reproduced to round-off (oracle O-B3-*)."""
import math

import numpy as np

from sothe_p2 import compat as C
from sothe_p2 import kit

WINDOWS = [(0, 5), (5, 10), (10, 20), (20, 30), (30, 40)]
KG_HP, KOP_HP = 0.3, 1.0                     # g4c: kappa = N S(k_op = 1) at kappa_g = 0.3 for the Hawking period


def propagate(N, d3, kind, P):
    """kind: base | half (dxi/2) | n2 (2n).  Returns the series as plain lists."""
    dxi, n = P["dxi"], P["n"]
    if kind == "half":
        dxi = dxi / 2
    elif kind == "n2":
        n = 2 * n
    X = kit.splitstep_loss(N, d3, P["xi_end"], dxi, n, P["L"], 0.0)
    return {"N": N, "d3": d3, "kind": kind, "dxi": dxi, "n": n, "L": P["L"], "xi_end": P["xi_end"],
            "xi": X["xi"].tolist(), "Q": X["Q"].tolist(), "tau_c": X["tau_c"].tolist(), "k_s": X["k_s"].tolist(),
            "tau_c_grid": X["tau_c_grid"].tolist(), "loss_rate": X["loss_rate"], "drift": X["drift"],
            "drift_grid": X["drift_grid"], "Q_ratio": X["Q_ratio"]}


def _polyfit1(x, y):
    p = np.polyfit(x, y, 1)
    return float(p[0]), float(p[1])


def reduce_point(base, half=None, n2=None):
    """g4c.run per-point reduction (same code path, same definitions)."""
    N, d3 = base["N"], base["d3"]
    xi, Q, tc, tg, ks = (np.asarray(base[k], dtype=float) for k in ("xi", "Q", "tau_c", "tau_c_grid", "k_s"))
    rec = {"N": N, "d3": d3, "Nd3": N * d3, "loss_rate": base["loss_rate"], "drift": base["drift"], "Q_ratio": base["Q_ratio"]}
    beta0 = 0.5 - d3
    S = abs(beta0) / math.sqrt(beta0 ** 2 + KG_HP ** 2)
    rec["kappa_hp"] = N * S
    HP = 2 * math.pi / rec["kappa_hp"]
    rec["hawking_period"] = HP
    wr, wd, wg, wm, wk = (np.zeros(len(WINDOWS)) for _ in range(5))
    for iw, (lo, hi) in enumerate(WINDOWS):
        sw = (xi >= lo) & (xi <= hi)
        wr[iw] = -_polyfit1(xi[sw], np.log(Q[sw]))[0]
        wd[iw] = _polyfit1(xi[sw], tc[sw])[0]
        wg[iw] = _polyfit1(xi[sw], tg[sw])[0]
        wk[iw] = _polyfit1(xi[sw], ks[sw])[0]
        wm[iw] = 0.5 * (lo + hi)
    rec["window_lo"] = [w[0] for w in WINDOWS]
    rec["window_hi"] = [w[1] for w in WINDOWS]
    rec["window_rate"] = wr.tolist()
    rec["window_drift"] = wd.tolist()
    rec["window_drift_grid"] = wg.tolist()
    rec["window_ks_rate"] = wk.tolist()
    rec["c_over_xi"] = float(np.mean(wr[2:5] * wm[2:5]))
    rec["loss_0_40"] = float(1 - Q[-1] / Q[0])
    rec["loss_first_hp"] = float(1 - C.interp1(xi, Q, 10.0 + HP) / C.interp1(xi, Q, 10.0))
    Neff = Q / 2
    vpred = (ks - 3 * d3 * ks ** 2) - d3 * Neff ** 2
    late = (xi >= 10) & (xi <= 38)
    rec["k_s_early"] = float(np.mean(ks[xi <= 5]))
    rec["k_s_late"] = float(np.mean(ks[xi >= 30]))
    for key, arr, suffix in (("tau_c", tc, ""), ("tau_c_grid", tg, "_grid")):
        dsm = C.conv_same(C.mgradient(arr, xi), np.ones(9) / 9)
        rec["recoil_max_dev" + suffix] = float(np.max(np.abs(dsm[late] - vpred[late])))
        rec["recoil_rel_dev" + suffix] = rec["recoil_max_dev" + suffix] / float(np.max(np.abs(dsm[late])))
    rec["drift_leading_order"] = -d3 * N ** 2
    if half is not None and n2 is not None:
        Qh0, Qhe = half["Q"][0], half["Q"][-1]
        Qn0, Qne = n2["Q"][0], n2["Q"][-1]
        rec["loss_rate_dxi_half"] = half["loss_rate"]
        rec["loss_rate_n2"] = n2["loss_rate"]
        rec["rel_dev"] = max(abs(half["loss_rate"] - base["loss_rate"]), abs(n2["loss_rate"] - base["loss_rate"])) / abs(base["loss_rate"])
        rec["loss_0_40_dxi_half"] = 1 - Qhe / Qh0
        rec["loss_0_40_n2"] = 1 - Qne / Qn0
        rec["rel_dev_0_40"] = max(abs(rec["loss_0_40_dxi_half"] - rec["loss_0_40"]), abs(rec["loss_0_40_n2"] - rec["loss_0_40"])) / rec["loss_0_40"]
        # drift convergence (sub-grid peak), window by window
        dh = [reduce_windows(half)[i] for i in range(len(WINDOWS))]
        dn = [reduce_windows(n2)[i] for i in range(len(WINDOWS))]
        rec["drift_conv_max_abs"] = float(max(max(abs(a - b) for a, b in zip(dh, wd)), max(abs(a - b) for a, b in zip(dn, wd))))
        rec["converged"] = rec["rel_dev"] <= 0.10
    else:
        rec["rel_dev"] = rec["rel_dev_0_40"] = rec["drift_conv_max_abs"] = float("nan")
        rec["converged"] = None
    return rec


def reduce_windows(X):
    xi, tc = np.asarray(X["xi"], dtype=float), np.asarray(X["tau_c"], dtype=float)
    out = []
    for lo, hi in WINDOWS:
        sw = (xi >= lo) & (xi <= hi)
        out.append(_polyfit1(xi[sw], tc[sw])[0])
    return out
