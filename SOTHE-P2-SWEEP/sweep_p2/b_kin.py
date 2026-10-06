"""B1 (kinematic layer over the grid) and B9 (extraction-bias parity tables, Supplement Tables 2 and 3).

kappa_extract() is suite.kinematic_kappa() of SOTHE-P2 1.0.0 with the grid size as a parameter; at n = 4096 it is
the same computation (checked by the oracle O-B1-identity).  eq7pp() is Eq. (7'') of the manuscript."""
import math

import numpy as np

from sothe_p2 import compat as C
from sothe_p2 import suite

PUBLISHED_TABLE2 = {"odd": {4095: -0.333263, 8191: -0.333316, 16383: -0.333329, 32767: -0.333332, 65535: -0.333333},
                    "even": {4096: -0.583065, 8192: -0.583266, 16384: -0.583317, 32768: -0.583329, 65536: -0.583332}}
PUBLISHED_TABLE3 = {"extracted": {1.0: 7.7027e-5, 2.0: 6.1600e-4, 3.0: 2.0778e-3},
                    "eq7pp": {1.0: 7.7036e-5, 2.0: 6.1628e-4, 3.0: 2.0800e-3}}


def kappa_extract(N, d3, kg, drive, k_op, n=4096):
    """suite.kinematic_kappa with an n-point grid on [-20, 20]."""
    Sv = float(suite.slow_light(k_op, d3, kg))
    c_match = max(abs(float(suite.vg0(k_op, d3))) * Sv, 1e-3)
    kappa0 = N * Sv
    tau_h = -(drive - 1.0) / (2.0 * N)
    tau = C.mlinspace(-20, 20, n)
    v = suite.flow_profile(tau, c_match, kappa0, tau_h)
    kappa, x_h = suite.surface_gravity(tau, v, c_match)
    return kappa, {"Sv": Sv, "c_match": c_match, "kappa0": kappa0, "tau_h": tau_h, "x_h": x_h}


def eq7pp(N, d3, kg, drive, k_op, n=4096):
    """Eq. (7''): kappa_num / kappa0 - 1 = -lambda^2 [dtau^2/3 + s (dtau - s)], s = crossing offset above the node below."""
    Sv = float(suite.slow_light(k_op, d3, kg))
    c_match = max(abs(float(suite.vg0(k_op, d3))) * Sv, 1e-3)
    kappa0 = N * Sv
    lam = kappa0 / c_match
    dt = 40.0 / (n - 1)
    tau_h = -(drive - 1.0) / (2.0 * N)
    j = math.floor((tau_h + 20.0) / dt)
    s = (tau_h + 20.0) - j * dt
    return -lam ** 2 * (dt ** 2 / 3.0 + s * (dt - s)), lam, dt, s


def kin_chunk(N, KOP, D3, KG, drive, kg_flow, n=4096):
    """All B1 rows for one amplitude N: every (k_op, delta3, kappa_g) point, and the kappa_g flows."""
    rows, flows = [], []
    for kop in KOP:
        for d3 in D3:
            for kg in KG:
                k1, aux = kappa_extract(N, d3, kg, 1.0, kop, n)
                k0 = aux["kappa0"]
                pred, lam, dt, s = eq7pp(N, d3, kg, 1.0, kop, n)
                kd = np.array([kappa_extract(N, d3, kg, d, kop, n)[0] for d in drive])
                mirror = [abs(kappa_extract(N, d3, kg, d, kop, n)[0] - kappa_extract(N, d3, kg, 2.0 - d, kop, n)[0]) for d in drive]
                preds = [eq7pp(N, d3, kg, d, kop, n)[0] for d in drive]
                rows.append({
                    "N": N, "k_op": kop, "d3": d3, "kg": kg, "S": aux["Sv"], "c_match": aux["c_match"], "two_c": 2 * aux["c_match"],
                    "lambda": lam, "kappa_closed": k0, "kappa_num": k1, "deficit_rel": k1 / k0 - 1.0 if k0 > 0 else float("nan"),
                    "eq7pp_rel": pred, "eq7pp_dev": (k1 / k0 - 1.0) - pred if k0 > 0 else float("nan"),
                    "drive": list(drive), "kappa_drive": kd.tolist(), "drive_spread_abs": float(np.max(kd) - np.min(kd)),
                    "drive_spread_rel": float((np.max(kd) - np.min(kd)) / k0) if k0 > 0 else float("nan"),
                    "band_rel": 0.25 * (lam * dt) ** 2, "band_abs": 0.25 * (lam * dt) ** 2 * k0,
                    "eq7pp_drive_max_dev": float(np.max(np.abs(kd / k0 - 1.0 - np.array(preds)))) if k0 > 0 else float("nan"),
                    "mirror_max": float(np.max(mirror)), "T_H": k0 / (2 * math.pi), "HP": 2 * math.pi / k0 if k0 > 0 else float("inf"),
                    "mu": 0.5 * N ** 2})
            fl = [kappa_extract(N, d3, g, 1.0, kop, n)[0] for g in kg_flow]
            cl = [N * float(suite.slow_light(kop, d3, g)) for g in kg_flow]
            flows.append({"N": N, "k_op": kop, "d3": d3, "kg": list(kg_flow), "kappa_num": fl, "kappa_closed": cl,
                          "monotone": bool(np.all(np.diff(fl) < 0)) if len(fl) > 1 else True})
    ident = None
    if n == 4096:                                      # identity with the SOTHE-P2 function at the operating point
        a = suite.kinematic_kappa(2.0, 0.05, 0.3, 1.0, 1.0)[0]
        b = kappa_extract(2.0, 0.05, 0.3, 1.0, 1.0, 4096)[0]
        ident = {"suite": a, "here": b, "equal": a == b}
    return {"N": N, "rows": rows, "flows": flows, "identity": ident}


def parity(cfg_par, op):
    """Supplement Table 2 (deficit / (lambda dtau)^2 at d = 1 for odd and even n_tau) and Table 3 (N = 1, 2, 3 at kappa_g = 0)."""
    N, d3, kg, kop = op["N"], op["d3"], op["kg"], op["k_op"]
    t2 = []
    for par, ns in (("odd", cfg_par["odd"]), ("even", cfg_par["even"])):
        for n in ns:
            k, aux = kappa_extract(N, d3, kg, 1.0, kop, n)
            pred, lam, dt, s = eq7pp(N, d3, kg, 1.0, kop, n)
            val = (k / aux["kappa0"] - 1.0) / (lam * dt) ** 2
            pub = PUBLISHED_TABLE2[par].get(n)
            t2.append({"parity": par, "n_tau": n, "deficit_over_lamdt2": val, "eq7pp_over_lamdt2": pred / (lam * dt) ** 2,
                       "s_over_dt": s / dt, "published": pub, "dev_published": (val - pub) if pub is not None else None})
    t3 = []
    for Nc in cfg_par["cont_N"]:
        k, aux = kappa_extract(Nc, d3, 0.0, 1.0, kop, 4096)
        pred, lam, dt, s = eq7pp(Nc, d3, 0.0, 1.0, kop, 4096)
        ext = Nc - k
        eqv = -pred * Nc
        t3.append({"N": Nc, "extracted": ext, "eq7pp": eqv, "relative": abs(ext - eqv) / eqv,
                   "published_extracted": PUBLISHED_TABLE3["extracted"].get(Nc),
                   "published_eq7pp": PUBLISHED_TABLE3["eq7pp"].get(Nc)})
    return {"table2": t2, "table3": t3}
