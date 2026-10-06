"""The run catalogue of SOTHE-P2-RAD 1.0.0.

Groups
  table     the radiative-loss protocol of SOTHE-P2 1.0.0 (launched N sech(N tau), laboratory frame, xi in [0, 40],
            dxi = 0.002, n = 8192, L = 400): the points of the loss table of the article, reproduced
  prepared  N = 1, delta3 raised adiabatically over xi_ramp = 100 (raised cosine) and then held; computational frame
            moving at -delta3 (the drift of the zero-momentum soliton); steady emission rate from the radiation behind
            the soliton
  recoil    N = 1, launched as in 'table' but followed to xi = 160 in the frame moving at -delta3
  scale     N = 2 launched runs to xi = 40, the scaled images of the N = 1 recoil runs (psi = N Psi(N tau, N^2 xi))
  check     numerical checks: dxi/2, 2n, a wider box, other ramp lengths, a spectral filter, the laboratory frame
"""

TABLE_POINTS = [(1.0, 0.02), (1.0, 0.03), (1.0, 0.04), (1.0, 0.05), (1.0, 0.06), (1.0, 0.08), (1.0, 0.10), (2.0, 0.025),
                (2.0, 0.05), (1.6, 0.05)]
PREPARED_STEADY = [0.035, 0.04, 0.045, 0.05, 0.055, 0.06, 0.065, 0.07, 0.08, 0.09]
PREPARED_STRONG = [0.10, 0.11, 0.12]
RECOIL = [0.06, 0.07, 0.08, 0.09, 0.10, 0.11, 0.12]
SCALE = [(2.0, 0.04), (2.0, 0.06)]


def _id(prefix, N, d3, extra=""):
    return "%s_N%s_d%s%s" % (prefix, ("%g" % N).replace(".", "p"), ("%.3f" % d3).replace(".", "p"), extra)


def catalogue():
    runs = []
    for N, d3 in TABLE_POINTS:
        runs.append({"id": _id("table", N, d3), "group": "table", "N": N, "d3": d3, "xi_end": 40.0, "sample": 0.5, "V": 0.0})
    for d3 in PREPARED_STEADY:
        runs.append({"id": _id("prep", 1.0, d3), "group": "prepared", "N": 1.0, "d3": d3, "schedule": "ramp", "xi_ramp": 100.0,
                     "xi_end": 200.0, "V": -d3})
    for d3 in PREPARED_STRONG:
        runs.append({"id": _id("prep", 1.0, d3), "group": "prepared", "N": 1.0, "d3": d3, "schedule": "ramp", "xi_ramp": 100.0,
                     "xi_end": 250.0, "V": -d3})
    for d3 in RECOIL:
        runs.append({"id": _id("recoil", 1.0, d3), "group": "recoil", "N": 1.0, "d3": d3, "xi_end": 160.0, "V": -d3})
    for N, d3 in SCALE:
        runs.append({"id": _id("scale", N, d3), "group": "scale", "N": N, "d3": d3, "xi_end": 40.0, "V": -d3 * N * N})
    chk = [
        ("prep", 0.05, {"dxi": 0.001}, "_dxi_half"),
        ("prep", 0.05, {"n": 16384}, "_n_double"),
        ("prep", 0.05, {"n": 12288, "L": 600.0}, "_L600"),
        ("prep", 0.05, {"xi_ramp": 50.0}, "_ramp50"),
        ("prep", 0.05, {"xi_ramp": 200.0, "xi_end": 300.0}, "_ramp200"),
        ("prep", 0.05, {"kfilter": 0.8}, "_kfilter"),
        ("prep", 0.04, {"dxi": 0.001}, "_dxi_half"),
        ("prep", 0.04, {"n": 16384}, "_n_double"),
        ("recoil", 0.10, {"dxi": 0.001}, "_dxi_half"),
        ("recoil", 0.10, {"n": 16384}, "_n_double"),
        ("recoil", 0.10, {"V": 0.0}, "_labframe"),
    ]
    for kind, d3, over, suffix in chk:
        if kind == "prep":
            r = {"id": _id("check_prep", 1.0, d3, suffix), "group": "check", "N": 1.0, "d3": d3, "schedule": "ramp", "xi_ramp": 100.0,
                 "xi_end": 200.0, "V": -d3}
        else:
            r = {"id": _id("check_recoil", 1.0, d3, suffix), "group": "check", "N": 1.0, "d3": d3, "xi_end": 160.0, "V": -d3}
        r.update(over)
        runs.append(r)
    ids = [r["id"] for r in runs]
    if len(set(ids)) != len(ids):
        raise RuntimeError("duplicate run identifiers")
    return runs


def select(runs, sets):
    if not sets or "all" in sets:
        return list(runs)
    return [r for r in runs if r["group"] in sets or r["id"] in sets]


def cost_estimate(r):
    """Relative cost: steps times grid size (for ordering the work, longest first)."""
    dxi = r.get("dxi", 0.002)
    n = r.get("n", 8192)
    return r["xi_end"] / dxi * n
