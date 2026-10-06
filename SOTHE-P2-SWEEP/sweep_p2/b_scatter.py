"""B7: mode-conversion scattering of the BdG operator on the soliton-plus-flow background (delta3 = 0).

REVIEW_G5 S1-2 / D-06: the protocol, numbers and checks behind the manuscript's sentence on the scattering solve.
INV5: scattering-derived alpha/beta are not Paper-2 claims; every number of this block is diagnostic (RECORDED).

Operator (frequency domain, u e^{-i w xi} + conj(v) e^{+i conj(w) xi}):
  M = [[H + F, -P], [P, -H + F]],  H = -1/2 d_tau^2 + mu - 2 P,  P = psi0^2 = N^2 sech^2(N tau),  mu = N^2 / 2,
  F = -(i/2)(V d_tau + d_tau V)  (the symmetrized advection by the kinematic flow V(tau)),
  V(tau) = c [1 + tanh(kappa0 (tau - tau_h) / c)],  kappa0 = N S(k_op; 0, kappa_g),  c = |vg0(k_op; 0)| S = k_op S.
Discretization: uniform grid on [-L, L], second-order stencils; outside the grid the coefficients are constant
(V = 0, P = 0 on the left; V = 2c, P = 0 on the right) and the solution is a sum of the exact discrete exterior modes
z^j, the two roots z of the u- and of the v-quadratic on each side.  Admissible exterior modes: outgoing (discrete
group velocity away from the grid) or decaying; the others (incoming, growing) are set by the right-hand side.  The
square sparse system (2n interior unknowns + 4 admissible amplitudes) is factorized once per omega (SuperLU) and
solved for each incoming channel with unit flux.  S entries are flux-normalized amplitudes; Krein-signed flux balance
sum_out s |S|^2 - s_in = 0 (s = +1 for u, -1 for v) is the unitarity check.

Channel labels: side (L, R) + component (u: positive norm, v: negative norm).  Reductions:
  r_conf = |S[Rv <- Ru]|^2 / |S[Ru <- Ru]|^2   (FINDINGS.md: the conversion ratio on the supersonic side)
  n_H    = |S[Lu <- Rv]|^2                     (upstream positive-norm quanta from the incoming negative-norm mode)
against the thermal forms exp(-2 pi w / kappa0) and 1 / (exp(2 pi w / kappa0) - 1)."""
import math

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

from sothe_p2 import suite

PROP_TOL = 1e-9


def kinematics(N, k_op, kg):
    """kappa0, c, S at delta3 = 0 (the kinematic layer of the manuscript with delta3 = 0)."""
    S = float(suite.slow_light(k_op, 0.0, kg))
    c = max(abs(float(suite.vg0(k_op, 0.0))) * S, 1e-3)
    return N * S, c, S


def flow(tau, c, kappa0, tau_h=0.0):
    return c * (1.0 + np.tanh(kappa0 * (tau - tau_h) / c))


def exterior_modes(omega, V, mu, h):
    """Roots z of the discrete u- and v-quadratics at constant V (P = 0): list of dicts with comp, z, k, prop,
    vg (discrete group velocity, propagating modes only), s (Krein sign)."""
    out = []
    cu = [-0.5 / h ** 2 - 0.5j * V / h, 1.0 / h ** 2 + (mu - omega), -0.5 / h ** 2 + 0.5j * V / h]
    cv = [0.5 / h ** 2 - 0.5j * V / h, -1.0 / h ** 2 - (mu + omega), 0.5 / h ** 2 + 0.5j * V / h]
    for comp, cc, s in (("u", cu, 1), ("v", cv, -1)):
        for z in np.roots(cc):
            z = complex(z)
            k = -1j * np.log(z) / h
            prop = abs(abs(z) - 1.0) < PROP_TOL
            if prop:
                kr = float(np.real(k))
                vg = (math.sin(kr * h) / h if comp == "u" else -math.sin(kr * h) / h) + V * math.cos(kr * h)
            else:
                vg = float("nan")
            out.append({"comp": comp, "z": z, "k": complex(k), "prop": bool(prop), "vg": vg, "s": s})
    return out


def _classify(md, side):
    for m in md:
        if m["prop"]:
            m["adm"] = (m["vg"] < 0) if side == "L" else (m["vg"] > 0)
        else:
            m["adm"] = (abs(m["z"]) > 1) if side == "L" else (abs(m["z"]) < 1)
    return md


def solve(N, c, kappa0, omega, L=25.0, h=0.01, tau_h=0.0, psi_on=True, flow_on=True, Vconst=None):
    """S-matrix at one omega.  Returns (S, info): S[i_out, j_in] flux-normalized; info has the channel labels
    ('Lu', 'Rv', ...), k and vg of each channel, the solve residual and the Krein flux balance per incoming channel."""
    mu = 0.5 * N ** 2
    n = int(round(2 * L / h)) + 1
    tau = -L + h * np.arange(n)
    if Vconst is not None:
        V = np.full(n, float(Vconst))
        VL = VR = float(Vconst)
    elif flow_on:
        V = flow(tau, c, kappa0, tau_h)
        VL, VR = 0.0, 2 * c
    else:
        V = np.zeros(n)
        VL = VR = 0.0
    P = (N / np.cosh(N * tau)) ** 2 if psi_on else np.zeros(n)
    V[0], V[-1] = VL, VR
    P[0], P[-1] = 0.0, 0.0
    mdL = _classify(exterior_modes(omega, VL, mu, h), "L")
    mdR = _classify(exterior_modes(omega, VR, mu, h), "R")
    allm = [("L", m) for m in mdL] + [("R", m) for m in mdR]
    na = len(allm)
    Vext = np.concatenate([[VL], V, [VR]])
    j = np.arange(n)
    Vj, Vm, Vp = Vext[1:n + 1], Vext[0:n], Vext[2:n + 2]
    fm = -0.5j * (-Vj - Vm) / (2 * h)
    fp = -0.5j * (Vj + Vp) / (2 * h)
    rows, cols, vals = [], [], []
    ghostL, ghostR = [], []                       # (row, coefficient, comp) of the stencil entries outside the grid
    for comp, off, sg in (("u", 0, 1.0), ("v", n, -1.0)):
        cm = sg * (-0.5 / h ** 2) + fm
        cp = sg * (-0.5 / h ** 2) + fp
        c0 = sg * (1.0 / h ** 2 + mu - 2 * P) - omega
        rows += list(off + j)
        cols += list(off + j)
        vals += list(c0)
        rows += list(off + j[1:])
        cols += list(off + j[:-1])
        vals += list(cm[1:])
        rows += list(off + j[:-1])
        cols += list(off + j[1:])
        vals += list(cp[:-1])
        ghostL.append((off + 0, cm[0], comp))
        ghostR.append((off + n - 1, cp[-1], comp))
        rows += list(off + j)                      # coupling: u-row -P v, v-row +P u
        cols += list((n if comp == "u" else 0) + j)
        vals += list(-P if comp == "u" else P)
    acol = {ia: 2 * n + ia for ia in range(na)}
    for r, coef, comp in ghostL:
        for ia, (sd, m) in enumerate(allm):
            if sd == "L" and m["comp"] == comp:
                rows.append(r); cols.append(acol[ia]); vals.append(coef / m["z"])
    for r, coef, comp in ghostR:
        for ia, (sd, m) in enumerate(allm):
            if sd == "R" and m["comp"] == comp:
                rows.append(r); cols.append(acol[ia]); vals.append(coef * m["z"])
    rr = 2 * n                                    # consistency: x_0 = sum_L a, x_{n-1} = sum_R a (per component)
    for comp, off in (("u", 0), ("v", n)):
        for sd, jn in (("L", 0), ("R", n - 1)):
            rows.append(rr); cols.append(off + jn); vals.append(1.0)
            for ia, (s2, m) in enumerate(allm):
                if s2 == sd and m["comp"] == comp:
                    rows.append(rr); cols.append(acol[ia]); vals.append(-1.0)
            rr += 1
    A = sp.csr_matrix((np.asarray(vals, dtype=complex), (np.asarray(rows), np.asarray(cols))), shape=(rr, 2 * n + na))
    adm = [ia for ia, (_, m) in enumerate(allm) if m["adm"]]
    inadm = [ia for ia, (_, m) in enumerate(allm) if not m["adm"]]
    keep = list(range(2 * n)) + [acol[ia] for ia in adm]
    Ak = A[:, keep].tocsc()
    if Ak.shape[0] != Ak.shape[1]:
        raise RuntimeError("non-square system %s at omega = %g (admissible %d, inadmissible %d)" % (Ak.shape, omega, len(adm), len(inadm)))
    lu = spla.splu(Ak)
    ins = [ia for ia in inadm if allm[ia][1]["prop"]]
    outs = [ia for ia in adm if allm[ia][1]["prop"]]
    S = np.zeros((len(outs), len(ins)), dtype=complex)
    resid, fluxbal = [], []
    for jj, ia in enumerate(ins):
        m = allm[ia][1]
        amp = 1.0 / math.sqrt(abs(m["vg"]))
        rhs = -A[:, acol[ia]].toarray().ravel() * amp
        x = lu.solve(rhs)
        full = np.zeros(2 * n + na, dtype=complex)
        full[keep] = x
        full[acol[ia]] = amp
        resid.append(float(np.max(np.abs(A @ full))))
        for ii, ib in enumerate(outs):
            S[ii, jj] = full[acol[ib]] * math.sqrt(abs(allm[ib][1]["vg"]))
        fluxbal.append(float(sum(allm[ib][1]["s"] * abs(S[ii, jj]) ** 2 for ii, ib in enumerate(outs)) - m["s"]))

    def lab(ia):
        sd, m = allm[ia]
        return sd + m["comp"]

    info = {"n": n, "mu": mu, "omega": omega, "L": L, "h": h, "VL": VL, "VR": VR,
            "ins": [lab(i) for i in ins], "outs": [lab(i) for i in outs],
            "k_in": [float(allm[i][1]["k"].real) for i in ins], "vg_in": [allm[i][1]["vg"] for i in ins],
            "k_out": [float(allm[i][1]["k"].real) for i in outs], "vg_out": [allm[i][1]["vg"] for i in outs],
            "resid": resid, "fluxbal": fluxbal}
    return S, info


def _entry(S2, info, out_lab, in_lab):
    if out_lab in info["outs"] and in_lab in info["ins"]:
        return float(S2[info["outs"].index(out_lab), info["ins"].index(in_lab)])
    return None


def reduce_omega(S, info, kappa0):
    """The per-omega record: channel sets, |S|^2 by label, r_conf, n_H, thermal forms."""
    S2 = np.abs(S) ** 2
    w = info["omega"]
    x = 2 * math.pi * w / kappa0
    rec = {"omega": w, "ins": info["ins"], "outs": info["outs"], "k_in": info["k_in"], "k_out": info["k_out"],
           "S2": {"%s<-%s" % (o, i): float(S2[a, b]) for a, o in enumerate(info["outs"]) for b, i in enumerate(info["ins"])},
           "fluxbal_max": float(max(abs(f) for f in info["fluxbal"])) if info["fluxbal"] else 0.0,
           "resid_max": float(max(info["resid"])) if info["resid"] else 0.0,
           "boltzmann": math.exp(-x), "planck": 1.0 / math.expm1(x),
           "two_sided": ("Lu" in info["outs"]) and ("Rv" in info["ins"])}
    a = _entry(S2, info, "Ru", "Ru")
    b = _entry(S2, info, "Rv", "Ru")
    rec["r_conf"] = (b / a) if (a is not None and b is not None and a > 0) else None
    nH = _entry(S2, info, "Lu", "Rv")
    rec["n_H"] = nH
    rec["n_H_over_planck"] = (nH / rec["planck"]) if nH is not None else None
    rec["r_conf_over_boltzmann"] = (rec["r_conf"] / rec["boltzmann"]) if rec["r_conf"] is not None else None
    return rec


def scatter_point(N, k_op, kg, omegas, L, h):
    """One grid point over a list of omega; returns the point record (errors per omega are recorded, not raised)."""
    kappa0, c, S = kinematics(N, k_op, kg)
    recs, errs = [], []
    for w in omegas:
        try:
            Sm, info = solve(N, c, kappa0, float(w), L=L, h=h)
            recs.append(reduce_omega(Sm, info, kappa0))
        except Exception as e:                    # recorded; the oracle counts them
            errs.append({"omega": float(w), "error": "%s: %s" % (type(e).__name__, e)})
    two = [r for r in recs if r["two_sided"] and r["n_H"] is not None and r["n_H"] > 0]
    fit = None
    if len(two) >= 3:
        ww = np.array([r["omega"] for r in two])
        yy = np.log(np.array([r["n_H"] for r in two]))
        p = np.polyfit(ww, yy, 1)
        fit = {"slope": float(p[0]), "intercept": float(p[1]),
               "kappa_fit": float(-2 * math.pi / p[0]) if p[0] < 0 else None, "kappa0": kappa0, "n": len(two)}
    rc = [r["r_conf"] for r in recs if r["r_conf"] is not None]
    summ = {"N": N, "k_op": k_op, "kg": kg, "kappa0": kappa0, "c": c, "S": S, "mu": 0.5 * N ** 2, "L": L, "h": h,
            "n_omega": len(omegas), "n_ok": len(recs), "n_err": len(errs),
            "frac_two_sided": float(np.mean([r["two_sided"] for r in recs])) if recs else 0.0,
            "fluxbal_max": max([r["fluxbal_max"] for r in recs] or [0.0]), "resid_max": max([r["resid_max"] for r in recs] or [0.0]),
            "r_conf_min": min(rc) if rc else None, "r_conf_max": max(rc) if rc else None,
            "r_conf_spread_rel": ((max(rc) - min(rc)) / float(np.mean(rc))) if rc else None,
            "n_H_over_planck_min": min([r["n_H_over_planck"] for r in two] or [None], key=lambda v: v if v is not None else 0),
            "n_H_over_planck_max": max([r["n_H_over_planck"] for r in two] or [None], key=lambda v: v if v is not None else 0),
            "thermal_fit_n_H": fit}
    return {"summary": summ, "omega_records": recs, "errors": errs}


def convergence_point(N, k_op, kg, omegas, L, h, conv_h, conv_L):
    """|S|^2 at (L, h), (L, conv_h) and (conv_L, h): max deviation over omega and channels, n_H relative deviation."""
    base = scatter_point(N, k_op, kg, omegas, L, h)
    fine = scatter_point(N, k_op, kg, omegas, L, conv_h)
    wide = scatter_point(N, k_op, kg, omegas, conv_L, h)

    def dev(a, b):
        worst, worst_nH = 0.0, 0.0
        for ra, rb in zip(a["omega_records"], b["omega_records"]):
            if ra["ins"] != rb["ins"] or ra["outs"] != rb["outs"]:
                continue
            for key, v in ra["S2"].items():
                worst = max(worst, abs(v - rb["S2"][key]))
            if ra["n_H"] and rb["n_H"]:
                worst_nH = max(worst_nH, abs(rb["n_H"] / ra["n_H"] - 1))
        same = sum(1 for ra, rb in zip(a["omega_records"], b["omega_records"]) if ra["ins"] == rb["ins"] and ra["outs"] == rb["outs"])
        return {"S2_max_abs_dev": worst, "n_H_max_rel_dev": worst_nH, "omega_same_channels": same, "n_omega": len(a["omega_records"])}

    return {"N": N, "k_op": k_op, "kg": kg, "L": L, "h": h, "conv_h": conv_h, "conv_L": conv_L,
            "h_half": dev(base, fine), "L_wide": dev(base, wide),
            "fluxbal_max": max(base["summary"]["fluxbal_max"], fine["summary"]["fluxbal_max"], wide["summary"]["fluxbal_max"])}


def kaup_null(omegas, hs, L=25.0):
    """Two exact checks of the solver.
    (1) Kaup: N = 1, V = 0 -- the static NLS soliton is reflectionless, |R|^2 = |S[Lu <- Lu]|^2 -> 0 (O(h^4)), |T|^2 -> 1.
    (2) No-soliton null: flow on, psi0 off -- u and v decouple, so every u<->v entry is exactly 0 (no mode conversion
        without the soliton); at (N, k_op) = (1, 1.5), kappa_g = 0.3."""
    kaup = []
    for h in hs:
        for w in omegas:
            S, info = solve(1.0, 1.0, 1.0, float(w), L=L, h=h, Vconst=0.0)
            S2 = np.abs(S) ** 2
            R = _entry(S2, info, "Lu", "Lu")
            T = _entry(S2, info, "Ru", "Lu")
            kaup.append({"h": h, "omega": float(w), "R2": R, "T2": T, "fluxbal_max": max(abs(f) for f in info["fluxbal"])})
    kappa0, c, _ = kinematics(1.0, 1.5, 0.3)
    null = []
    for w in (0.2, 0.8, 1.4):
        S, info = solve(1.0, c, kappa0, w, L=L, h=hs[-1], psi_on=False)
        S2 = np.abs(S) ** 2
        mix = [float(S2[a, b]) for a, o in enumerate(info["outs"]) for b, i in enumerate(info["ins"]) if o[1] != i[1]]
        null.append({"omega": w, "ins": info["ins"], "outs": info["outs"], "uv_mix_max": max(mix) if mix else 0.0,
                     "fluxbal_max": max(abs(f) for f in info["fluxbal"]) if info["fluxbal"] else 0.0})
    return {"kaup": kaup, "null": null}


def findings_check(N, k_op, kg, omegas, ratios, L, h):
    """FINDINGS.md table (N = 2, k_op = 3, delta3 = 0, kappa_g = 0.3): r_conf against the recorded |b/a|^2."""
    P = scatter_point(N, k_op, kg, omegas, L, h)
    rows = []
    for r, ref in zip(P["omega_records"], ratios):
        rows.append({"omega": r["omega"], "r_conf": r["r_conf"], "findings": ref,
                     "rel_dev": (r["r_conf"] / ref - 1) if r["r_conf"] is not None else None,
                     "boltzmann": r["boltzmann"], "ins": r["ins"], "outs": r["outs"], "fluxbal_max": r["fluxbal_max"]})
    return {"N": N, "k_op": k_op, "kg": kg, "L": L, "h": h, "kappa0": P["summary"]["kappa0"], "rows": rows,
            "max_rel_dev": max(abs(x["rel_dev"]) for x in rows if x["rel_dev"] is not None) if rows else None,
            "note": "FINDINGS.md reports unitarity to ~0.4%; this solver's Krein flux balance is at round-off"}
