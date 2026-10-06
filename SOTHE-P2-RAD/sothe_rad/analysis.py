"""Analysis of a SOTHE-P2-RAD results folder: emission rates, the scaled rate function, the recoil model.

Estimators of the relative emission rate R = -d ln Q / dxi of the soliton:
  flux   the radiation behind the soliton (nearest bin of the tail window, distance d): F = I |v_g0(k_b) - V| with I the
         band-filtered intensity, k_b the band wavenumber and V the soliton velocity; it was emitted at the retarded
         time xi_e = xi - d / |v_g0(k_b) - V|, and R(xi_e) = F / Q(xi_e).  It resolves rates far below the oscillation of
         the core norm (a steady emission is a plane wave of uniform intensity behind the soliton).
  core   local slope of -ln Q over xi +- 2.5 (instantaneous; usable where R exceeds ~1e-6).
  box    slope of -ln Qbox over a hold interval (the norm left in the region without absorption; steady state only).
The scaled rate Rs = R / N^2 is a function of the scaled third-order coefficient delta_tilde of the local soliton
(model.py), with N = Q / (2 sqrt(D2)) and D2 = 1 - 6 delta3 k_s, k_s the momentum centroid of the core.

Steps: (1) rates of the prepared solitons over the hold phase; a run is steady when its rate varies by less than 5
percent (relative standard deviation) over the hold phase; (2) the rate function Rs = C dt^p exp(g dt - pi q_r(dt)),
fitted to every hold-phase rate of the prepared solitons (the calibration set: pointwise, equal weight per run, so that
slowly evolving states enter with their instantaneous delta_tilde), with two-parameter forms on the steady runs for
comparison; (3) the launched solitons (N = 1, 1.6, 2), after the launch transient, as the test set: their rates
against the rate function; (4) the recoil model (ode.py) with this rate function and no further parameter, started from
the state of each launched run at xi = 10 and compared with the run (to xi = 40 or 160), and two variants with a
multiplier calibrated on the measured rate at xi = 40; (5) window exponents c of -d ln Q/dxi ~ c/xi against the model;
(6) the loss table of the article reproduced; (7) numerical checks.
"""
import collections
import csv
import json
import math
import os

import numpy as np

from . import __version__, model, ode

PUBLISHED = {(1.0, 0.02): (1.4e-4, -0.020), (1.0, 0.03): (3.1e-4, -0.030), (1.0, 0.04): (5.7e-4, -0.041), (1.0, 0.05): (9.3e-4, -0.052),
             (2.0, 0.025): (9.4e-4, -0.104), (1.0, 0.06): (1.4e-3, -0.064), (1.0, 0.08): (5.5e-3, -0.103), (1.0, 0.10): (2.3e-2, -0.200),
             (2.0, 0.05): (3.5e-2, -0.526)}
STEADY_RELSTD = 0.05        # a prepared run is steady when its flux rate varies by less than this (relative std) over the hold phase
TRANSIENT_END = 20.0        # launched runs: radiation emitted (retarded time) before this belongs to the launch transient
XI0_MODEL = 10.0            # start of the recoil model on the launched runs
WINDOWS_C = [10.0, 20.0, 40.0, 80.0, 160.0]


# ------------------------------------------------------------------------------------------------ helpers
def local_slope(x, y, half=2.5):
    out = np.full(x.size, np.nan)
    for i in range(x.size):
        m = (x >= x[i] - half) & (x <= x[i] + half)
        if np.count_nonzero(m) >= 5:
            out[i] = np.polyfit(x[m], y[m], 1)[0]
    return out


def local_mean(x, y, half=1.0):
    out = np.full(x.size, np.nan)
    for i in range(x.size):
        m = (x >= x[i] - half) & (x <= x[i] + half)
        out[i] = np.mean(y[m])
    return out


def fmt(x, n=3):
    if x is None or not np.isfinite(x):
        return "nan"
    return ("%%.%de" % (n - 1)) % x


def load(folder):
    runs = {}
    with open(os.path.join(folder, "catalogue.json"), encoding="utf-8") as f:
        cat = json.load(f)
    for r in cat:
        d = os.path.join(folder, "runs", r["id"])
        if not os.path.isdir(d):
            continue
        s = dict(np.load(os.path.join(d, "series.npz")))
        with open(os.path.join(d, "meta.json"), encoding="utf-8") as f:
            meta = json.load(f)
        runs[r["id"]] = {"s": s, "meta": meta, "cfg": meta["cfg"], "group": r["group"], "id": r["id"]}
    return runs


def derive(run):
    s, c = run["s"], run["cfg"]
    xi = s["xi"]
    d3s = s["d3"]
    Q = s["Q"]
    ks = local_mean(xi, s["ks_mom"], 1.0)
    D = 1.0 - 6.0 * d3s * ks
    N = Q / (2.0 * np.sqrt(D))
    dt = N * d3s / D ** 1.5
    tc_lab = s["tc"] + c["V"] * xi
    vlab = local_slope(xi, tc_lab, 2.5)
    ph = np.unwrap(s["phase"])
    phrate = local_slope(xi, ph, 2.5)
    Rcore = -local_slope(xi, np.log(Q), 2.5)
    kb = s["kband"]
    In = s["Iband_prof"][:, -1]
    dn = s["Iband_dist"][:, -1]
    with np.errstate(invalid="ignore"):
        vrel = np.abs(model.vg0(kb, d3s) - vlab)
        F = In * vrel
        xie = xi - dn / vrel
    ok = np.isfinite(F) & np.isfinite(xie) & (xie >= xi[0])
    Qe = np.where(ok, np.interp(np.where(ok, xie, xi[0]), xi, Q), np.nan)
    kse = np.where(ok, np.interp(np.where(ok, xie, xi[0]), xi, ks), np.nan)
    d3e = np.where(ok, np.interp(np.where(ok, xie, xi[0]), xi, d3s), np.nan)
    with np.errstate(invalid="ignore", divide="ignore"):
        De = 1.0 - 6.0 * d3e * kse
        Ne = Qe / (2.0 * np.sqrt(De))
        dte = Ne * d3e / De ** 1.5
        Rflux = F / Qe
    kmeas = np.full(xi.size, np.nan)
    for i in range(xi.size):
        if d3s[i] > 0 and np.isfinite(phrate[i]) and np.isfinite(vlab[i]):
            try:
                kmeas[i] = model.k_res_measured(phrate[i], vlab[i], d3s[i])
            except ValueError:
                pass
    return {"xi": xi, "d3": d3s, "Q": Q, "ks": ks, "D2": D, "N": N, "dt": dt, "vlab": vlab, "phrate": phrate, "Rcore": Rcore,
            "kband": kb, "kmeas": kmeas, "vrel": vrel, "F": F, "xie": xie, "Qe": Qe, "kse": kse, "d3e": d3e, "Ne": Ne, "dte": dte,
            "Rflux": Rflux, "Qbox": s["Qbox"], "amp": s["amp"], "E": s["E"], "Mbox": s["Mbox"], "tc_lab": tc_lab}


def rs_fit_function(C, p, g=0.0, variant="drift"):
    return lambda dt: C * dt ** p * math.exp(g * dt - math.pi * model.q_res(dt, variant))


def hold_mask(d, c):
    """Hold phase of a prepared run: xi >= xi_ramp + 40, radiation emitted after xi_ramp + 5 (flux estimator defined)."""
    xr = c["xi_ramp"]
    return (d["xi"] >= xr + 40.0) & np.isfinite(d["Rflux"]) & (d["xie"] >= xr + 5.0)


def fit_rate(dt, rs, variant="drift", nterms=3, w=None):
    """Least squares of ln Rs + pi q_r(dt) = ln C + p ln dt [+ g dt]; returns (C, p, g, residuals of ln Rs)."""
    dt, rs = np.asarray(dt, float), np.asarray(rs, float)
    cols = [np.ones_like(dt), np.log(dt)] + ([dt] if nterms == 3 else [])
    A = np.array(cols).T
    b = np.log(rs) + math.pi * np.array([model.q_res(x, variant) for x in dt])
    sw = np.sqrt(np.asarray(w, float)) if w is not None else np.ones_like(dt)
    coef, *_ = np.linalg.lstsq(A * sw[:, None], b * sw, rcond=None)
    res = b - A @ coef
    return float(math.exp(coef[0])), float(coef[1]), float(coef[2]) if nterms == 3 else 0.0, res


def _stats(v):
    v = np.asarray(v, float)
    return {"n": int(v.size), "median": float(np.median(v)), "p05": float(np.percentile(v, 5)), "p95": float(np.percentile(v, 95)),
            "within_5pc": float(np.mean(np.abs(v) <= math.log(1.05))), "within_10pc": float(np.mean(np.abs(v) <= math.log(1.10)))}


# ------------------------------------------------------------------------------------------------ main
def analyse(folder):
    runs = load(folder)
    if not runs:
        raise SystemExit("no runs in %s" % folder)
    out_dir = os.path.join(folder, "analysis")
    os.makedirs(out_dir, exist_ok=True)
    D = {rid: derive(r) for rid, r in runs.items()}
    summary = {"package": "SOTHE-P2-RAD", "version": __version__, "results": os.path.basename(folder)}

    # (1) rates of the prepared solitons over the hold phase ---------------------------------------------------
    steady = []
    for rid, r in sorted(runs.items(), key=lambda kv: kv[1]["cfg"]["d3"]):
        if r["group"] != "prepared":
            continue
        c, d = r["cfg"], D[rid]
        xr, xe = c["xi_ramp"], c["xi_end"]
        hold = hold_mask(d, c)
        Rf = d["Rflux"][hold]
        mq = (d["xi"] >= xr + 20.0)
        Rc = -np.polyfit(d["xi"][mq], np.log(d["Q"][mq]), 1)[0]
        Rb = -np.polyfit(d["xi"][mq], np.log(d["Qbox"][mq]), 1)[0]
        dtm = float(np.nanmean(d["dte"][hold]))
        Nm = float(np.nanmean(d["Ne"][hold]))
        ksm = float(np.nanmean(d["kse"][hold]))
        Rm = float(np.mean(Rf))
        rb = model.rate_born(Nm, ksm, c["d3"], "drift")
        relstd = float(np.std(Rf) / Rm)
        steady.append({"id": rid, "d3": c["d3"], "N": c["N"], "Nd3": c["N"] * c["d3"], "delta_tilde": dtm,
                       "delta_tilde_min": float(np.nanmin(d["dte"][hold])), "delta_tilde_max": float(np.nanmax(d["dte"][hold])),
                       "N_local": Nm, "ks": ksm, "R_flux": Rm, "R_flux_relstd": relstd,
                       "R_flux_rel_change_over_hold": float((Rf[-1] - Rf[0]) / Rm) if Rf.size > 1 else float("nan"),
                       "R_core_hold": float(Rc), "R_box_hold": float(Rb), "k_band": float(np.nanmean(d["kband"][hold])),
                       "k_measured": float(np.nanmean(d["kmeas"][hold])), "k_model": model.k_res_lab(Nm, ksm, c["d3"], "drift"),
                       "k_model_nodrift": model.k_res_lab(Nm, ksm, c["d3"], "nodrift"),
                       "V": float(np.nanmean(d["vlab"][hold])), "phase_rate": float(np.nanmean(d["phrate"][hold])),
                       "amp": float(np.nanmean(d["amp"][hold])), "R_born": rb, "ratio_to_born": Rm / rb,
                       "loss_ramp_end": float(1 - np.interp(xr, d["xi"], d["Q"]) / d["Q"][0]), "loss_end": float(1 - d["Q"][-1] / d["Q"][0]),
                       "hold": [float(xr + 40.0), float(xe)], "steady": relstd < STEADY_RELSTD})
    _write_csv(os.path.join(out_dir, "steady_rates.csv"), steady)
    summary["steady"] = steady

    # (2) rate function: every hold-phase rate of the prepared solitons ------------------------------------------
    cal = []
    for rid, r in sorted(runs.items(), key=lambda kv: kv[1]["cfg"]["d3"]):
        if r["group"] != "prepared":
            continue
        d, c = D[rid], r["cfg"]
        for i in np.flatnonzero(hold_mask(d, c) & (d["Rflux"] > 0) & np.isfinite(d["dte"])):
            cal.append((rid, float(d["dte"][i]), float(d["Rflux"][i] / d["Ne"][i] ** 2)))
    cnt = collections.Counter(x[0] for x in cal)
    wts = np.array([1.0 / cnt[x[0]] for x in cal])
    dt_c, rs_c = np.array([x[1] for x in cal]), np.array([x[2] for x in cal])
    C, p, g, res = fit_rate(dt_c, rs_c, "drift", 3, wts)
    per_run = {}
    for rid in cnt:
        v = np.array([res[j] for j, x in enumerate(cal) if x[0] == rid])
        per_run[rid] = {"n": int(v.size), "median": float(np.median(v)), "p05": float(np.percentile(v, 5)), "p95": float(np.percentile(v, 95)),
                        "delta_tilde_min": float(min(x[1] for x in cal if x[0] == rid)), "delta_tilde_max": float(max(x[1] for x in cal if x[0] == rid))}
    fitinfo = {"form": "Rs(dt) = C dt^p exp(g dt - pi q_r(dt)), q_r the drift-corrected root",
               "data": "every hold-phase rate (flux estimator) of the prepared solitons, pointwise, equal weight per run",
               "C": C, "p": p, "g": g, "points": int(dt_c.size), "runs": len(cnt),
               "delta_tilde_range": [float(dt_c.min()), float(dt_c.max())],
               "weighted_rms_log_residual": float(math.sqrt(np.sum(wts * res ** 2) / np.sum(wts))),
               "residual": _stats(res), "per_run": per_run,
               "note": "an empirical interpolation formula over its delta_tilde range; C, p and g are strongly correlated and have no separate meaning"}
    sm = [x for x in steady if x["steady"]]
    dts_s = np.array([x["delta_tilde"] for x in sm])
    rss_s = np.array([x["R_flux"] / x["N_local"] ** 2 for x in sm])
    rs_fit = rs_fit_function(C, p, g)
    fitinfo["steady_runs"] = [x["id"] for x in sm]
    fitinfo["steady_residuals"] = {x["id"]: float(math.log(x["R_flux"] / x["N_local"] ** 2 / rs_fit(x["delta_tilde"]))) for x in sm}
    alt = {}
    for nm, var, nt in (("two_parameter_drift", "drift", 2), ("two_parameter_nodrift", "nodrift", 2), ("three_parameter_drift", "drift", 3),
                        ("three_parameter_nodrift", "nodrift", 3)):
        c_, p_, g_, r_ = fit_rate(dts_s, rss_s, var, nt)
        alt[nm] = {"C": c_, "p": p_, "g": g_, "rms_log_residual": float(np.sqrt(np.mean(r_ ** 2))), "max_abs_log_residual": float(np.max(np.abs(r_)))}
    fitinfo["fits_to_the_steady_runs"] = alt
    br = np.array([math.log(rs / model.rate_born_scaled(dt, "drift")) for dt, rs in zip(dts_s, rss_s)])
    cb = np.polyfit(np.log(dts_s), br, 1)
    fitinfo["ratio_to_born_law"] = {"form": "Rs / R_Born = a dt^b (steady runs)", "a": float(math.exp(cb[1])), "b": float(cb[0]),
                                    "rms_log_residual": float(np.sqrt(np.mean((br - np.polyval(cb, np.log(dts_s))) ** 2))),
                                    "ratio_min": float(np.exp(br.min())), "ratio_max": float(np.exp(br.max()))}
    x_inv = np.array([1.0 / x["Nd3"] for x in sm])
    y_ln = np.array([math.log(x["R_flux"]) for x in sm])
    sl = np.polyfit(x_inv, y_ln, 1)
    fitinfo["lnR_vs_inv_Nd3"] = {"slope": float(sl[0]), "intercept": float(sl[1]), "pi_over_2": math.pi / 2,
                                  "note": "ln R against 1/(N delta3) over the steady prepared runs; -pi/2 is the asymptotic slope of the exponent alone"}
    summary["rate_function"] = fitinfo
    with open(os.path.join(out_dir, "rate_function.json"), "w", encoding="utf-8") as f:
        json.dump(fitinfo, f, indent=1)

    # (3) every measured rate against the rate function: calibration (prepared) and test (launched) -------------
    upts = []
    for rid, r in sorted(runs.items()):
        grp = r["group"]
        if grp not in ("prepared", "recoil", "scale", "table"):
            continue
        d, c = D[rid], r["cfg"]
        if grp == "prepared":
            m, role = hold_mask(d, c), "calibration"
        else:
            m, role = (d["xi"] >= 15.0) & np.isfinite(d["Rflux"]) & (d["xie"] >= TRANSIENT_END), "test"
        m = m & (d["Rflux"] > 0) & np.isfinite(d["dte"]) & (d["dte"] >= 0.03)
        for i in np.flatnonzero(m):
            upts.append({"id": rid, "group": grp, "role": role, "estimator": "flux", "N": c["N"], "xi": float(d["xi"][i]), "xi_e": float(d["xie"][i]),
                         "d3": float(d["d3e"][i]), "delta_tilde": float(d["dte"][i]), "R": float(d["Rflux"][i]),
                         "Rs": float(d["Rflux"][i] / d["Ne"][i] ** 2), "Rs_fit": rs_fit(float(d["dte"][i]))})
        if grp == "prepared":         # the slope of the core norm, where it resolves the rate (cross-check)
            mc = (d["xi"] >= c["xi_ramp"] + 20.0) & np.isfinite(d["Rcore"]) & (d["Rcore"] > 1e-6) & (d["dt"] >= 0.03)
            for i in np.flatnonzero(mc):
                upts.append({"id": rid, "group": grp, "role": "cross-check", "estimator": "core", "N": c["N"], "xi": float(d["xi"][i]),
                             "xi_e": float(d["xi"][i]), "d3": float(d["d3"][i]), "delta_tilde": float(d["dt"][i]), "R": float(d["Rcore"][i]),
                             "Rs": float(d["Rcore"][i] / d["N"][i] ** 2), "Rs_fit": rs_fit(float(d["dt"][i]))})
    for u in upts:
        u["log_ratio_to_fit"] = math.log(u["Rs"] / u["Rs_fit"])
    _write_csv(os.path.join(out_dir, "universal_points.csv"), upts)
    lo_c, hi_c = fitinfo["delta_tilde_range"]
    coll = {"by_run": {}, "by_group": {}}
    for rid in sorted({u["id"] for u in upts if u["role"] == "test"}):
        v = [u for u in upts if u["id"] == rid and u["role"] == "test"]
        st_ = _stats([u["log_ratio_to_fit"] for u in v])
        st_.update(N=v[0]["N"], d3=runs[rid]["cfg"]["d3"], delta_tilde_min=min(u["delta_tilde"] for u in v), delta_tilde_max=max(u["delta_tilde"] for u in v),
                   xi_e_min=min(u["xi_e"] for u in v), xi_e_max=max(u["xi_e"] for u in v))
        coll["by_run"][rid] = st_
    for grp in ("table", "recoil", "scale"):
        v = [u["log_ratio_to_fit"] for u in upts if u["role"] == "test" and u["group"] == grp]
        if v:
            coll["by_group"][grp] = _stats(v)
    allt = [u for u in upts if u["role"] == "test"]
    coll["test_all"] = _stats([u["log_ratio_to_fit"] for u in allt])
    inr = [u["log_ratio_to_fit"] for u in allt if lo_c <= u["delta_tilde"] <= hi_c]
    coll["test_within_calibration_range"] = _stats(inr) if inr else None
    coll["calibration_range"] = [lo_c, hi_c]
    coll["cross_check_core_prepared"] = _stats([u["log_ratio_to_fit"] for u in upts if u["role"] == "cross-check"])
    summary["collapse"] = coll

    # (4) recoil model on the launched runs -------------------------------------------------------------
    comp, tracks, cexp = [], {}, []
    rs_born = lambda dt: model.rate_born_scaled(dt, "drift")
    rs_exp = lambda dt: math.exp(-math.pi * model.q_res(dt, "nodrift"))
    for rid, r in sorted(runs.items(), key=lambda kv: (kv[1]["group"], kv[1]["cfg"]["N"] * kv[1]["cfg"]["d3"], kv[1]["cfg"]["N"])):
        if r["group"] not in ("recoil", "table", "scale"):
            continue
        d, c = D[rid], r["cfg"]
        d3 = c["d3"]
        i0 = int(np.argmin(np.abs(d["xi"] - XI0_MODEL)))
        Q0, ks0, Qlaunch = float(d["Q"][i0]), float(d["ks"][i0]), float(d["Q"][0])
        if float(d["dt"][i0]) < lo_c:            # the rate function is not defined below its calibration range
            continue
        xe = c["xi_end"]
        A_ = ode.integrate(d3, Q0, ks0, d["xi"][i0], xe, rs_fit, "drift", "norm")
        tracks[rid] = {"xi": d["xi"], "Q": d["Q"], "ks": d["ks"], "vlab": d["vlab"], "Rcore": d["Rcore"], "A": A_, "Qlaunch": Qlaunch, "Q0": Q0}
        B_ = C_ = None
        if r["group"] == "recoil":
            # variants: the first-order rate law and a pure exponential law (no-drift root, amplitude depleted with the norm),
            # each with a multiplier fixed by the measured (flux) rate emitted at xi = 40 +- 2, started from the state at xi = 40
            iref = int(np.argmin(np.abs(d["xi"] - 40.0)))
            sel = np.isfinite(d["Rflux"]) & (np.abs(d["xie"] - 40.0) <= 2.0)
            Rref = float(np.median(d["Rflux"][sel])) if np.any(sel) else float("nan")
            Qr, ksr = float(d["Q"][iref]), float(d["ks"][iref])
            if math.isfinite(Rref) and Rref > 0:
                sB = ode.calibrate_scale(Rref, Qr, ksr, d3, rs_born, "drift", "norm")
                B_ = ode.integrate(d3, Qr, ksr, d["xi"][iref], xe, rs_born, "drift", "norm", sB)
                Nr = Qr / (2.0 * math.sqrt(1 - 6 * d3 * ksr))
                sC = Rref / (Nr ** 2 * rs_exp(Nr * d3 / (1 - 6 * d3 * ksr) ** 1.5))
                C_ = ode.integrate(d3, Qr, ksr, d["xi"][iref], xe, rs_exp, "nodrift", "deplete", sC, N0=Nr)
                tracks[rid].update(B=B_, C=C_, Qref=Qr, Rref=Rref)
        for xq in (20.0, 40.0, 80.0, 120.0, 160.0):
            if xq > xe + 1e-9:
                continue
            Qg, Qm = float(np.interp(xq, d["xi"], d["Q"])), float(np.interp(xq, A_["xi"], A_["Q"]))
            g_inc, a_inc = 1 - Qg / Q0, 1 - Qm / Q0
            row = {"id": rid, "group": r["group"], "N": c["N"], "d3": d3, "Nd3": c["N"] * d3, "xi": xq,
                   "loss_gnlse": 1 - Qg / Qlaunch, "loss_model": 1 - Qm / Qlaunch,
                   "loss_since_10_gnlse": g_inc, "loss_since_10_model": a_inc, "loss_since_10_rel_err": (a_inc - g_inc) / g_inc if g_inc > 0 else float("nan"),
                   "ks_gnlse": float(np.interp(xq, d["xi"], d["ks"])), "ks_model": float(np.interp(xq, A_["xi"], A_["ks"])),
                   "V_gnlse": float(np.interp(xq, d["xi"], d["vlab"])), "V_model": float(np.interp(xq, A_["xi"], A_["V"]))}
            if B_ is not None and xq >= 40.0:
                for nm, T in (("B", B_), ("C", C_)):
                    row["loss_model" + nm] = float(1 - np.interp(xq, T["xi"], T["Q"]) / Qlaunch)
                    row["ks_model" + nm] = float(np.interp(xq, T["xi"], T["ks"]))
            comp.append(row)
        if r["group"] == "recoil":
            g_c = ode.window_exponents(d["xi"], d["Q"], WINDOWS_C)
            a_c = ode.window_exponents(A_["xi"], A_["Q"], WINDOWS_C)
            for j in range(len(WINDOWS_C) - 1):
                mid = math.sqrt(WINDOWS_C[j] * WINDOWS_C[j + 1])
                cexp.append({"id": rid, "d3": d3, "window": "[%g,%g]" % (WINDOWS_C[j], WINDOWS_C[j + 1]), "c_gnlse": g_c[j], "c_model": a_c[j],
                             "c_inst_model_mid": float(np.interp(mid, A_["xi"], A_["c_inst"]))})
    # momentum balance along the launched runs: dk_s/dxi against -(k_b - k_s) R (core rate, band wavenumber)
    mb = []
    for rid, r in sorted(runs.items(), key=lambda kv: kv[1]["cfg"]["d3"]):
        if r["group"] != "recoil":
            continue
        d = D[rid]
        dks = local_slope(d["xi"], d["ks"], 2.5)
        pred = -(d["kband"] - d["ks"]) * d["Rcore"]
        m = (d["xi"] >= 20.0) & (d["xi"] <= d["xi"][-1] - 5.0) & np.isfinite(dks) & np.isfinite(pred)
        if np.count_nonzero(m) > 10:
            ratio = float(np.median(dks[m] / pred[m]))
            mb.append({"id": rid, "d3": r["cfg"]["d3"], "median_ratio_dks_to_balance": ratio,
                       "p10": float(np.percentile(dks[m] / pred[m], 10)), "p90": float(np.percentile(dks[m] / pred[m], 90)),
                       "ks_change_10_end": float(d["ks"][-1] - np.interp(10.0, d["xi"], d["ks"]))})
    _write_csv(os.path.join(out_dir, "momentum_balance.csv"), mb)
    summary["momentum_balance"] = mb
    _write_csv(os.path.join(out_dir, "recoil_model_vs_gnlse.csv"), comp)
    _write_csv(os.path.join(out_dir, "window_exponents.csv"), cexp)
    summary["recoil"] = comp
    summary["window_exponents"] = cexp
    np.savez_compressed(os.path.join(out_dir, "model_tracks.npz"), **{"%s__%s__%s" % (rid, k, kk): np.asarray(vv)
                                                                         for rid, T in tracks.items() for k, v in T.items()
                                                                         for kk, vv in (v.items() if isinstance(v, dict) else [("v", v)])})

    # (5) the loss table of the article ----------------------------------------------------------------
    tab = []
    for rid, r in runs.items():
        if r["group"] != "table":
            continue
        d, c = D[rid], r["cfg"]
        key = (float(c["N"]), round(float(c["d3"]), 3))
        m = (d["xi"] >= 30.0) & (d["xi"] <= 40.0)
        drift = float(np.polyfit(d["xi"][m], d["tc_lab"][m], 1)[0])
        loss = float(1 - d["Q"][-1] / d["Q"][0])
        pub = PUBLISHED.get(key)
        tab.append({"id": rid, "N": c["N"], "d3": c["d3"], "Nd3": c["N"] * c["d3"], "loss_0_40": loss, "drift_30_40": drift,
                    "printed_loss": pub[0] if pub else None, "printed_drift": pub[1] if pub else None,
                    "loss_matches_print": (float("%.1e" % loss) == pub[0]) if pub else None,
                    "drift_matches_print": (round(drift, 3) == pub[1]) if pub else None})
    tab.sort(key=lambda x: (x["Nd3"], x["N"]))
    _write_csv(os.path.join(out_dir, "loss_table_reproduction.csv"), tab)
    summary["loss_table"] = tab

    # (6) scaling pairs and numerical checks -------------------------------------------------------------
    scal = []
    for (n2, d2_, n1id) in ((2.0, 0.04, "recoil_N1_d0p080"), (2.0, 0.05, "recoil_N1_d0p100"), (2.0, 0.06, "recoil_N1_d0p120")):
        a2 = [rid for rid, r in runs.items() if r["group"] in ("scale", "table") and abs(r["cfg"]["N"] - n2) < 1e-12 and abs(r["cfg"]["d3"] - d2_) < 1e-12]
        if a2 and n1id in D:
            l2 = float(1 - D[a2[0]]["Q"][-1] / D[a2[0]]["Q"][0])
            l1 = float(1 - np.interp(160.0, D[n1id]["xi"], D[n1id]["Q"]) / D[n1id]["Q"][0])
            scal.append({"N2_run": a2[0], "N1_run": n1id, "loss_N2_xi40": l2, "loss_N1_xi160": l1, "rel_diff": (l2 - l1) / l1})
    summary["scaling"] = scal
    checks = []
    base = {x["id"]: x for x in steady}
    for rid, r in runs.items():
        if r["group"] != "check":
            continue
        d, c = D[rid], r["cfg"]
        if c.get("schedule") == "ramp":
            bid = "prep_N1_d%s" % ("%.3f" % c["d3"]).replace(".", "p")
            xr = c["xi_ramp"]
            hold = (d["xi"] >= xr + 40.0) & np.isfinite(d["Rflux"]) & (d["xie"] >= xr + 5.0)
            Rm = float(np.mean(d["Rflux"][hold]))
            ref = base[bid]["R_flux"] if bid in base else float("nan")
            checks.append({"id": rid, "quantity": "steady rate (flux)", "value": Rm, "base": ref, "rel_diff": (Rm - ref) / ref,
                           "k_band": float(np.nanmean(d["kband"][hold]))})
        else:
            bid = "recoil_N1_d%s" % ("%.3f" % c["d3"]).replace(".", "p")
            if bid in D:
                v = float(1 - d["Q"][-1] / d["Q"][0])
                ref = float(1 - D[bid]["Q"][-1] / D[bid]["Q"][0])
                q_diff = float(np.max(np.abs(np.interp(D[bid]["xi"], d["xi"], d["Q"]) - D[bid]["Q"])))
                checks.append({"id": rid, "quantity": "loss over [0, 160]", "value": v, "base": ref, "rel_diff": (v - ref) / ref,
                               "max_abs_Q_diff": q_diff})
    _write_csv(os.path.join(out_dir, "numerical_checks.csv"), checks)
    summary["checks"] = checks
    # conservation diagnostic: Hamiltonian drift of the table runs before radiation reaches the absorber is not meaningful
    # with the absorber; the unit tests check conservation without it.
    with open(os.path.join(out_dir, "summary.json"), "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=1, default=float)
    _summary_md(folder, summary)
    print("analysis written to %s" % out_dir)
    return summary


def _write_csv(path, rows):
    if not rows:
        with open(path, "w", encoding="utf-8") as f:
            f.write("")
        return
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", encoding="utf-8", newline="") as f:
        w = csv.DictWriter(f, fieldnames=keys, lineterminator="\n")
        w.writeheader()
        for r in rows:
            w.writerow({k: (("%.10g" % v) if isinstance(v, float) else v) for k, v in r.items()})


def _summary_md(folder, S):
    L = ["# SOTHE-P2-RAD results: %s" % S["results"], "", "Package SOTHE-P2-RAD %s. All numbers below are computed from the series of this folder by "
         "`python rad.py analyse`; the tables are in `analysis/`. Notation: R is the relative emission rate -d ln Q/dxi of the soliton "
         "(flux estimator unless stated), Rs = R/N^2, dt the scaled third-order coefficient delta_tilde of the local soliton, "
         "log ratio = ln(Rs measured / Rs of the rate function)." % S["version"], ""]
    L += ["## Prepared solitons (N = 1): rates over the hold phase", "",
          "| delta3 | delta_tilde (min-max) | R (flux) | rel. std | steady | R (core) | R (box) | k band | k measured | k model | k model, no drift | R first order | ratio |",
          "|---|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for x in S["steady"]:
        L.append("| %.3f | %.4f (%.4f-%.4f) | %s | %.1e | %s | %s | %s | %.4f | %.4f | %.4f | %.4f | %s | %.1f |" % (
            x["d3"], x["delta_tilde"], x["delta_tilde_min"], x["delta_tilde_max"], fmt(x["R_flux"], 4), x["R_flux_relstd"], "yes" if x["steady"] else "no",
            fmt(x["R_core_hold"]), fmt(x["R_box_hold"]), x["k_band"], x["k_measured"], x["k_model"], x["k_model_nodrift"], fmt(x["R_born"]), x["ratio_to_born"]))
    L += ["", "R (core) resolves rates above about 1e-9 here, R (box) above about 1e-12 (round-off drift of the norm, about 1e-13 per unit xi). "
          "A run is steady when its flux rate varies by less than %g percent (relative standard deviation) over the hold phase." % (100 * STEADY_RELSTD)]
    F = S["rate_function"]
    L += ["", "## Rate function", "",
          "Rs(dt) = C dt^p exp(g dt - pi q_r(dt)), q_r the drift-corrected root: C = %.5g, p = %.4f, g = %.3f, fitted to %d hold-phase rates of %d prepared "
          "runs (pointwise, equal weight per run), dt = %.4f to %.4f; weighted rms log residual %.4f, 5%% to 95%% of the residuals %.4f to %.4f. "
          "An empirical interpolation formula: C, p and g are strongly correlated." % (
              F["C"], F["p"], F["g"], F["points"], F["runs"], F["delta_tilde_range"][0], F["delta_tilde_range"][1], F["weighted_rms_log_residual"],
              F["residual"]["p05"], F["residual"]["p95"]), "",
          "Residual (log ratio) by prepared run:", "", "| run | n | dt range | median | 5% | 95% |", "|---|---|---|---|---|---|"]
    for rid, v in sorted(F["per_run"].items()):
        L.append("| %s | %d | %.4f-%.4f | %+.4f | %+.4f | %+.4f |" % (rid, v["n"], v["delta_tilde_min"], v["delta_tilde_max"], v["median"], v["p05"], v["p95"]))
    L += ["", "Fits to the mean rates of the steady runs only (%s):" % ", ".join(F["steady_runs"]), ""]
    for nm, v in F["fits_to_the_steady_runs"].items():
        L.append("- %s: C = %.5g, p = %.4f, g = %.3f; rms log residual %.4f, largest %.4f" % (nm.replace("_", " "), v["C"], v["p"], v["g"], v["rms_log_residual"], v["max_abs_log_residual"]))
    rb = F["ratio_to_born_law"]
    L += ["", "Ratio of the steady rates to the first-order rate: %.1f to %.1f; as a power law a dt^b: a = %.4g, b = %.3f (rms log residual %.3f). "
          "ln R against 1/(N delta3) over the steady runs: slope %.4f (the exponent alone gives -pi/2 = -%.4f at small dt)." % (
              rb["ratio_min"], rb["ratio_max"], rb["a"], rb["b"], rb["rms_log_residual"], F["lnR_vs_inv_Nd3"]["slope"], math.pi / 2), ""]
    Cl = S["collapse"]
    L += ["## Launched solitons against the rate function (test set)", "",
          "Flux rates emitted after xi = %g (retarded time), i.e. after the launch transient. Calibration range of dt: %.4f to %.4f." % (
              TRANSIENT_END, Cl["calibration_range"][0], Cl["calibration_range"][1]), "",
          "| run | N | delta3 | n | dt range | xi_e range | median | 5% | 95% | within 5% | within 10% |", "|---|---|---|---|---|---|---|---|---|---|---|"]
    for rid, v in Cl["by_run"].items():
        L.append("| %s | %g | %.3f | %d | %.4f-%.4f | %.0f-%.0f | %+.3f | %+.3f | %+.3f | %.2f | %.2f |" % (
            rid, v["N"], v["d3"], v["n"], v["delta_tilde_min"], v["delta_tilde_max"], v["xi_e_min"], v["xi_e_max"], v["median"], v["p05"], v["p95"],
            v["within_5pc"], v["within_10pc"]))
    L += [""]
    for nm, v in [("all test points", Cl["test_all"]), ("test points within the calibration range of dt", Cl["test_within_calibration_range"]),
                  ("core-norm slope of the prepared runs (cross-check)", Cl["cross_check_core_prepared"])] + [("group " + k, v) for k, v in Cl["by_group"].items()]:
        if v:
            L.append("- %s: n = %d, median %+.4f, 5%% %+.4f, 95%% %+.4f, within 5%%: %.2f, within 10%%: %.2f" % (nm, v["n"], v["median"], v["p05"], v["p95"], v["within_5pc"], v["within_10pc"]))
    L += ["", "## Recoil model against the launched runs (started at xi = %g from the state of the run, no free parameter)" % XI0_MODEL, "",
          "Loss since xi = %g: 1 - Q(xi)/Q(%g); total loss: 1 - Q(xi)/Q(0). Runs whose delta_tilde at xi = %g lies below the calibration "
          "range of the rate function are not modelled." % (XI0_MODEL, XI0_MODEL, XI0_MODEL), "",
          "| run | xi | loss since 10, GNLSE | model | rel. err | total loss, GNLSE | model | k_s GNLSE | k_s model | V GNLSE | V model |",
          "|---|---|---|---|---|---|---|---|---|---|---|"]
    for x in S["recoil"]:
        L.append("| %s | %g | %s | %s | %+.3f | %s | %s | %.4f | %.4f | %.4f | %.4f |" % (
            x["id"], x["xi"], fmt(x["loss_since_10_gnlse"], 4), fmt(x["loss_since_10_model"], 4), x["loss_since_10_rel_err"], fmt(x["loss_gnlse"], 4),
            fmt(x["loss_model"], 4), x["ks_gnlse"], x["ks_model"], x["V_gnlse"], x["V_model"]))
    vb = [x for x in S["recoil"] if "loss_modelB" in x and x["xi"] == 160.0]
    if vb:
        L += ["", "Calibrated variants at xi = 160 (multiplier fixed by the flux rate emitted at xi = 40, started at xi = 40): total loss of the run, "
              "of the first-order rate law (B) and of the exponential law with the no-drift root and depleted amplitude (C):", ""]
        for x in vb:
            L.append("- %s: GNLSE %s, B %s, C %s; k_s GNLSE %.4f, B %.4f, C %.4f" % (x["id"], fmt(x["loss_gnlse"], 4), fmt(x["loss_modelB"], 4),
                                                                                   fmt(x["loss_modelC"], 4), x["ks_gnlse"], x["ks_modelB"], x["ks_modelC"]))
    L += ["", "## Window exponents c of -d ln Q/dxi ~ c/xi", "", "| run | window | c GNLSE | c model | 1/K model (mid) |", "|---|---|---|---|---|"]
    for x in S["window_exponents"]:
        L.append("| %s | %s | %s | %s | %s |" % (x["id"], x["window"], fmt(x["c_gnlse"]), fmt(x["c_model"]), fmt(x["c_inst_model_mid"])))
    L += ["", "## Momentum balance along the launched runs", "", "Median ratio of dk_s/dxi to -(k_band - k_s) R_core over xi >= 20 (10% and 90% quantiles):", ""]
    for x in S["momentum_balance"]:
        L.append("- %s: %.3f (%.3f to %.3f); change of k_s from xi = 10 to the end %.4f" % (x["id"], x["median_ratio_dks_to_balance"], x["p10"], x["p90"], x["ks_change_10_end"]))
    L += ["", "## Loss table of the article, reproduced", "", "| N | delta3 | loss [0,40] | printed | drift [30,40] | printed |", "|---|---|---|---|---|---|"]
    for x in S["loss_table"]:
        L.append("| %g | %.3f | %s | %s | %.4f | %s |" % (x["N"], x["d3"], fmt(x["loss_0_40"], 4), x["printed_loss"], x["drift_30_40"], x["printed_drift"]))
    L += ["", "## Scaling pairs (exact symmetry psi = N Psi(N tau, N^2 xi); different discretizations)", ""]
    for x in S["scaling"]:
        L.append("- %s (xi = 40) against %s (xi = 160): loss %s against %s, relative difference %.2e" % (x["N2_run"], x["N1_run"], fmt(x["loss_N2_xi40"], 5),
                                                                                                 fmt(x["loss_N1_xi160"], 5), x["rel_diff"]))
    L += ["", "## Numerical checks", ""]
    for x in S["checks"]:
        L.append("- %s: %s %s (base %s), relative difference %+.2e" % (x["id"], x["quantity"], fmt(x["value"], 5), fmt(x["base"], 5), x["rel_diff"]))
    with open(os.path.join(folder, "SUMMARY.md"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(L) + "\n")
