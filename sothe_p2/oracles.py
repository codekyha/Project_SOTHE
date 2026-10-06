"""Stage S4 (T-14): ports of run_oracles.m 0.1.0 (29 oracles) and run_oracles_v020.m 0.2.0 (35 oracles).

Both ledgers run on one physics pass (suite.compute_all).  Rows keep the MATLAB ids, names, classes,
tolerances and reference values.  RECORDED only (HR-3)."""
import math

import numpy as np

from . import kit
from .par import pmap
from .suite import vg0


def _mono_dec(x):
    return bool(np.all(np.diff(np.asarray(x, dtype=float)) < 0))


def _mono_inc(x):
    return bool(np.all(np.diff(np.asarray(x, dtype=float)) > 0))


def _bnd(x, lo, hi):
    return (x >= lo) and (x <= hi)


def _bdg(R, d3):
    for b in R["bdg"]:
        if (d3 == 0.0 and b["d3"] == 0.0) or (d3 != 0.0 and abs(b["d3"] - d3) < 1e-12):
            return b["out"]
    raise KeyError(d3)


def run_oracles(P, R, verbose=True):
    """run_oracles.m 0.1.0: the 29-oracle acceptance ledger."""
    T = P["T"]
    L = []

    def add(i, name, ok, val, ref, note):
        L.append({"id": i, "name": name, "pass": bool(ok), "val": float(val), "ref": float(ref), "note": note})

    red = R["red"]
    add("P5.1", "dS_tot > 0 & value", abs(red["dS_tot"] - T["dS_tot"]) < P["tol_det"], red["dS_tot"], T["dS_tot"], "frozen trajectory")
    add("P5.2", "<eta_GSL> in [.52,.62] & val", _bnd(red["mean_eta_GSL"], 0.52, 0.62) and abs(red["mean_eta_GSL"] - T["eta_GSL"]) < P["tol_det"],
        red["mean_eta_GSL"], T["eta_GSL"], "coarse-grained GSL")
    add("P5.3", "entropy partition exact", red["partition_res"] < P["tol_exact"], red["partition_res"], 0.0, "S_tot=S_hor+S_rad")
    add("P5.4", "norm drift negligible", red["norm_drift"] < P["tol_det"], red["norm_drift"], 0.0, "|P/P0-1|")
    dc = float(np.max(np.abs(R["kappa_cont"] - np.asarray(P["cont_N"], dtype=float))))
    add("Q1.1", "kappa(kg->0)=N (N=1,2,3)", dc < 3e-3, dc, 0.0, "CQG continuity")
    add("Q1.2", "drive spread ~9.1e-5", abs(R["kappa_drive_spread"] - T["kappa_drive_spread"]) < 0.15 * T["kappa_drive_spread"],
        R["kappa_drive_spread"], T["kappa_drive_spread"], "kappa flat vs drive")
    add("Q1.3", "monotone gapless flow", _mono_dec(R["kappa_flow"]), R["kappa_flow"][-1], R["kappa_flow"][0], "kappa decreasing in kg")
    add("Q1.4", "gapless endpoint < base", R["kappa_flow"][-1] < 0.5 * R["kappa_flow"][0], R["kappa_flow"][-1], R["kappa_flow"][0], "kg=1 vs kg=0")
    add("P2.1", "flux residual ~0", R["flux_residual"] < 1e-12, R["flux_residual"], T["flux_residual"], "|a2-b2-1|")
    add("P2.2", "thermality R^2 = 1", abs(R["R2"] - 1.0) < P["tol_det"], R["R2"], T["thermality_R2"], "ln(ratio) linear")
    add("P2.3", "slope = -2pi/kappa", abs(R["slope"] - (-2 * math.pi / R["kappa_kin"])) < P["tol_eig"], R["slope"], T["thermality_slope"], "detailed balance")
    add("P2.4", "kappa_fit = kappa_kin", abs(R["kappa_fit"] - R["kappa_kin"]) < P["tol_eig"], R["kappa_fit"], R["kappa_kin"], "consistency")
    add("P2.5", "fit exact on synthetic", R["synfit_slope_err"] < 1e-10 and abs(R["synfit_R2"] - 1) < 1e-12, R["synfit_slope_err"], 0.0, "exp(m0 w) recovery")
    G = R["Gamma"]
    add("P3.1", "Gamma in [0,1]", bool(np.all((G >= 0) & (G <= 1))), np.min(G), 0.0, "transmission bound")
    add("P3.2", "Gamma monotone in omega", _mono_inc(G), G[-1], G[0], "barrier opening")
    add("P3.3", "low-omega suppression", abs(np.min(G) - T["greybody_lo"]) < 5e-3 and G[0] < G[-1], np.min(G), T["greybody_lo"], "IR suppression")
    dx = float(np.max(np.abs(R["E_N"] - R["E_N_cross"])))
    add("QI.1", "two closed forms agree", dx < 1e-12, dx, 0.0, "2asinh|b|=2ln(|a|+|b|)")
    add("QI.2", "E_N monotone decreasing", _mono_dec(R["E_N"]), R["E_N"][-1], R["E_N"][0], "decay with omega")
    add("QI.3", "NPT for all omega", bool(np.all(R["nu_minus"] < 1.0)), np.max(R["nu_minus"]), T["PT_nu_minus_max"], "nu_- < 1")
    add("QI.4", "E_N(3kappa) value", abs(R["E_N_at_3kappa"] - T["E_N_at_3kappa"]) < 1e-9, R["E_N_at_3kappa"], T["E_N_at_3kappa"], "tail anchor")
    add("QI.5", "low-omega asymptote", abs(R["E_N"][0] - R["EN_asymptote_lo"][0]) <= 1e-3, R["E_N"][0], R["EN_asymptote_lo"][0], "ln(2k/pi w)")
    add("P4.1", "null point (0,1)=0", R["null_a"] < P["tol_exact"], R["null_a"], 0.0, "no coupling")
    add("P4.2", "null point (.5,0)=0", R["null_b"] < P["tol_exact"], R["null_b"], 0.0, "no drive")
    add("P4.3", "off-null = tanh(.12)>0", abs(R["offnull"] - math.tanh(0.12)) < 1e-12 and R["offnull"] > 0, R["offnull"], math.tanh(0.12), "spurious signal")
    o0, o2 = _bdg(R, 0.0), _bdg(R, 0.02)
    ka = o0["krein"][o0["ann"]]
    add("P1.1", "integrable stable (d3=0)", o0["maxIm"] < 1e-9, o0["maxIm"], T["bdg_maxIm_d3_0"], "max|Im| ~ 0")
    add("P1.2", "gap/mu (d3=0)", abs(o0["gap"] / o0["mu"] - T["bdg_gap_over_mu_d3_0"]) < P["tol_eig"], o0["gap"] / o0["mu"], T["bdg_gap_over_mu_d3_0"], "mass gap")
    add("P1.3", "Krein balance (annulus)", np.sum(ka > 0) == np.sum(ka < 0), np.sum(ka > 0), np.sum(ka < 0), "n_+ = n_- off zero-mode")
    add("P1.4", "rate(d3=0.02)", abs(o2["maxIm"] - T["bdg_maxIm_d3_002"]) < P["tol_eig"], o2["maxIm"], T["bdg_maxIm_d3_002"], "resonant radiation")
    add("P1.5", "rate monotone in d3", _mono_inc(R["rate_table"]), R["rate_table"][-1], R["rate_table"][0], "growth with TOD")
    n_total, n_pass = len(L), sum(r["pass"] for r in L)
    if verbose:
        print("\n================ PHASE-2 ORACLE LEDGER (%d/%d) ================" % (n_pass, n_total))
        print("%-6s %-30s %-6s %-16s %-16s %-11s" % ("ID", "ORACLE", "STATE", "VALUE", "REFERENCE", "|DELTA|"))
        print("-" * 92)
        for r in L:
            print("%-6s %-30s %-6s %-16.8g %-16.8g %-11.3g" % (r["id"], r["name"], "PASS" if r["pass"] else "FAIL", r["val"], r["ref"], abs(r["val"] - r["ref"])))
        print("-" * 92)
        print("RESULT: %d/%d oracles PASS" % (n_pass, n_total))
        if n_pass < n_total:
            print("FAILED: " + " ".join(r["id"] for r in L if not r["pass"]))
        print("==============================================================\n")
    return n_pass, n_total, L


def _tb_task(name, fast):
    from .g4a import _task
    return _task(name, fast)


def run_oracles_v020(P, R, csv_path, fast=False, verbose=True):
    """run_oracles_v020.m 0.2.0: retired P1.4, P1.5, P3.3; modified P1.1, P1.3, P5.x, Q1.1, Q1.2;
    added TB-1..TB-4, PR-1, GS-1, KP-1.  Every row carries a class: ACCURACY, STRUCTURAL or REGRESSION."""
    T = P["T"]
    L = []
    X = {}

    def add(i, name, cls, change, ok, val, ref, note):
        L.append({"id": i, "name": name, "class": cls, "change": change, "pass": bool(ok), "val": float(val),
                  "ref": float(ref), "note": note})

    # the independent expensive pieces first, in parallel: TOD solitons + BdG about them, psi0 regression, PR-1
    names = ["T04", "tod_002", "tod_005", "psi0_bg"]
    tb = dict(zip(names, pmap(_tb_task, [(n, fast) for n in names])))
    G = kit.gsl_eta15(csv_path)
    X["gsl"] = G
    red = R["red"]
    add("P5.1", "DeltaS_tot of the archive", "REGRESSION", "kept", abs(red["dS_tot"] - T["dS_tot"]) < P["tol_det"], red["dS_tot"], T["dS_tot"], "archived entropy_trajectory.csv")
    add("P5.2", "running mean of eta_GSL", "REGRESSION", "modified", abs(red["mean_eta_GSL"] - T["eta_GSL"]) < P["tol_det"], red["mean_eta_GSL"], T["eta_GSL"],
        "regression number only (OD-P1); not the GSL measure")
    add("P5.3", "entropy partition exact", "STRUCTURAL", "kept", red["partition_res"] < P["tol_exact"], red["partition_res"], 0.0, "S_tot = S_hor + S_rad")
    add("P5.4", "norm drift negligible", "ACCURACY", "kept", red["norm_drift"] < P["tol_det"], red["norm_drift"], 0.0, "|P/P0 - 1|")
    add("P5.5", "Eq. (15) window xi >= 2: eta", "REGRESSION", "added", abs(G["eta_xi_ge_2"] - T["eta_xi_ge_2"]) < P["tol_det"], G["eta_xi_ge_2"], T["eta_xi_ge_2"], "window value of Eq. (15)")
    add("P5.6", "Eq. (15) window xi >= 2: DeltaS", "REGRESSION", "added", abs(G["DeltaS_xi_ge_2"] - T["DeltaS_xi_ge_2"]) < P["tol_det"], G["DeltaS_xi_ge_2"], T["DeltaS_xi_ge_2"], "net production after xi = 2")
    add("GS-1", "Eq. (15) value", "ACCURACY", "added", abs(G["eta15"] - T["eta15"]) <= 5e-4, G["eta15"], T["eta15"], "|eta15 - 0.5452| <= 5e-4")
    dtau = P["tau_kin_dtau"]
    v0 = abs(float(vg0(P["k_op"], P["d3"])))
    law = T["nsafe_factor"] * (dtau / v0) ** 2
    Ns = np.asarray(P["cont_N"], dtype=float)
    q11 = (Ns - R["kappa_cont"]) / Ns ** 3
    X["nsafe"] = {"N": Ns, "kappa": R["kappa_cont"], "deficit_over_N3": q11, "law": law, "dtau": dtau, "vg0": v0}
    add("Q1.1", "N-safe law (N-kappa)/N^3", "ACCURACY", "modified", np.max(np.abs(q11 - law)) / law <= 1e-2, np.max(q11), law,
        "(7/12)(dtau/|v_g0|)^2, kg -> 0, N = 1,2,3; rel. tol. 1e-2")
    k0 = R["kin_aux"]["kappa0"]
    lam = k0 / R["kin_aux"]["c_match"]
    env = 0.25 * k0 * (lam * dtau) ** 2
    X["envelope"] = {"kappa0": k0, "lambda": lam, "dtau": dtau, "envelope": env, "spread": R["kappa_drive_spread"]}
    add("Q1.2", "drive spread within envelope", "ACCURACY", "modified", _bnd(R["kappa_drive_spread"], 0, env), R["kappa_drive_spread"], env,
        "(1/4) kappa0 (lambda dtau)^2, lambda = kappa0/c_match")
    add("Q1.3", "monotone gapless flow", "STRUCTURAL", "kept", _mono_dec(R["kappa_flow"]), R["kappa_flow"][-1], R["kappa_flow"][0], "kappa decreasing in kg")
    add("Q1.4", "gapless endpoint < base", "STRUCTURAL", "kept", R["kappa_flow"][-1] < 0.5 * R["kappa_flow"][0], R["kappa_flow"][-1], R["kappa_flow"][0], "kg = 1 vs kg = 0")
    add("P2.1", "flux residual ~0", "ACCURACY", "kept", R["flux_residual"] < 1e-12, R["flux_residual"], T["flux_residual"], "|a2 - b2 - 1|")
    add("P2.2", "thermality R^2 = 1", "ACCURACY", "kept", abs(R["R2"] - 1.0) < P["tol_det"], R["R2"], T["thermality_R2"], "ln(ratio) linear")
    add("P2.3", "slope = -2pi/kappa", "ACCURACY", "kept", abs(R["slope"] - (-2 * math.pi / R["kappa_kin"])) < P["tol_eig"], R["slope"], T["thermality_slope"], "detailed balance")
    add("P2.4", "kappa_fit = kappa_kin", "ACCURACY", "kept", abs(R["kappa_fit"] - R["kappa_kin"]) < P["tol_eig"], R["kappa_fit"], R["kappa_kin"], "consistency")
    add("P2.5", "fit exact on synthetic", "ACCURACY", "kept", R["synfit_slope_err"] < 1e-10 and abs(R["synfit_R2"] - 1) < 1e-12, R["synfit_slope_err"], 0.0, "exp(m0 w) recovery")
    Gm = R["Gamma"]
    add("P3.1", "Gamma in [0,1]", "STRUCTURAL", "kept", bool(np.all((Gm >= 0) & (Gm <= 1))), np.min(Gm), 0.0, "transmission bound")
    add("P3.2", "Gamma monotone in omega", "STRUCTURAL", "kept", _mono_inc(Gm), Gm[-1], Gm[0], "barrier opening")
    dx = float(np.max(np.abs(R["E_N"] - R["E_N_cross"])))
    add("QI.1", "two closed forms agree", "ACCURACY", "kept", dx < 1e-12, dx, 0.0, "2 asinh|b| = 2 ln(|a| + |b|)")
    add("QI.2", "E_N monotone decreasing", "STRUCTURAL", "kept", _mono_dec(R["E_N"]), R["E_N"][-1], R["E_N"][0], "decay with omega")
    add("QI.3", "NPT for all omega", "STRUCTURAL", "kept", bool(np.all(R["nu_minus"] < 1.0)), np.max(R["nu_minus"]), T["PT_nu_minus_max"], "nu_- < 1")
    add("QI.4", "E_N(3kappa) value", "ACCURACY", "kept", abs(R["E_N_at_3kappa"] - T["E_N_at_3kappa"]) < 1e-9, R["E_N_at_3kappa"], T["E_N_at_3kappa"], "tail anchor")
    add("QI.5", "low-omega asymptote", "ACCURACY", "kept", abs(R["E_N"][0] - R["EN_asymptote_lo"][0]) <= 1e-3, R["E_N"][0], R["EN_asymptote_lo"][0], "ln(2k/pi w)")
    add("P4.1", "null point (0,1) = 0", "STRUCTURAL", "kept", R["null_a"] < P["tol_exact"], R["null_a"], 0.0, "no coupling")
    add("P4.2", "null point (.5,0) = 0", "STRUCTURAL", "kept", R["null_b"] < P["tol_exact"], R["null_b"], 0.0, "no drive")
    add("P4.3", "off-null = tanh(.12) > 0", "STRUCTURAL", "kept", abs(R["offnull"] - math.tanh(0.12)) < 1e-12 and R["offnull"] > 0, R["offnull"], math.tanh(0.12), "spurious signal")
    o0 = _bdg(R, 0.0)
    ka = o0["krein"][o0["ann"]]
    add("P1.1", "integrable limit real (d3 = 0)", "STRUCTURAL", "modified", o0["maxIm"] < 1e-9, o0["maxIm"], T["bdg_maxIm_d3_0"], "max|Im| ~ 0; relabelled STRUCTURAL")
    add("P1.2", "gap/mu (d3 = 0)", "ACCURACY", "kept", abs(o0["gap"] / o0["mu"] - T["bdg_gap_over_mu_d3_0"]) < P["tol_eig"], o0["gap"] / o0["mu"], T["bdg_gap_over_mu_d3_0"],
        "the accuracy oracle of the BdG block")
    add("P1.3", "Krein balance (annulus)", "STRUCTURAL", "modified", np.sum(ka > 0) == np.sum(ka < 0), np.sum(ka > 0), np.sum(ka < 0), "n+ = n- off the zero mode; relabelled STRUCTURAL")
    s2, s5 = tb["tod_002"], tb["tod_005"]
    conv = [s2["converged"], s5["converged"]]
    newton = [s2["newton_res"], s5["newton_res"]]
    small = [s2["small"], s5["small"]]
    annIm = [s2["maxIm"], s5["maxIm"]]
    reg = abs(tb["psi0_bg"]["maxIm"] - _bdg(R, 0.02)["maxIm"])
    X["tod"] = {"delta3": [0.02, 0.05], "converged": conv, "newton_res": newton, "zero_sector": small, "annulus_maxIm": annIm, "regression_psi0": reg}
    add("TB-1", "TOD soliton converges", "ACCURACY", "added", all(conv) and all(v <= T["tod_newton_tol"] for v in newton), max(newton), T["tod_newton_tol"],
        "max|F| at (1,0.02,500,16), (1,0.05,500,16)")
    add("TB-2", "zero sector about phi", "ACCURACY", "added", all(v <= T["zero_sector_tol"] for v in small), max(small), T["zero_sector_tol"], "max|w_small|, both delta3")
    add("TB-3", "annulus reality about phi", "ACCURACY", "added", all(v <= T["annulus_tol"] for v in annIm), max(annIm), T["annulus_tol"], "maxIm on 0.1 mu < |Re w| <= 3 mu, both delta3")
    add("TB-4", "bg routine regression at psi0", "STRUCTURAL", "added", reg <= T["bg_regression_tol"], reg, 0.0, "|maxIm_bg(psi0) - maxIm_suite| at delta3 = 0.02")
    Y = tb["T04"]
    if fast:
        band = [1e-5, 1e-3]
        pnote = "FAST smoke run (xi_end 20, dxi 0.004, n 4096); band (1e-5, 1e-3]; never canonical"
    else:
        band = T["loss_band"]
        pnote = "window fit over xi in [10, 40]; T-04 protocol"
    X["loss"] = {"loss_rate": Y["loss_rate"], "drift": Y["drift"], "Q_ratio": Y["Q_ratio"], "band": band, "fast": bool(fast)}
    add("PR-1", "operating-point radiative loss", "ACCURACY", "added", _bnd(Y["loss_rate"], band[0], band[1]), Y["loss_rate"], float(np.mean(band)), pnote)
    K = kit.kaup_residual(0.7, 8.0, (0.2, 0.1, 0.05, 0.025))
    X["kaup"] = K
    add("KP-1", "Kaup continuum residual order", "ACCURACY", "added", np.min(K["order"]) >= T["kaup_min_order"], np.min(K["order"]), T["kaup_min_order"],
        "res %.3e %.3e %.3e %.3e; recorded 4.121e-3 2.771e-4 1.765e-5 1.108e-6" % tuple(K["res"]))
    n_total, n_pass = len(L), sum(r["pass"] for r in L)
    if verbose:
        print("\n================ PHASE-2 ORACLE LEDGER v0.2.0 (%d/%d) ================" % (n_pass, n_total))
        print("%-6s %-32s %-10s %-8s %-6s %-16s %-16s" % ("ID", "ORACLE", "CLASS", "CHANGE", "STATE", "VALUE", "REFERENCE"))
        print("-" * 100)
        for r in L:
            print("%-6s %-32s %-10s %-8s %-6s %-16.8g %-16.8g" % (r["id"], r["name"], r["class"], r["change"], "PASS" if r["pass"] else "FAIL", r["val"], r["ref"]))
        print("-" * 100)
        msg = "RESULT: %d/%d oracles PASS" % (n_pass, n_total)
        if n_pass < n_total:
            msg += "   FAILED: " + " ".join(r["id"] for r in L if not r["pass"])
        print(msg)
        print("retired: P1.4, P1.5, P3.3 (pins in P.T_retired)")
    return n_pass, n_total, L, X
