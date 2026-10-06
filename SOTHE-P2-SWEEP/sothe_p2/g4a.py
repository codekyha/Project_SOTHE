"""Stage S1 (C-05, gate G4a): port of run_g4a_canonical.m 0.1.1 (matlab/kit).

Writes g4a_report.json (same fields, same order as the MATLAB kit script) into the working folder,
which must hold entropy_trajectory.csv.  The five expensive pieces (T-01 spectrum, the two stationary
TOD solitons with their BdG spectra, the psi0 regression spectrum and the T-04 propagation) run in
parallel.  RECORDED only (HR-3)."""
import os
import platform
import sys
import time

import numpy as np

from . import jsonio, kit
from . import compat as C
from .par import pmap
from .suite import bdg_spectrum

KAPPA_EXTRACTED = 1.663588045080       # G3u extracted kappa (T-06 regression)
KAPPA_CLOSED = 1.664100588676          # closed-form kappa, main-text Eq. (7')


def runtime_string():
    return "Python %s, NumPy %s (%s)" % (platform.python_version(), np.__version__, sys.platform)


def _task(name, fast):
    if name == "T01":
        b = bdg_spectrum(1.0, 0.02, 500, 16.0)
        return {"maxIm": b["maxIm"]}
    if name in ("tod_002", "tod_005"):
        d3 = 0.02 if name == "tod_002" else 0.05
        S = kit.tod_soliton(1.0, d3, 500, 16.0)
        B = kit.bdg_spectrum_bg(S["phi"], 1.0, d3, 500, 16.0)
        return {"res_psi0": S["res_psi0"], "newton_res": S["res"][-1], "converged": S["converged"],
                "maxIm": B["maxIm"], "small": float(np.max(np.abs(B["w_small"]))),
                "tau": S["tau"], "phi": S["phi"], "newton_history": np.array(S["res"]), "kmean": S["kmean"],
                "w": B["w"], "krein": B["krein"], "w_small": B["w_small"]}
    if name == "psi0_bg":
        h = 2.0 * 16.0 / 500
        tau = -16.0 + np.arange(500) * h
        P2 = kit.bdg_spectrum_bg(1.0 * C.sech(1.0 * tau), 1.0, 0.02, 500, 16.0)
        return {"maxIm": P2["maxIm"], "w": P2["w"], "krein": P2["krein"], "w_small": P2["w_small"]}
    if name == "T04":
        if fast:
            X = kit.splitstep_loss(2.0, 0.05, 20.0, 0.004, 4096, 400.0, 0.0)
        else:
            X = kit.splitstep_loss(2.0, 0.05, 40.0, 0.002, 8192, 400.0, 0.0)
        return X
    raise ValueError(name)


def run(workdir, fast=False):
    """Run G4a in WORKDIR (must contain entropy_trajectory.csv); returns the report dict."""
    t0 = time.time()
    names = ["T04", "T01", "tod_002", "tod_005", "psi0_bg"]
    out = dict(zip(names, pmap(_task, [(n, fast) for n in names])))
    R = {}
    R["T01_maxIm_frozen"] = out["T01"]["maxIm"]
    R["T01_pass"] = abs(out["T01"]["maxIm"] - 0.09572703807379634) < 1e-9
    s2, s5 = out["tod_002"], out["tod_005"]
    R["T02_res_psi0"] = [s2["res_psi0"], s5["res_psi0"]]
    R["T02_newton_res"] = [s2["newton_res"], s5["newton_res"]]
    R["T02_pass"] = s2["converged"] and s5["converged"] and all(v <= 1e-12 for v in R["T02_newton_res"])
    R["T03_regression_psi0"] = abs(out["psi0_bg"]["maxIm"] - out["T01"]["maxIm"])
    R["T03_true_maxIm"] = [s2["maxIm"], s5["maxIm"]]
    R["T03_true_small"] = [s2["small"], s5["small"]]
    R["T03_pass"] = (R["T03_regression_psi0"] < 1e-10 and all(v <= 1e-10 for v in R["T03_true_maxIm"])
                     and all(v <= 1e-6 for v in R["T03_true_small"]))
    X = out["T04"]
    R["T04_loss_rate"] = X["loss_rate"]
    R["T04_drift"] = X["drift"]
    R["T04_drift_grid"] = X["drift_grid"]          # same fit with the grid-maximum peak (MATLAB 0.1.1 definition)
    R["T04_peak"] = X["peak"]
    R["T04_Q_ratio"] = X["Q_ratio"]
    R["T04_loss_per_hawking_period"] = X["loss_rate"] * 2 * np.pi / KAPPA_CLOSED
    if fast:
        R["T04_pass"] = (X["loss_rate"] > 1e-5) and (X["loss_rate"] <= 1e-3)
    else:
        R["T04_pass"] = (X["loss_rate"] >= 2e-4) and (X["loss_rate"] <= 6e-4)
    G = kit.gsl_eta15(os.path.join(workdir, "entropy_trajectory.csv"))
    R["T05"] = G
    R["T05_pass"] = abs(G["eta15"] - 0.5452) < 5e-4 and abs(G["eta_mean_running"] - 0.6093) < 5e-4
    Gx = kit.gauss_layer(KAPPA_EXTRACTED)
    Gc = kit.gauss_layer(KAPPA_CLOSED)
    R["T06_regression"] = [Gx["EN_max"], Gx["nu_max"], Gx["nbar_peak"], Gx["nbar_band"]]
    R["T06_regression_pass"] = (abs(Gx["EN_max"] - 3.0539) < 5e-5 and abs(Gx["nu_max"] - 0.9071) < 5e-5
                                and abs(Gx["nbar_peak"] - 10.10) < 5e-3 and abs(Gx["nbar_band"] - 0.0512) < 5e-5)
    R["T06_closed_EN_max"] = Gc["EN_max"]
    R["T06_closed_nu_max"] = Gc["nu_max"]
    R["T06_closed_nbar_peak"] = Gc["nbar_peak"]
    R["T06_closed_nbar_band"] = Gc["nbar_band"]
    R["T06_closed_slope"] = Gc["slope"]
    R["T06_closed_EN_at_3kappa"] = Gc["EN_at_3kappa"]
    R["T06_closed_IR_dev"] = Gc["IR_dev"]
    R["T06_tokens"] = Gc["tokens"]
    R["T06_pass"] = R["T06_regression_pass"] and Gc["wronskian_residual"] < 1e-12 and Gc["closed_form_agreement"] < 1e-12
    R["wall_s"] = time.time() - t0
    R["runtime"] = runtime_string()
    R["G4a_all_pass"] = all(R["T0%d_pass" % i] for i in range(1, 7))
    print("\n==== G4a replication (Python port of run_g4a_canonical 0.1.1; RECORDED) ====")
    print("T-01 frozen maxIm           %.14f  pass=%d" % (R["T01_maxIm_frozen"], R["T01_pass"]))
    print("T-02 max|F(psi0)|           %.3e %.3e ; Newton %.1e %.1e  pass=%d" % (tuple(R["T02_res_psi0"]) + tuple(R["T02_newton_res"]) + (R["T02_pass"],)))
    print("T-03 true-bg maxIm          %.2e %.2e ; zero sector %.1e %.1e ; regression %.1e  pass=%d" % (
        tuple(R["T03_true_maxIm"]) + tuple(R["T03_true_small"]) + (R["T03_regression_psi0"], R["T03_pass"])))
    print("T-04 loss %.3e /xi (%.2e per Hawking period), drift %.4f (grid-maximum peak %.4f), Q ratio %.4f  pass=%d" % (
        R["T04_loss_rate"], R["T04_loss_per_hawking_period"], R["T04_drift"], R["T04_drift_grid"], R["T04_Q_ratio"],
        R["T04_pass"]))
    print("T-05 eta15 %.4f ; mean running %.4f ; DeltaS(xi>=2) %.4f ; eta(xi>=2) %.4f  pass=%d" % (
        G["eta15"], G["eta_mean_running"], G["DeltaS_xi_ge_2"], G["eta_xi_ge_2"], R["T05_pass"]))
    print("T-06 regression at kappa=1.663588: E_N^max %.4f nu_max %.4f nbar %.2f / %.4f (pass=%d); at kappa=1.664101: "
          "E_N^max %s nu_max %s nbar %s / %s  pass=%d" % (tuple(R["T06_regression"]) + (R["T06_regression_pass"],) + (
              Gc["tokens"]["G4A_EN_MAX"], Gc["tokens"]["G4A_NU_MAX"], Gc["tokens"]["G4A_NBAR_PEAK"],
              Gc["tokens"]["G4A_NBAR_BAND"], R["T06_pass"])))
    print("G4a all pass: %d   (wall %.0f s, %s)" % (R["G4a_all_pass"], R["wall_s"], R["runtime"]))
    jsonio.write(os.path.join(workdir, "g4a_report.json"), R, pretty=False)
    # raw data for analysis (not part of the MATLAB report)
    write_series_csv(os.path.join(workdir, "T04_series.csv"), X)
    np.savez_compressed(os.path.join(workdir, "bdg_backgrounds.npz"),
                        tau=s2["tau"], phi_d3_0p02=s2["phi"], phi_d3_0p05=s5["phi"],
                        newton_history_d3_0p02=s2["newton_history"], newton_history_d3_0p05=s5["newton_history"],
                        kmean_d3_0p02=s2["kmean"], kmean_d3_0p05=s5["kmean"],
                        w_phi_d3_0p02=s2["w"], krein_phi_d3_0p02=s2["krein"], w_small_phi_d3_0p02=s2["w_small"],
                        w_phi_d3_0p05=s5["w"], krein_phi_d3_0p05=s5["krein"], w_small_phi_d3_0p05=s5["w_small"],
                        w_psi0_d3_0p02=out["psi0_bg"]["w"], krein_psi0_d3_0p02=out["psi0_bg"]["krein"],
                        w_small_psi0_d3_0p02=out["psi0_bg"]["w_small"])
    return R


def write_series_csv(path, X):
    """xi, Q, tau_c, k_s, tau_c_grid of a propagation, one row per record (every 0.5 in xi).
    tau_c is the peak the tokens use (sub-grid); tau_c_grid is the grid maximum (MATLAB 0.1.1)."""
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("xi,Q,tau_c,k_s,tau_c_grid\n")
        for row in zip(X["xi"], X["Q"], X["tau_c"], X["k_s"], X["tau_c_grid"]):
            f.write("%.17g,%.17g,%.17g,%.17g,%.17g\n" % row)
