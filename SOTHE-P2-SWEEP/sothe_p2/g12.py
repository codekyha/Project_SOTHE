"""Stage S3 (G-12 / T-19): where entropy_trajectory.csv comes from.

Port of stage_S3_g12.m (P2_altay_pack 1.0.0).  The MATLAB stage ran Ref. [25]'s MATLAB solver; this one runs
Ref. [25]'s own Python solver (SOTHE_pkg v1.0.0 validation/gnlse_dimensionless.py, unchanged, sothe_p2/ref25/),
the arithmetic closest to the generator of the archive, at
  A  the settings recorded by sothe-phase0/results/validation_run.json: N_sol = 3.5, delta3 = 0.02,
     Nt = 2^13, 4000 steps, xi_max = 12, 200 saves, window 20, sech launch, mask width 3, order 10;
  B  the same physics on Nt = 2^14, 6000 steps (the supplement's convergence claim);
  C  the settings the manuscript declares (N = 2, delta3 = 0.05), stated grid;
  D  C on the refined grid;
writes each run as a CSV in the archive's column format, evaluates Eq. (15) on it and compares it
sample by sample with the archive.  Pass: run A reproduces the archive, max |dS| < 1e-4.
RECORDED references: NumPy 2.2.6 (PI machine) A 1.8142850 / 0.5452 (max|dS_tot| 4.5e-8), B 1.9302016 / 0.5410,
C 1.1240017 / 0.5542, D 1.1764210 / 0.5550."""
import os
import time

import numpy as np

from . import compat as C          # (also installs np.trapezoid on NumPy 1.x before the validation port is imported)
from . import jsonio, kit
from .par import pmap
from .suite import read_csv_numeric

RUNS = [("A_archived_settings", 3.5, 0.02, 8192, 4000),
        ("B_archived_refined", 3.5, 0.02, 16384, 6000),
        ("C_declared_settings", 2.0, 0.05, 8192, 4000),
        ("D_declared_refined", 2.0, 0.05, 16384, 6000)]
SOLVER = ("SOTHE_pkg v1.0.0 validation/gnlse_dimensionless.py (Python port of Ref. [25], unchanged; "
          "sothe_p2/ref25/gnlse_dimensionless.py)")
SRC_SOLVER = "sothe_p2/ref25/src.py gnlse_dimensionless (statement-by-statement port of Ref. [25] src/gnlse_dimensionless.m)"
MASK_RULE = ("M_S(k) = exp{-[(k - k_c)/w]^m}, w = 3, m = 10, k_c = int k P dk / int P dk with P = |psi~(k)|^2; "
             "M_R = 1 - M_S (a smooth window, not a disjoint partition); p_X = M_X P / int M_X P dk; "
             "S_X = -int p_X ln p_X dk over p_X > 1e-20 (trapezoidal rule on the fftshifted k grid) "
             "[SOTHE_pkg v1.0.0: soliton_mask, spectral_entropy]")
DECISION = ("PI (G-12/T-19): either restate the trajectory parameters in Sec. III F and supplement Sec. 2 (runs A and B), "
            "or regenerate the trajectory at the declared parameters (runs C and D) and update the frozen GSL values; "
            "the supplement's convergence sentence must follow the run B/D numbers.")


def _dispatch(fn, args):
    return fn(*args)


def _run(tag, Nsol, d3, Nt, ns):
    from .ref25 import gnlse_dimensionless as g
    t0 = time.time()
    o = g.gnlse_dimensionless(N_sol=Nsol, delta3=d3, xi_max=12.0, n_steps=ns, Nt=Nt, tau_window=20.0, n_save=200,
                              pulse_shape="sech", mask_width=3.0, mask_order=10, verbose=False)
    Pn = np.trapezoid(o.spec_hist, o.k, axis=0)
    eta = np.array(o.eta_gsl, dtype=float)
    eta[0] = np.nan                                     # archive convention
    M = np.column_stack([o.xi, o.S_hor, o.S_rad, o.S_tot, eta, Pn])
    return {"M": M, "DeltaS_tot": float(o.Delta_S_tot), "photon_number": np.asarray(o.photon_number, dtype=float),
            "wall_s": time.time() - t0}


def _run_src(Nsol, d3, Nt, ns):
    """Run A with the port of the MATLAB solver (ref25/src.py): the same physics through MATLAB's statements."""
    from .ref25 import src
    t0 = time.time()
    o = src.gnlse_dimensionless(N_sol=Nsol, delta3=d3, xi_max=12.0, n_steps=ns, Nt=Nt, tau_window=20.0, n_save=200,
                                pulse_shape="sech", mask_width=3.0, mask_order=10, verbose=False)
    Pn = C.mtrapz(o["k"], o["spec_hist"])
    eta = np.array(o["eta_gsl"], dtype=float)
    eta[0] = np.nan
    M = np.column_stack([o["xi"], o["S_hor"], o["S_rad"], o["S_tot"], eta, Pn])
    return {"M": M, "DeltaS_tot": float(o["Delta_S_tot"]), "wall_s": time.time() - t0}


def _write_csv(path, M):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("xi,S_hor,S_rad,S_tot,eta_GSL,P_norm\n")
        for row in M:
            f.write(",".join("nan" if np.isnan(v) else "%.18e" % v for v in row) + "\n")


def run(outdir, pack_root, fast=False):
    ref = read_csv_numeric(os.path.join(pack_root, "data", "entropy_trajectory.csv"))
    runs = [RUNS[0], RUNS[2]] if fast else RUNS
    jobs = [(_run, r) for r in runs] + ([] if fast else [(_run_src, RUNS[0][1:])])
    allres = pmap(_dispatch, jobs)
    res = allres[:len(runs)]
    rep = {"reference_csv": "data/entropy_trajectory.csv", "reference_DeltaS_tot": float(ref[-1, 3] - ref[0, 3]),
           "solver": SOLVER, "runs": []}
    for (tag, Nsol, d3, Nt, ns), o in zip(runs, res):
        M = o["M"]
        name = "traj_%s_N%.1f_d%.2f_Nt%d_ns%d.csv" % (tag, Nsol, d3, Nt, ns)
        csv = os.path.join(outdir, name)
        _write_csv(csv, M)
        G = kit.gsl_eta15(csv)
        r = {"tag": tag, "N_sol": Nsol, "delta3": d3, "Nt": Nt, "n_steps": ns, "xi_max": 12.0, "tau_window": 20.0,
             "n_save": 200, "mask_width": 3.0, "mask_order": 10, "csv": name, "wall_s": o["wall_s"],
             "DeltaS_tot": o["DeltaS_tot"], "eta15": G["eta15"], "eta_mean_running": G["eta_mean_running"],
             "DeltaS_xi_ge_2": G["DeltaS_xi_ge_2"], "eta_xi_ge_2": G["eta_xi_ge_2"], "rise_first4": G["rise_first4"],
             "photon_number_0": float(o["photon_number"][0]),
             "photon_drift_max": float(np.max(np.abs(o["photon_number"] / o["photon_number"][0] - 1))),
             "P_norm_drift_max": float(np.max(np.abs(M[:, 5] / M[0, 5] - 1)))}
        n = min(M.shape[0], ref.shape[0])
        dev = {"xi": float(np.max(np.abs(M[:n, 0] - ref[:n, 0]))),
               "S_hor": float(np.max(np.abs(M[:n, 1] - ref[:n, 1]))),
               "S_rad": float(np.max(np.abs(M[:n, 2] - ref[:n, 2]))),
               "S_tot": float(np.max(np.abs(M[:n, 3] - ref[:n, 3]))),
               "eta_GSL": float(np.max(np.abs(M[1:n, 4] - ref[1:n, 4]))),
               "P_norm_rel": float(np.max(np.abs(M[:n, 5] / ref[:n, 5] - 1)))}
        r["max_abs_dev_vs_archive"] = dev
        mx = max(dev["S_hor"], dev["S_rad"], dev["S_tot"])
        r["reproduces_archive"] = mx < 1e-4
        if mx < 1e-7:
            r["agreement"] = "same to 1e-7 (same FFT arithmetic as the generator)"
        elif mx < 1e-4:
            r["agreement"] = "reproduced to FFT-rounding level (1e-7 .. 1e-4)"
        else:
            r["agreement"] = "NOT the archived trajectory"
        rep["runs"].append(r)
        print("%-20s N=%.1f d3=%.2f Nt=%5d ns=%4d  DeltaS_tot=%.7f  eta15=%.4f  max|dS_tot| vs archive=%.2e  (%.0f s)" % (
            tag, Nsol, d3, Nt, ns, r["DeltaS_tot"], r["eta15"], dev["S_tot"], r["wall_s"]))
    if not fast:
        o = allres[-1]
        M = o["M"]
        n = min(M.shape[0], ref.shape[0])
        rep["src_port_check"] = {
            "solver": SRC_SOLVER, "run": "A_archived_settings", "DeltaS_tot": o["DeltaS_tot"], "wall_s": o["wall_s"],
            "max_abs_dev_vs_archive": {"S_hor": float(np.max(np.abs(M[:n, 1] - ref[:n, 1]))),
                                       "S_rad": float(np.max(np.abs(M[:n, 2] - ref[:n, 2]))),
                                       "S_tot": float(np.max(np.abs(M[:n, 3] - ref[:n, 3])))},
            "max_abs_dev_vs_validation_port": float(np.max(np.abs(M[:, 3] - res[0]["M"][:, 3]))),
            "note": "not a pass criterion: MATLAB's trapz order and linspace in place of the validation port's"}
        print("%-20s src.py port of the MATLAB solver: DeltaS_tot=%.7f  max|dS_tot| vs archive=%.2e, vs the validation "
              "port %.2e  (%.0f s)" % ("A (MATLAB-src port)", o["DeltaS_tot"], rep["src_port_check"]["max_abs_dev_vs_archive"]["S_tot"],
                                        rep["src_port_check"]["max_abs_dev_vs_validation_port"], o["wall_s"]))
    A = rep["runs"][0]
    rep["archive_reproduced_by_A"] = A["reproduces_archive"]
    if A["reproduces_archive"]:
        rep["verdict"] = ("entropy_trajectory.csv is reproduced by Ref. [25]'s solver at N_sol = 3.5, delta3 = 0.02, Nt = 2^13, "
                          "4000 steps, xi_max = 12, window 20, mask width 3, order 10 (max|dS_tot| = %.1e); the manuscript "
                          "(Sec. III F, supplement Sec. 2) states N = 2, delta3 = 0.05." % A["max_abs_dev_vs_archive"]["S_tot"])
    else:
        rep["verdict"] = "run A does not reproduce the archive on this runtime; compare the RECORDED values in the header"
    rep["mask_rule"] = MASK_RULE
    rep["decision_needed"] = DECISION
    jsonio.write(os.path.join(outdir, "g12_report.json"), rep)
    return rep
