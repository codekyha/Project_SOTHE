"""Stage S6: port of run_phase2_all.m v0.1.0 (sothe_phase2_matlab), with p2_write_csv.m and p2_write_json.m.

One physics pass (suite.compute_all at the full stability grid, 13 x 13 at n_tau = 360 unless PHASE2_STAB_*
say otherwise), the 29-oracle ledger v0.1.0, and the data products of the MATLAB suite:
  <out>/paper2_data/thermality_ratio.csv          omega, ln_ratio, ratio, alpha2, beta2
                    kappa_drive_independence.csv  drive, kappa_kin
                    kinematic_kappa_flow.csv      kappa_g, kappa_kin, T_H
                    greybody.csv                  omega, Gamma
                    partner_log_negativity.csv    omega, E_N, E_N_cross, nu_minus
                    radiation_rate_table.csv      delta3, maxImOmega
                    bdg_spectrum_d3_000.csv       Re_omega, Im_omega, krein   (|Re w| <= 3 mu)
                    bdg_spectrum_d3_020.csv       Re_omega, Im_omega, krein
                    false_positive_map.csv        kappa_g, drive, beta_abs     (column-major meshgrid flatten)
                    stability_map.csv             N, delta3, maxImOmega        (column-major meshgrid flatten)
                    phase2_summary.json           the MATLAB field order and number format (+ port, status)
  <out>/paper2_figures/fig1_thermality.{pdf,eps,png} ... fig6_partner_entanglement.{pdf,eps,png}, captions.md, captions.tex
                    (suite_figs.py: publication style, no titles or annotations; PNG at 600 dpi)
  <out>/phase2_bundle.npz  (MATLAB: phase2_bundle.mat, -v7); phase2_bundle.mat as well when SciPy is present
  <out>/parity_vs_matlab.json  comparison with reference/matlab_R2025b_U1/ (the MATLAB R2025b products), if present

Caveat carried by the products (G3u audit, S1-1): the delta3 > 0 BdG rates (radiation_rate_table.csv,
bdg_spectrum_d3_020.csv, stability_map.csv, fig2 (b), fig5) linearize about psi0 = N sech(N tau), which is not
stationary once delta3 != 0; oracles P1.4 and P1.5 that pinned them were retired in v0.2.0.  They are
regression products of suite v0.1.0, not results.  RECORDED only (HR-3)."""
import json
import math
import os
import time

import numpy as np

from . import compat as C
from . import oracles, suite

SUITE_PROVENANCE = ("Standalone MATLAB/Octave reconstruction of the SOTHE Phase-2 pipeline; reproduces "
                    "phase2_standalone.py v0.2.0 physics and pins all load-bearing numbers via inline analytic "
                    "oracles. NOT sothe.bdg (private engine, commit 774c7f8). INV1-INV5 enforced: kappa is kinematic "
                    "only; no global-emission-temperature claim.")
CSV_FILES = ("thermality_ratio.csv", "kappa_drive_independence.csv", "kinematic_kappa_flow.csv", "greybody.csv",
             "partner_log_negativity.csv", "radiation_rate_table.csv", "bdg_spectrum_d3_000.csv",
             "bdg_spectrum_d3_020.csv", "false_positive_map.csv", "stability_map.csv")
ARTIFACT_NOTE = ("delta3 > 0 BdG rates about psi0 = N sech(N tau) (radiation_rate_table.csv, bdg_spectrum_d3_020.csv, "
                 "stability_map.csv, fig2 (b), fig5): psi0-linearization artifact (G3u audit, S1-1; oracles P1.4, P1.5 "
                 "retired in v0.2.0); regression products of suite v0.1.0, not results")


# ============================================================ writers (p2_write_csv.m, p2_write_json.m)
def _g17(x):
    """sprintf('%.17g', x) as MATLAB spells it (NaN, Inf, -Inf)."""
    return C.mfmt("%.17g", float(x))


def write_csv(path, header, M):
    """p2_write_csv.m: header row, then every row at %.17g, comma-separated, LF line ends."""
    M = np.asarray(M, dtype=float)
    if M.ndim == 1:
        M = M.reshape(-1, 1)
    with open(path, "w", encoding="ascii", newline="\n") as f:
        f.write(",".join(header) + "\n")
        for row in M:
            f.write(",".join(_g17(v) for v in row) + "\n")


def _num_scalar(x):
    x = float(x)
    if math.isnan(x) or math.isinf(x):
        return "null"
    if x == math.floor(x) and abs(x) < 1e15:
        return "%d" % x
    return "%.17g" % x


def _json_val(v, ind):
    """p2_write_json.m json_val: struct -> object (2-space indent per level), char -> string (only '"' escaped,
    as MATLAB's strrep does), logical -> true/false, numeric scalar -> %d or %.17g, numeric array -> [a, b]
    (column-major flatten), anything else -> null."""
    pad = "  " * ind
    padi = "  " * (ind + 1)
    if isinstance(v, dict):
        if not v:
            return "{}"
        lines = ['%s"%s": %s' % (padi, k, _json_val(x, ind + 1)) for k, x in v.items()]
        return "{\n" + ",\n".join(lines) + "\n" + pad + "}"
    if isinstance(v, str):
        return '"' + v.replace('"', '\\"') + '"'
    if isinstance(v, (bool, np.bool_)):
        return "true" if v else "false"
    if isinstance(v, (int, float, np.integer, np.floating)):
        return _num_scalar(v)
    if isinstance(v, (list, tuple, np.ndarray)):
        a = np.asarray(v)
        if a.dtype == bool:                              # MATLAB: if v (all elements true, not empty)
            return "true" if a.size and bool(np.all(a)) else "false"
        if not np.issubdtype(a.dtype, np.number):
            return "null"
        if a.size == 1:                                  # MATLAB isscalar: a 1 x 1 array prints as a scalar
            return _num_scalar(a.ravel()[0])
        return "[" + ", ".join(_num_scalar(x) for x in a.ravel(order="F")) + "]"
    return "null"


def write_json(path, s):
    """p2_write_json.m: pretty JSON of a struct, newline-terminated."""
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(_json_val(s, 0) + "\n")


def json_text(s):
    return _json_val(s, 0) + "\n"


# ============================================================ the products
def engine_string():
    from .runtime import blas_info
    import platform
    b = blas_info().get("blas") or "BLAS n/a"
    return "Python %s, NumPy %s (%s)" % (platform.python_version(), np.__version__, b.strip())


def write_products(ddir, R):
    """The ten CSV files of run_phase2_all.m, in its order; returns the file names."""
    os.makedirs(ddir, exist_ok=True)
    out = []

    def w(name, header, cols):
        write_csv(os.path.join(ddir, name), header, np.column_stack([np.asarray(c, dtype=float).ravel(order="F") for c in cols]))
        out.append(name)

    w("thermality_ratio.csv", ["omega", "ln_ratio", "ratio", "alpha2", "beta2"],
      [R["omega"], np.log(R["ratio"]), R["ratio"], R["alpha2"], R["beta2"]])
    w("kappa_drive_independence.csv", ["drive", "kappa_kin"], [R["drive"], R["kappa_drive"]])
    w("kinematic_kappa_flow.csv", ["kappa_g", "kappa_kin", "T_H"], [R["kg_flow"], R["kappa_flow"], R["kappa_flow"] / (2 * np.pi)])
    w("greybody.csv", ["omega", "Gamma"], [R["omega"], R["Gamma"]])
    w("partner_log_negativity.csv", ["omega", "E_N", "E_N_cross", "nu_minus"], [R["omega"], R["E_N"], R["E_N_cross"], R["nu_minus"]])
    w("radiation_rate_table.csv", ["delta3", "maxImOmega"], [R["rate_d3"], R["rate_table"]])
    for b in R["bdg"]:
        o = b["out"]
        wv = o["w"][o["win"]]
        kr = o["krein"][o["win"]]
        tag = "%03d" % int(C.mround(b["d3"] * 1000))
        w("bdg_spectrum_d3_%s.csv" % tag, ["Re_omega", "Im_omega", "krein"], [wv.real, wv.imag, kr])
    KG, DR = np.meshgrid(R["fp_kg"], R["fp_drive"])                   # rows = drive, cols = kg (MATLAB meshgrid)
    w("false_positive_map.csv", ["kappa_g", "drive", "beta_abs"], [KG, DR, R["fp_beta"]])
    NN, DD = np.meshgrid(R["stab_N"], R["stab_d3"])                   # rows = delta3, cols = N
    w("stability_map.csv", ["N", "delta3", "maxImOmega"], [NN, DD, R["stab_map"]])
    return out


def summary(P, R, n_pass, n_total, engine):
    """The struct S of run_phase2_all.m (same fields, same order), then two fields of this port."""
    o1 = R["bdg"][0]["out"]
    o2 = R["bdg"][1]["out"]
    S = {"suite": "sothe_phase2_matlab", "version": "0.1.0", "provenance": SUITE_PROVENANCE, "engine": engine,
         "oracles_passed": n_pass, "oracles_total": n_total,
         "kappa_kinematic": R["kappa_kin"], "T_H": R["T_H"], "thermality_R2": R["R2"], "thermality_slope": R["slope"],
         "kappa_fit": R["kappa_fit"], "flux_residual": R["flux_residual"],
         "dS_tot": R["red"]["dS_tot"], "eta_GSL": R["red"]["mean_eta_GSL"],
         "entropy_partition_residual": R["red"]["partition_res"], "norm_drift": R["red"]["norm_drift"],
         "E_N_max_nats": R["E_N_max"], "E_N_at_3kappa": R["E_N_at_3kappa"], "PT_nu_minus_max": float(np.max(R["nu_minus"])),
         "kappa_drive_spread": R["kappa_drive_spread"], "kappa_cont_kg0": np.asarray(R["kappa_cont"]).ravel(),
         "greybody_range": np.array([np.min(R["Gamma"]), np.max(R["Gamma"])]),
         "bdg_gap_over_mu_d3_0": o1["gap"] / o1["mu"], "bdg_maxIm_d3_0": o1["maxIm"], "bdg_maxIm_d3_002": o2["maxIm"],
         "radiation_rate_table": np.asarray(R["rate_table"]).ravel(),
         "config": {"N": P["N"], "delta3": P["d3"], "kappa_g": P["kg"], "k_op": P["k_op"], "n_omega": len(P["omega"]),
                    "n_drive": len(P["drive"]), "bdg_n_tau": P["bdg_n_tau"], "bdg_tau_max": P["bdg_tau_max"],
                    "stab_grid": np.array([len(P["stab_N"]), len(P["stab_d3"])])}}
    from . import PACKAGE, __version__
    S["port"] = "%s %s: Python port of run_phase2_all.m v0.1.0" % (PACKAGE, __version__)
    S["status"] = "RECORDED (HR-3: only MATLAB R2025b output is canonical)"
    return S


def save_bundle(outroot, R, P, name="phase2_bundle", mat=True):
    """<name>.npz: every array and scalar of the physics pass (flat keys); <name>.mat with the structs R and P as
    MATLAB's run_phase2_all saves them (-v7), when MAT is true and SciPy is available.  Returns file names."""
    keys = ("omega", "alpha2", "beta2", "ratio", "Gamma", "E_N", "E_N_cross", "nu_minus", "EN_asymptote_lo", "drive",
            "kappa_drive", "kg_flow", "kappa_flow", "kappa_cont", "fp_kg", "fp_drive", "fp_beta", "rate_d3", "rate_table",
            "stab_N", "stab_d3", "stab_map")
    d = {k: np.asarray(R[k]) for k in keys}
    for k in ("kappa_kin", "T_H", "slope", "intercept", "R2", "kappa_fit", "flux_residual", "E_N_max", "E_N_at_3kappa",
              "kappa_drive_spread", "null_a", "null_b", "offnull", "synfit_slope_err", "synfit_R2"):
        d[k] = np.float64(R[k])
    for k, v in R["red"].items():
        d["red_" + k] = np.float64(v)
    for k, v in R["kin_aux"].items():
        d["kin_aux_" + k] = np.float64(v)
    for b in R["bdg"]:
        tag = "d3_%03d" % int(C.mround(b["d3"] * 1000))
        o = b["out"]
        d["bdg_tau"] = o["tau"]
        for f in ("w", "krein", "win", "ann"):
            d["bdg_%s_%s" % (f, tag)] = o[f]
        for f in ("mu", "maxIm", "gap"):
            d["bdg_%s_%s" % (f, tag)] = np.float64(o[f])
    np.savez_compressed(os.path.join(outroot, name + ".npz"), **d)
    files = [name + ".npz"]
    if not mat:
        return files
    try:
        import scipy.io
    except Exception:
        return files
    try:
        scipy.io.savemat(os.path.join(outroot, name + ".mat"), {"R": _mat_R(R), "P": _mat_P(P)},
                         do_compression=True, oned_as="column")
        files.append(name + ".mat")
    except Exception as e:                              # the .mat copy is a convenience; never fail the stage for it
        print("%s.mat not written: %s" % (name, e))
    return files


def _row(x):
    return np.asarray(x, dtype=float).reshape(1, -1)


def _col(x):
    return np.asarray(x).reshape(-1, 1)


def _mat_R(R):
    M = {"kappa_kin": R["kappa_kin"], "kin_aux": dict(R["kin_aux"]), "T_H": R["T_H"],
         "omega": _col(R["omega"]), "alpha2": _col(R["alpha2"]), "beta2": _col(R["beta2"]), "ratio": _col(R["ratio"]),
         "flux_residual": R["flux_residual"], "slope": R["slope"], "intercept": R["intercept"], "R2": R["R2"],
         "kappa_fit": R["kappa_fit"], "synfit_slope_err": R["synfit_slope_err"], "synfit_R2": R["synfit_R2"],
         "Gamma": _col(R["Gamma"]), "E_N": _col(R["E_N"]), "E_N_cross": _col(R["E_N_cross"]), "nu_minus": _col(R["nu_minus"]),
         "E_N_max": R["E_N_max"], "E_N_at_3kappa": R["E_N_at_3kappa"], "EN_asymptote_lo": _col(R["EN_asymptote_lo"]),
         "drive": _col(R["drive"]), "kappa_drive": _col(R["kappa_drive"]), "kappa_drive_spread": R["kappa_drive_spread"],
         "kg_flow": _col(R["kg_flow"]), "kappa_flow": _col(R["kappa_flow"]), "kappa_cont": _col(R["kappa_cont"]),
         "fp_kg": _row(R["fp_kg"]), "fp_drive": _row(R["fp_drive"]), "fp_beta": np.asarray(R["fp_beta"]),
         "null_a": R["null_a"], "null_b": R["null_b"], "offnull": R["offnull"]}
    bdg = np.empty((1, len(R["bdg"])), dtype=[("d3", "O"), ("out", "O")])
    for i, b in enumerate(R["bdg"]):
        o = b["out"]
        bdg[0, i] = (b["d3"], {"w": _col(o["w"]), "krein": _col(o["krein"]), "mu": o["mu"], "tau": _row(o["tau"]),
                               "maxIm": o["maxIm"], "gap": o["gap"], "win": _col(o["win"]), "ann": _col(o["ann"])})
    M["bdg"] = bdg
    M.update(rate_d3=_col(R["rate_d3"]), rate_table=_col(R["rate_table"]), stab_N=_row(R["stab_N"]),
             stab_d3=_row(R["stab_d3"]), stab_map=np.asarray(R["stab_map"]), red=dict(R["red"]))
    return M


def _mat_P(P):
    M = {}
    for k, v in P.items():
        if k in ("T", "T_retired"):
            continue
        if isinstance(v, dict):
            M[k] = dict(v)
        elif isinstance(v, (list, tuple, np.ndarray)):
            M[k] = _row(v)
        else:
            M[k] = v
    M["T"] = {k: (_row(v) if isinstance(v, (list, tuple, np.ndarray)) else v) for k, v in P["T"].items()}
    return M


# ============================================================ parity with the MATLAB R2025b products
def _read_csv(path):
    with open(path, encoding="ascii") as f:
        head = f.readline().strip().split(",")
    M = np.genfromtxt(path, delimiter=",", skip_header=1, dtype=float)
    if M.ndim == 1:
        M = M.reshape(-1, len(head))
    return head, M


def _match_spectra(a, b):
    """Symmetric nearest-neighbour distance between two eigenvalue sets (the order of eig is not portable)."""
    if a.size == 0 or b.size == 0:
        return float("nan")
    D = np.abs(a.reshape(-1, 1) - b.reshape(1, -1))
    return float(max(np.max(np.min(D, axis=1)), np.max(np.min(D, axis=0))))


def compare_with_reference(ddir, refdir):
    """Per-file comparison of the ten CSV products (and phase2_summary.json) with a reference folder."""
    rep = {"reference": refdir, "files": {}}
    worst = 0.0
    for name in CSV_FILES:
        p, q = os.path.join(ddir, name), os.path.join(refdir, name)
        if not (os.path.exists(p) and os.path.exists(q)):
            rep["files"][name] = {"compared": False}
            continue
        h1, A = _read_csv(p)
        h2, B = _read_csv(q)
        r = {"compared": True, "header_equal": h1 == h2, "rows": [int(A.shape[0]), int(B.shape[0])],
             "byte_identical": open(p, "rb").read() == open(q, "rb").read()}
        if name.startswith("bdg_spectrum"):
            wa = A[:, 0] + 1j * A[:, 1]
            wb = B[:, 0] + 1j * B[:, 1]
            r["eig_set_distance"] = _match_spectra(wa, wb)
            r["krein_counts"] = {"this": [int(np.sum(A[:, 2] > 0)), int(np.sum(A[:, 2] < 0)), int(np.sum(A[:, 2] == 0))],
                                 "reference": [int(np.sum(B[:, 2] > 0)), int(np.sum(B[:, 2] < 0)), int(np.sum(B[:, 2] == 0))]}
            # the annulus 0.1 mu < |Re w| (mu = 1/2 for N = 1): the zero-mode sector is Krein-neutral, its signs are round-off
            aa, ab = np.abs(A[:, 0]) > 0.05, np.abs(B[:, 0]) > 0.05
            r["krein_counts_annulus"] = {"this": [int(np.sum(A[aa, 2] > 0)), int(np.sum(A[aa, 2] < 0))],
                                         "reference": [int(np.sum(B[ab, 2] > 0)), int(np.sum(B[ab, 2] < 0))]}
            r["eig_set_distance_annulus"] = _match_spectra(wa[aa], wb[ab])
            worst = max(worst, r["eig_set_distance"]) if math.isfinite(r["eig_set_distance"]) else worst
        elif A.shape == B.shape:
            cols = {}
            for j, h in enumerate(h2):
                d = np.abs(A[:, j] - B[:, j])
                scale = np.maximum(np.abs(B[:, j]), max(1e-12 * float(np.max(np.abs(B[:, j]))), 1e-300))  # zeros: column scale
                cols[h] = {"identical": bool(np.array_equal(A[:, j], B[:, j])), "max_abs": float(np.max(d)),
                           "max_rel": float(np.max(d / scale))}
                worst = max(worst, cols[h]["max_abs"])
            r["columns"] = cols
        rep["files"][name] = r
    js, jr = os.path.join(ddir, "phase2_summary.json"), os.path.join(refdir, "phase2_summary.json")
    if os.path.exists(js) and os.path.exists(jr):
        a = json.load(open(js))
        b = json.load(open(jr))
        fields = {}
        for k, vb in b.items():
            va = a.get(k)
            if isinstance(vb, (int, float)) and isinstance(va, (int, float)):
                fields[k] = {"this": va, "reference": vb, "abs_diff": abs(va - vb)}
            elif isinstance(vb, list) and isinstance(va, list) and len(va) == len(vb):
                fields[k] = {"max_abs_diff": max(abs(x - y) for x, y in zip(va, vb))}
        rep["summary_json"] = {"reference_engine": b.get("engine"), "fields": fields}
    rep["worst_abs_dev"] = worst
    rep["files_compared"] = sum(1 for v in rep["files"].values() if v.get("compared"))
    rep["byte_identical"] = sorted(k for k, v in rep["files"].items() if v.get("byte_identical"))
    return rep


# ============================================================ driver
def run_phase2_all(outroot, csv_path, P=None, figures=True, refdir=None, verbose=True):
    """run_phase2_all.m v0.1.0: physics -> 29 oracles -> CSVs -> figures -> summary; returns a report dict."""
    ddir = os.path.join(outroot, "paper2_data")
    fdir = os.path.join(outroot, "paper2_figures")
    os.makedirs(ddir, exist_ok=True)
    os.makedirs(fdir, exist_ok=True)
    from . import PACKAGE, __version__
    engine = engine_string()
    print("SOTHE Paper-2 standalone suite v0.1.0, %s %s port  (%s)" % (PACKAGE, __version__, engine))
    print("output root: %s" % outroot)
    t0 = time.time()
    P = suite.phase2_params() if P is None else P
    R = suite.compute_all(P, csv_path)
    t_phys = time.time() - t0
    print("physics pass: %.1f s (stability map %d x %d at n_tau = %d)" % (t_phys, len(P["stab_N"]), len(P["stab_d3"]), P["stab_n_tau"]))
    bundle = save_bundle(outroot, R, P)
    n_pass, n_total, L = oracles.run_oracles(P, R, verbose=verbose)
    files = write_products(ddir, R)
    figs, errs, ffiles = [], {}, []
    if figures:
        from . import suite_figs
        FR = suite_figs.render_all(R, P, fdir)
        figs, errs, ffiles = FR["figures"], FR["errors"], FR["files"]
        print("figures written: %d/6 to %s (PDF, EPS, PNG 600 dpi; captions.md, captions.tex)" % (len(figs), fdir))
        for k, v in errs.items():
            print("WARNING: %s render failed: %s" % (k, v))
    S = summary(P, R, n_pass, n_total, engine)
    write_json(os.path.join(ddir, "phase2_summary.json"), S)
    print("summary written: %s" % os.path.join(ddir, "phase2_summary.json"))
    parity = None
    if refdir and os.path.isdir(refdir):
        parity = compare_with_reference(ddir, refdir)
        with open(os.path.join(outroot, "parity_vs_matlab.json"), "w", encoding="utf-8", newline="\n") as f:
            json.dump(parity, f, indent=1)
        print("parity with %s: %d files compared, byte-identical: %s; worst |dev| %.2e" % (
            refdir, parity["files_compared"], ", ".join(parity["byte_identical"]) or "none", parity["worst_abs_dev"]))
    print("\nDONE  |  oracles %d/%d  |  total %.1f s" % (n_pass, n_total, time.time() - t0))
    if n_pass < n_total:
        print("WARNING: %d oracle(s) failed" % (n_total - n_pass))
    return {"n_pass": n_pass, "n_total": n_total, "ledger": L, "summary": S, "csv": files, "figures": figs,
            "figure_files": ffiles, "figure_errors": errs, "bundle": bundle, "parity": parity, "physics_pass_s": t_phys,
            "stab_grid": [len(P["stab_N"]), len(P["stab_d3"])], "stab_n_tau": P["stab_n_tau"]}
