"""Stage drivers S0..S7 and the in-process runner (port of matlab/pack/stage_S*.m and p2pack_run_stage.m).

Run folder layout:
  <run>/S1_g4a/ S2_g4c/ S3_g12/ S4_oracles/ S5_figures/ S6_phase2/ S7_sweep/
        each with stage_status.json, stage_exit.json, provenance.json, console.log and the stage outputs
  <run>/SUMMARY.md summary.json tokens_proposed.json         written by postprocess.py
A failing stage is recorded (ok = false, exit_code = 1) and the run goes on with the next stage.
No stage writes into the process environment: parameters are passed explicitly."""
import json
import os
import shutil
import sys
import time
import traceback

from . import PACK_ROOT, jsonio, par
from .runtime import Status, Tee

STAGE_DIRS = {"S0": "S0_probe", "S1": "S1_g4a", "S2": "S2_g4c", "S3": "S3_g12", "S4": "S4_oracles", "S5": "S5_figures",
              "S6": "S6_phase2", "S7": "S7_sweep"}
TITLES = {"S0": "runtime probe",
          "S1": "C-05 G4a (port of run_g4a_canonical 0.1.1)",
          "S2": "C-05b T-10 + T-12 (port of run_g4c_loss_scan 0.1.2, sub-grid peak)",
          "S3": "G-12/T-19 provenance of entropy_trajectory.csv (Ref. [25] Python solver)",
          "S4": "T-14 oracle ledger v0.2.0 (and v0.1.0)",
          "S5": "T-13 + T-18 manuscript figures (port of make_paper2_figs 0.2.1)",
          "S6": "phase-2 data products (port of run_phase2_all v0.1.0: 10 CSV, summary JSON, 6 figures)",
          "S7": "Ref. [25] robustness sweep, 300 configurations (port of run_robustness_sweep)"}
DEFAULT_STAGES = ("S1", "S2", "S3", "S4", "S5", "S6", "S7")


def stage_S0(outdir, pack, fast):
    from . import kit, suite
    st = Status("S0", TITLES["S0"], fast)
    try:
        P = st["provenance"]
        print("engine   %s %s (%s)" % (P["engine"], P["version"], P["release"]))
        print("versions %s" % P["versions"])
        print("blas     %s\nlapack   %s\nthreads  %s\nhost     %s" % (P["blas"], P["lapack"], P["threads"], P["host"]))
        t0 = time.time()
        b = suite.bdg_spectrum(1.0, 0.02, 200, 16.0)
        K = kit.kaup_residual(0.7, 8.0, (0.2, 0.1))
        st["checks"] = {"bdg_n200_maxIm": b["maxIm"], "kaup_res": K["res"], "smoke_s": time.time() - t0,
                        "matplotlib": P["versions"].get("matplotlib"), "scipy": P["versions"].get("scipy"),
                        "eligible_canonical": False}
        missing = [m for m in ("matplotlib",) if not P["versions"].get(m)]
        if missing:
            st["notes"].append("missing: %s (stages S5 and S6 need matplotlib)" % ", ".join(missing))
        if not P["versions"].get("scipy"):
            st["notes"].append("SciPy absent: S6 writes phase2_bundle.npz only (no .mat copy); nothing else needs it")
        st["notes"].append("RECORDED only (HR-3): Python port, never canonical")
        st["ok"] = True
        st["pass"] = not missing
        jsonio.write(os.path.join(outdir, "probe.json"), st["checks"])
    except Exception:
        st["error"] = traceback.format_exc()
    st.finish(outdir)
    print("S0 smoke test %s  (%s; RECORDED only)" % ("OK" if st["ok"] and st["pass"] else "FAILED", st["provenance"]["runtime_short"]))
    return st


def stage_S1(outdir, pack, fast):
    from . import g4a
    st = Status("S1", TITLES["S1"], fast)
    try:
        shutil.copyfile(os.path.join(pack, "data", "entropy_trajectory.csv"), os.path.join(outdir, "entropy_trajectory.csv"))
        if fast:
            st["notes"].append("FAST: short T-04 propagation with a looser band")
        R = g4a.run(outdir, fast)
        st["checks"] = {"T01": R["T01_pass"], "T02": R["T02_pass"], "T03": R["T03_pass"], "T04": R["T04_pass"],
                        "T05": R["T05_pass"], "T06": R["T06_pass"], "G4a_all_pass": R["G4a_all_pass"]}
        st["pass"] = bool(R["G4a_all_pass"])
        st["outputs"] = ["g4a_report.json", "T04_series.csv", "bdg_backgrounds.npz"]
        st["ok"] = True
    except Exception:
        st["error"] = traceback.format_exc()
    return st.finish(outdir)


def stage_S2(outdir, pack, fast):
    from . import g4c
    st = Status("S2", TITLES["S2"], fast)
    try:
        if fast:
            st["notes"].append("FAST: no convergence runs")
        R = g4c.run(outdir, fast)
        pts = R["points"]
        op = [p for p in pts if p["N"] == 2 and abs(p["d3"] - 0.05) < 1e-12][0]
        conv = [bool(p["converged"]) for p in pts]
        sc = R["scaling"]["d3_050"]["ratio"]
        st["checks"] = {"all_converged": all(conv), "n_converged": sum(conv), "n_points": len(pts), "scaling_ratio": sc,
                        "scaling_within_5pct": abs(sc - 1) <= 0.05, "recoil_rel_dev_op": op["recoil_rel_dev"],
                        "recoil_within_15pct": op["recoil_rel_dev"] <= 0.15, "recoil_rel_dev_op_grid_peak": op["recoil_rel_dev_grid"],
                        "peak_estimator": R.get("peak_estimator"), "tokens_changed_by_peak": R.get("tokens_changed_by_peak"),
                        "tokens": R["tokens"]}
        c = st["checks"]
        st["pass"] = c["all_converged"] and c["scaling_within_5pct"] and c["recoil_within_15pct"]
        st["outputs"] = ["g4c_loss_report.json", "series/", "newton_solutions.npz"]
        st["ok"] = True
    except Exception:
        st["error"] = traceback.format_exc()
    return st.finish(outdir)


def stage_S3(outdir, pack, fast):
    from . import g12
    st = Status("S3", TITLES["S3"], fast)
    try:
        if fast:
            st["notes"].append("FAST: runs A and C only")
        rep = g12.run(outdir, pack, fast)
        A = rep["runs"][0]
        st["checks"] = {"A_reproduces_archive": A["reproduces_archive"], "A_max_abs_dev_S_tot": A["max_abs_dev_vs_archive"]["S_tot"],
                        "A_DeltaS_tot": A["DeltaS_tot"], "A_eta15": A["eta15"]}
        if rep.get("src_port_check"):
            st["checks"]["A_src_port_max_abs_dev_S_tot"] = rep["src_port_check"]["max_abs_dev_vs_archive"]["S_tot"]
        st["pass"] = bool(A["reproduces_archive"])
        st["outputs"] = ["g12_report.json"] + [r["csv"] for r in rep["runs"]]
        st["ok"] = True
    except Exception:
        st["error"] = traceback.format_exc()
    return st.finish(outdir)


def stage_S4(outdir, pack, fast):
    from . import oracles, phase2, suite
    st = Status("S4", TITLES["S4"], fast)
    try:
        csv = os.path.join(pack, "data", "entropy_trajectory.csv")
        P1 = suite.phase2_params(suite.SMALL_STAB)        # the stability map feeds no oracle: 2 x 2 at n_tau = 48
        t0 = time.time()
        R = suite.compute_all(P1, csv)
        t_phys = time.time() - t0
        print("physics pass (compute_all, suite v0.1.0 port): %.1f s" % t_phys)
        n1, m1, L1 = oracles.run_oracles(P1, R)
        jsonio.write(os.path.join(outdir, "oracle_ledger_v010.json"), {"suite": "0.1.0", "n_pass": n1, "n_total": m1, "rows": L1})
        P2 = suite.phase2_params_v020(suite.SMALL_STAB)
        n2, m2, L2, X2 = oracles.run_oracles_v020(P2, R, csv, fast)
        jsonio.write(os.path.join(outdir, "oracle_ledger_v020.json"), {
            "suite": "0.2.0", "n_pass": n2, "n_total": m2, "rows": L2, "retired_ids": ["P1.4", "P1.5", "P3.3"],
            "retired_pins": P2["T_retired"], "diagnostics": X2, "physics_pass_s": t_phys})
        phase2.save_bundle(outdir, R, P1, name="physics_pass", mat=False)
        st["checks"] = {"v010_pass": n1, "v010_total": m1, "v020_pass": n2, "v020_total": m2,
                        "v020_failed": [r["id"] for r in L2 if not r["pass"]]}
        st["pass"] = n2 == m2
        if fast:
            st["notes"].append("PR-1 on the FAST smoke propagation")
        st["outputs"] = ["oracle_ledger_v010.json", "oracle_ledger_v020.json", "physics_pass.npz"]
        st["ok"] = True
    except Exception:
        st["error"] = traceback.format_exc()
    return st.finish(outdir)


def stage_S5(outdir, pack, fast):
    from . import figures
    st = Status("S5", TITLES["S5"], fast)
    try:
        rep = os.environ.get("P2_G4C_REPORT", "") or os.path.join(os.path.dirname(outdir), "S2_g4c", "g4c_loss_report.json")
        if not os.path.exists(rep):
            st["notes"].append("no g4c_loss_report.json at %s: fig5_radiative_loss not rendered" % rep)
            rep = ""
        M = figures.make_paper2_figs(outdir, pack, rep, st["provenance"])
        names = [f["file"] for f in M["figures"]]
        need = ["fig1_thermality.png", "fig3_gapless_flow.png", "fig2_bdg_krein.png", "fig5_radiative_loss.png", "fig6_partner_entanglement.png"]
        st["checks"] = {"figures_written": names, "missing": [n for n in need if n not in names]}
        st["pass"] = all(n in names for n in need)
        st["notes"].append("RECORDED renders (HR-3); publication style (IOP guidelines), no titles or annotations; "
                           "captions in captions.md / captions.tex; provenance in figs_manifest.json")
        st["outputs"] = (names + [f["vector"] for f in M["figures"]] + [f["eps"] for f in M["figures"]]
                         + M.get("captions", []) + ["figs_manifest.json", "figure_data.npz"])
        st["ok"] = True
    except Exception:
        st["error"] = traceback.format_exc()
    return st.finish(outdir)


def stage_S6(outdir, pack, fast):
    from . import phase2, suite
    st = Status("S6", TITLES["S6"], fast)
    try:
        P = suite.phase2_params()                          # full stability grid unless PHASE2_STAB_* are set
        if fast:
            P = suite.phase2_params({"PHASE2_STAB_NTAU": "96", "PHASE2_STAB_GRIDN": "7", "PHASE2_STAB_GRIDD": "7"})
            st["notes"].append("FAST: stability map 7 x 7 at n_tau = 96 (stability_map.csv and fig5 differ from the full grid)")
        refdir = os.path.join(pack, "reference", "matlab_R2025b_U1")
        rep = phase2.run_phase2_all(outdir, os.path.join(pack, "data", "entropy_trajectory.csv"), P=P,
                                    refdir=refdir if os.path.isdir(refdir) else None)
        with open(os.path.join(outdir, "oracle_ledger_v010_fullgrid.json"), "w", encoding="utf-8", newline="\n") as f:
            json.dump(jsonio.jsonable({"suite": "0.1.0", "n_pass": rep["n_pass"], "n_total": rep["n_total"], "rows": rep["ledger"]}), f, indent=1)
        par_ = rep["parity"] or {}
        st["checks"] = {"oracles_v010": [rep["n_pass"], rep["n_total"]], "csv_written": len(rep["csv"]),
                        "figures_written": len(rep["figures"]), "figure_errors": rep["figure_errors"],
                        "stab_grid": rep["stab_grid"], "stab_n_tau": rep["stab_n_tau"], "bundle": rep["bundle"],
                        "parity_files_compared": par_.get("files_compared"), "parity_byte_identical": par_.get("byte_identical"),
                        "parity_worst_abs_dev": par_.get("worst_abs_dev")}
        st["pass"] = rep["n_pass"] == rep["n_total"] and len(rep["csv"]) == 10 and len(rep["figures"]) == 6
        st["notes"].append(phase2.ARTIFACT_NOTE)
        st["outputs"] = (["paper2_data/" + n for n in rep["csv"]] + ["paper2_data/phase2_summary.json"]
                         + ["paper2_figures/" + n for n in rep["figure_files"]] + rep["bundle"]
                         + (["parity_vs_matlab.json"] if rep["parity"] else []) + ["oracle_ledger_v010_fullgrid.json"])
        st["ok"] = True
    except Exception:
        st["error"] = traceback.format_exc()
    return st.finish(outdir)


def stage_S7(outdir, pack, fast):
    from .ref25 import src
    import numpy as np
    st = Status("S7", TITLES["S7"], fast)
    try:
        wm = (lambda fn, args: par.pmap(fn, args, progress="sweep")) if par.workers() > 1 else None
        out = src.run_robustness_sweep(fast=fast, workers_map=wm)
        chk = src.ref25_checks(fast=fast, workers_map=par.pmap if par.workers() > 1 else None)
        n1, pm, dm = chk["nominal"], chk["phase_matching"], chk["dimensional"]
        print("\nRef. [25] checks (%s):" % ("FAST" if fast else "validate_paper_claims.py settings"))
        print("  nominal N=3.5, d3=0.02: DeltaS_tot %.3f, eta %.3f, dS<0 fraction %.2f, photon drift %.1e  -> %s" % (
            n1["Delta_S_tot"], n1["eta_GSL"], n1["neg_frac"], n1["drift"], "PASS" if n1["passed"] else "FAIL"))
        print("  phase matching d3 = 0.06..0.10: k_RR %s; c_RR %.3f (1.01), R^2 %.3f (0.93), mean k_RR d3 %.3f (1.02), spread %.1f%% "
              "(<5%%)  -> %s" % (" ".join("%.2f" % v for v in pm["k_RR"]), pm["c_RR"], pm["R2"], pm["invariant_mean"],
                                100 * pm["invariant_spread"], "PASS" if pm["passed"] else "FAIL"))
        print("  gnlse_dimensional defaults: DeltaS_tot %.3f, eta %.3f, photon drift %.1e" % (dm["Delta_S_tot"], dm["eta_GSL"], dm["photon_drift"]))
        jsonio.write(os.path.join(outdir, "ref25_checks.json"), chk)
        np.savez_compressed(os.path.join(outdir, "robustness_sweep.npz"), grid_N=out["grid_N"], grid_d3=out["grid_d3"],
                            shapes=np.array(out["shapes"]), mask_widths=out["mask_widths"], mask_orders=out["mask_orders"],
                            DeltaS=out["DeltaS"], eta=out["eta"], photon_drift_max=out["photon_drift_max"])
        with open(os.path.join(outdir, "robustness_table.csv"), "w", encoding="ascii", newline="\n") as f:
            f.write("N_sol,delta3,pulse_shape,mask_width,mask_order,DeltaS_tot,eta_GSL_final,pass,photon_drift_max\n")
            for r in out["rows"]:
                f.write("%.17g,%.17g,%s,%.17g,%d,%.17g,%.17g,%d,%.3e\n" % (r["N_sol"], r["delta3"], r["pulse_shape"], r["mask_width"],
                                                                           r["mask_order"], r["DeltaS_tot"], r["eta_GSL_final"],
                                                                           r["pass"], r["photon_drift_max"]))
        by = {}
        for key, idx in (("pulse_shape", 2), ("mask_scheme", 3)):
            vals = out["shapes"] if key == "pulse_shape" else ["w%.1f_m%d" % (w, m) for w, m in zip(out["mask_widths"], out["mask_orders"])]
            by[key] = {}
            for i, v in enumerate(vals):
                D = np.take(out["DeltaS"], i, axis=idx)
                E = np.take(out["eta"], i, axis=idx)
                by[key][v] = {"DeltaS_range": [float(D.min()), float(D.max())], "eta_range": [float(E.min()), float(E.max())],
                              "pass": int(np.sum((D > 0) & (E > 0.5))), "n": int(D.size)}
        rep = {k: out[k] for k in ("grid_N", "grid_d3", "shapes", "mask_widths", "mask_orders", "DeltaS_range", "eta_range",
                                   "pass_count", "pass_rate", "runtime_s", "settings", "config_order")}
        rep["claims_of_ref25"] = {"DeltaS_range": [0.82, 2.42], "eta_range": [0.52, 0.62], "pass_rate": 1.0,
                                  "source": "Ref. [25] (Class. Quantum Grav. 43, 135014), Sec. 4.2, Eq. (14)"}
        rep["by"] = by
        rep["rows"] = out["rows"]
        rep["note"] = ("run_robustness_sweep.m varies the width and order of the centroid-tracking mask (2.5/12, 3/10, 4/8); "
                       "the per-configuration protocol of Ref. [25] is in its supplementary material")
        jsonio.write(os.path.join(outdir, "robustness_sweep.json"), rep)
        st["checks"] = {"n_configs": len(out["rows"]), "pass_count": out["pass_count"], "pass_rate": out["pass_rate"],
                        "DeltaS_range": out["DeltaS_range"], "eta_range": out["eta_range"], "settings": out["settings"],
                        "nominal_passed": n1["passed"], "phase_matching_passed": pm["passed"], "c_RR": pm["c_RR"],
                        "dimensional_DeltaS_tot": dm["Delta_S_tot"]}
        st["pass"] = len(out["rows"]) == 300
        if fast:
            st["notes"].append("FAST: Nt = 2^12, 1500 steps, xi_max = 6 (the 'fast' flag of run_robustness_sweep.m)")
        st["notes"].append("pass = the 300 configurations ran; the GSL criteria per configuration are reported, not required")
        st["outputs"] = ["robustness_sweep.json", "robustness_sweep.npz", "robustness_table.csv", "ref25_checks.json"]
        st["ok"] = True
    except Exception:
        st["error"] = traceback.format_exc()
    return st.finish(outdir)


DRIVERS = {"S0": stage_S0, "S1": stage_S1, "S2": stage_S2, "S3": stage_S3, "S4": stage_S4, "S5": stage_S5, "S6": stage_S6,
           "S7": stage_S7}


def run_stage(S, rundir, fast=False, pack=PACK_ROOT):
    """Run one stage into <rundir>/<stage dir>, tee its output into console.log, write stage_exit.json."""
    outdir = os.path.join(rundir, STAGE_DIRS[S])
    os.makedirs(outdir, exist_ok=True)
    t0 = time.time()
    out, err = sys.stdout, sys.stderr
    log = open(os.path.join(outdir, "console.log"), "a", encoding="utf-8")
    sys.stdout, sys.stderr = Tee(log, out), Tee(log, err)
    code = 1
    try:
        print("[%s] stage %s -> %s  (workers %d, fast %s)" % (time.strftime("%H:%M:%S", time.gmtime()), S, outdir, par.workers(), fast))
        st = DRIVERS[S](outdir, pack, fast)
        code = 0 if st["ok"] else 1
    except Exception:
        traceback.print_exc()
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        sys.stdout, sys.stderr = out, err
        log.close()
    jsonio.write(os.path.join(outdir, "stage_exit.json"), {"stage": S, "exit_code": code, "wall_s": round(time.time() - t0),
                                                           "runner": os.environ.get("P2_RUNNER", "")})
    return code


def run_all(rundir, stages=DEFAULT_STAGES, fast=False, workers=1, post=True):
    """All stages in this process (sharing one worker pool and the physics-pass cache), then postprocess."""
    par.set_workers(workers)
    codes = {}
    try:
        for S in stages:
            codes[S] = run_stage(S, rundir, fast)
    finally:
        par.shutdown()
    if post:
        from . import postprocess
        try:
            postprocess.main(["--run", rundir])
        except Exception:
            traceback.print_exc()
    return codes
