"""Run summary: SUMMARY.md, summary.json and tokens_proposed.json for one run folder.

Usage: python -m sothe_p2.postprocess --run RUNDIR

Reads the stage outputs of RUNDIR, states which criteria held, lists the key numbers and the token values
each report supplies, and indexes the raw-data files kept for analysis.  Nothing here edits a manuscript.
Every value is RECORDED (HR-3).  Standard library only."""
import json
import math
import os
import sys

sys.dont_write_bytecode = True

STAGES = [("S1", "S1_g4a", "C-05 gate G4a"), ("S2", "S2_g4c", "C-05b G4c tokens"),
          ("S3", "S3_g12", "G-12 trajectory provenance"), ("S4", "S4_oracles", "T-14 oracle ledger v0.2.0"),
          ("S5", "S5_figures", "T-13 + T-18 manuscript figures"), ("S6", "S6_phase2", "phase-2 data products (suite v0.1.0)"),
          ("S7", "S7_sweep", "Ref. [25] robustness sweep")]


def load(path, default=None):
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError):
        return default


def as_list(x):
    if x is None:
        return []
    return x if isinstance(x, list) else [x]


def mx(x):
    v = [float(a) for a in as_list(x) if a is not None]
    return max(v) if v else None


def sci2(x, up=False):
    """Two significant figures as LaTeX; up=True rounds up (bound tokens)."""
    if x is None:
        return None
    if x == 0:
        return "0"
    e = int(math.floor(math.log10(abs(x))))
    q = 10.0 ** (e - 1)
    m = (math.ceil(abs(x) / q - 1e-9) if up else round(abs(x) / q)) * q
    if m >= 10.0 ** (e + 1):
        e += 1
    mant = m / 10.0 ** e
    return "%.1f\\times10^{%d}" % (math.copysign(mant, x), e)


def fmt(x):
    if x is None:
        return "n/a"
    if isinstance(x, bool):
        return str(x)
    if isinstance(x, (int, float)):
        return "%.6g" % x
    if isinstance(x, list):
        return "[" + ", ".join(fmt(a) for a in x) + "]"
    return str(x)


def main(argv=None):
    a = list(sys.argv[1:] if argv is None else argv)
    run = a[a.index("--run") + 1]
    S = {"run": os.path.basename(os.path.normpath(run)), "stages": {}, "decisions": []}
    md = []
    prov = None
    for sid, d, title in STAGES:
        st = load(os.path.join(run, d, "stage_status.json"))
        ex = load(os.path.join(run, d, "stage_exit.json"), {})
        if st is None:
            S["stages"][sid] = {"title": title, "state": "not run" if not os.path.isdir(os.path.join(run, d)) else "no status (crashed or timed out)",
                                "exit_code": ex.get("exit_code")}
            continue
        prov = prov or st.get("provenance")
        err = (st.get("error") or "").strip().splitlines()
        S["stages"][sid] = {"title": title, "ran": st.get("ok"), "pass": st.get("pass"), "exit_code": ex.get("exit_code"),
                            "wall_s": st.get("wall_s"), "error": err[-1] if err else "", "notes": st.get("notes")}
    S["provenance"] = prov or {}
    S["canonical_eligible"] = False
    md.append("# SOTHE-P2 run %s" % S["run"])
    md.append("")
    if prov:
        th = prov.get("threads") or {}
        md.append("Runtime: %s %s (%s; Matplotlib %s) on %s; BLAS %s; workers %s; job %s; %s %s; runner %s." % (
            prov.get("engine"), prov.get("version"), prov.get("release"), (prov.get("versions") or {}).get("matplotlib"),
            prov.get("host"), prov.get("blas") or "n/a", th.get("workers"), prov.get("slurm_job_id") or "local",
            prov.get("package") or "SOTHE-P2", prov.get("pack_version"), prov.get("runner") or "n/a"))
        md.append("")
    md.append("**Status: RECORDED (HR-3).** Every number and figure below is RECORDED evidence of a Python runtime; "
              "only MATLAB R2025b output is canonical.")
    md.append("")
    md.append("| Stage | Item | Ran | Criteria met | Exit | Wall (s) |")
    md.append("|---|---|---|---|---|---|")
    for sid, d, title in STAGES:
        x = S["stages"][sid]
        if "ran" not in x:
            md.append("| %s | %s | %s | | %s | |" % (sid, title, x["state"], fmt(x.get("exit_code"))))
        else:
            md.append("| %s | %s | %s | %s | %s | %.0f |" % (sid, title, x["ran"], x["pass"], fmt(x["exit_code"]), x["wall_s"] or 0))
    md.append("")
    fast = bool(prov and prov.get("fast_mode"))
    if fast:
        md.append("Fast run: short propagations and no convergence runs, so S2 cannot meet its criteria and T-04 / PR-1 use the smoke band; "
                  "S6 uses a 7 x 7 stability map and S7 the fast sweep grid.")
        md.append("")
    tokens = {}

    # ---- S1 ------------------------------------------------------------------------------------------
    g4a = load(os.path.join(run, "S1_g4a", "g4a_report.json"))
    if g4a:
        md.append("## S1: gate G4a (C-05)")
        md.append("")
        md.append("G4a all pass: **%s** (T-01 %s, T-02 %s, T-03 %s, T-04 %s, T-05 %s, T-06 %s)." % tuple(
            [g4a.get("G4a_all_pass")] + [g4a.get("T0%d_pass" % i) for i in range(1, 7)]))
        md.append("")
        md.append("| Field | Value |")
        md.append("|---|---|")
        for k in ("T01_maxIm_frozen", "T02_res_psi0", "T02_newton_res", "T03_true_maxIm", "T03_true_small", "T03_regression_psi0",
                  "T04_loss_rate", "T04_drift", "T04_drift_grid", "T04_Q_ratio", "T04_loss_per_hawking_period", "T06_closed_EN_max",
                  "T06_closed_nu_max", "T06_closed_nbar_peak", "T06_closed_nbar_band"):
            md.append("| %s | %s |" % (k, fmt(g4a.get(k))))
        for k in ("eta15", "eta_mean_running", "DeltaS_xi_ge_2", "eta_xi_ge_2"):
            md.append("| T05.%s | %s |" % (k, fmt((g4a.get("T05") or {}).get(k))))
        md.append("")
        for k, v in (g4a.get("T06_tokens") or {}).items():
            tokens["[" + k.replace("_", "-") + "]"] = {"value": v, "source": "g4a_report.json T06_tokens (closed-form kappa)"}
        tokens["[G4A-NEWTON-RES]"] = {"value": sci2(mx(g4a.get("T02_newton_res"))),
                                      "source": "g4a_report.json max(T02_newton_res) over delta3 in {0.02, 0.05}, 2 s.f."}
        tokens["[G4A-ZERO-SECTOR]"] = {"value": sci2(mx(g4a.get("T03_true_small")), up=True),
                                       "source": "g4a_report.json max(T03_true_small), 2 s.f., rounded up (bound)"}
        tokens["[G4A-ANNULUS-IM]"] = {"value": sci2(mx(g4a.get("T03_true_maxIm")), up=True),
                                      "source": "g4a_report.json max(T03_true_maxIm), 2 s.f., rounded up (bound)"}

    # ---- S2 ------------------------------------------------------------------------------------------
    g4c = load(os.path.join(run, "S2_g4c", "g4c_loss_report.json"))
    if g4c:
        pts = as_list(g4c.get("points"))
        op = [p for p in pts if p.get("N") == 2 and abs(p.get("d3", 0) - 0.05) < 1e-9]
        op = op[0] if op else {}
        sc = (g4c.get("scaling") or {}).get("d3_050", {})
        md.append("## S2: G4c (C-05b, T-10 + T-12)")
        md.append("")
        md.append("Converged points: %d of %d; exact-scaling ratio %s (within 5%%: %s); recoil deviation at the operating point %s "
                  "(<= 15%%: %s); fast run: %s." % (
                      sum(1 for p in pts if p.get("converged")), len(pts), fmt(sc.get("ratio")),
                      sc.get("ratio") is not None and abs(sc["ratio"] - 1) <= 0.05, fmt(op.get("recoil_rel_dev")),
                      op.get("recoil_rel_dev") is not None and op["recoil_rel_dev"] <= 0.15, g4c.get("fast")))
        md.append("")
        peak = g4c.get("peak_estimator") or "grid"
        md.append("Peak position tau_c: %s." % (
            "sub-grid (vertex of the parabola through ln|psi|^2 at the grid maximum and its two neighbours; SOTHE-P2). "
            "The grid-maximum columns and tokens are the definition of the MATLAB kit (p2_splitstep_loss.m 0.1.1)"
            if peak == "subgrid" else "grid maximum (p2_splitstep_loss.m 0.1.1)"))
        md.append("")
        md.append("| N | delta3 | 1-Q(40)/Q(0) | c | first-HP loss | drift [30,40] | drift [30,40], grid peak | recoil dev | recoil dev, grid peak | rel_dev | converged |")
        md.append("|---|---|---|---|---|---|---|---|---|---|---|")
        for p in pts:
            md.append("| %d | %.3f | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
                p["N"], p["d3"], fmt(p.get("loss_0_40")), fmt(p.get("c_over_xi")), fmt(p.get("loss_first_hp")),
                fmt((p.get("window_drift") or [None] * 5)[4]), fmt((p.get("window_drift_grid") or [None] * 5)[4]),
                fmt(p.get("recoil_rel_dev")), fmt(p.get("recoil_rel_dev_grid")), fmt(p.get("rel_dev")), p.get("converged")))
        md.append("")
        chg = as_list(g4c.get("tokens_changed_by_peak"))
        tg = g4c.get("tokens_grid_peak") or {}
        if tg:
            md.append("Tokens whose text depends on the peak estimator: %s." % (", ".join("`[%s]`" % k.replace("_", "-") for k in chg) or "none"))
            md.append("")
            if chg:
                md.append("| Token | this run (sub-grid peak) | grid-maximum peak (MATLAB 0.1.1 definition) |")
                md.append("|---|---|---|")
                for k in chg:
                    a, b = (g4c.get("tokens") or {}).get(k, ""), tg.get(k, "")
                    if k == "G4C_LOSS_TABLE":
                        def col(t):
                            return " ".join(line.split("&")[3].strip().split("\\")[0].strip() for line in t.splitlines() if line.count("&") == 3 and "N\\dtres" not in line)
                        a, b = "drift column: " + col(a), "drift column: " + col(b)
                    md.append("| `[%s]` | `%s` | `%s` |" % (k.replace("_", "-"), a, b))
                md.append("")
            S["peak_estimator"] = peak
            S["tokens_changed_by_peak"] = {k: {"subgrid": (g4c.get("tokens") or {}).get(k), "grid": tg.get(k)} for k in chg}
            if chg:
                S["decisions"].append("Peak estimator (sub-grid, SOTHE-P2): %s differ from the grid-maximum definition of the MATLAB kit "
                                      "(p2_splitstep_loss.m 0.1.1); a canonical MATLAB run reproduces the grid-maximum values unless the "
                                      "kit adopts the sub-grid peak (PI decision)" % ", ".join("[%s]" % k.replace("_", "-") for k in chg))
        for k, v in (g4c.get("tokens") or {}).items():
            if k == "G4C_LOSS_TABLE":
                tokens["[G4C-LOSS-TABLE]"] = {"value": "(LaTeX array, see g4c_loss_report.json tokens.G4C_LOSS_TABLE)",
                                              "source": "g4c_loss_report.json tokens"}
            else:
                tokens["[" + k.replace("_", "-") + "]"] = {"value": v, "source": "g4c_loss_report.json tokens"}
        for k in as_list(g4c.get("tokens_changed_by_peak")):
            key = "[" + k.replace("_", "-") + "]"
            if key in tokens:
                tokens[key]["source"] += " (sub-grid peak; grid-maximum peak gives %s)" % (
                    "a different drift column" if k == "G4C_LOSS_TABLE" else "`%s`" % (g4c.get("tokens_grid_peak") or {}).get(k))
                tokens[key]["value_grid_peak"] = (g4c.get("tokens_grid_peak") or {}).get(k)

    # ---- S3 ------------------------------------------------------------------------------------------
    g12 = load(os.path.join(run, "S3_g12", "g12_report.json"))
    if g12:
        md.append("## S3: provenance of entropy_trajectory.csv (G-12, T-19)")
        md.append("")
        md.append("| Run | N_sol | delta3 | Nt / steps | DeltaS_tot | eta15 | eta(xi>=2) | DeltaS(xi>=2) | max dS_tot vs archive |")
        md.append("|---|---|---|---|---|---|---|---|---|")
        for r in as_list(g12.get("runs")):
            md.append("| %s | %.1f | %.2f | %d / %d | %.7f | %.4f | %.4f | %.4f | %.2e |" % (
                r["tag"], r["N_sol"], r["delta3"], r["Nt"], r["n_steps"], r["DeltaS_tot"], r["eta15"], r["eta_xi_ge_2"],
                r["DeltaS_xi_ge_2"], r["max_abs_dev_vs_archive"]["S_tot"]))
        md.append("")
        md.append("Verdict: %s" % g12.get("verdict"))
        md.append("")
        md.append("Mask rule (source of `[G4C-MASK-RULE]`): %s" % g12.get("mask_rule"))
        md.append("")
        S["g12"] = {"verdict": g12.get("verdict"), "A_reproduces": g12.get("archive_reproduced_by_A")}
        S["decisions"].append(g12.get("decision_needed") or "G-12/T-19: see g12_report.json")
        tokens["[G4C-MASK-RULE]"] = {"value": "(PI wording, T-19)", "source": "S3 g12_report.json mask_rule; supplement Sec. 2"}

    # ---- S4 ------------------------------------------------------------------------------------------
    o1 = load(os.path.join(run, "S4_oracles", "oracle_ledger_v010.json"))
    o2 = load(os.path.join(run, "S4_oracles", "oracle_ledger_v020.json"))
    if o2:
        rows = as_list(o2.get("rows"))
        md.append("## S4: oracle ledger v0.2.0 (T-14)")
        md.append("")
        md.append("v0.1.0: %s/%s; v0.2.0: %s/%s; retired: %s." % (
            (o1 or {}).get("n_pass"), (o1 or {}).get("n_total"), o2.get("n_pass"), o2.get("n_total"), ", ".join(as_list(o2.get("retired_ids")))))
        md.append("")
        md.append("| ID | Oracle | Class | Change | Pass | Value | Reference |")
        md.append("|---|---|---|---|---|---|---|")
        for r in rows:
            if r.get("change") != "kept" or not r.get("pass"):
                md.append("| %s | %s | %s | %s | %s | %s | %s |" % (r["id"], r["name"], r["class"], r["change"], r["pass"], fmt(r["val"]), fmt(r["ref"])))
        md.append("")
        md.append("(Rows kept unchanged and passing are in oracle_ledger_v020.json.)")
        md.append("")
        S["oracles"] = {"v010": [(o1 or {}).get("n_pass"), (o1 or {}).get("n_total")], "v020": [o2.get("n_pass"), o2.get("n_total")]}
        S["decisions"].append("T-14: record the oracle changes (retired P1.4, P1.5, P3.3 and their pins) in a DELTA note")

    # ---- S5 ------------------------------------------------------------------------------------------
    fm = load(os.path.join(run, "S5_figures", "figs_manifest.json"))
    if fm:
        md.append("## S5: manuscript figures (T-13, T-18)")
        md.append("")
        md.append("RECORDED renders (HR-3) in the publication style of pubstyle.py (IOP guidelines), no watermark. Captions: "
                  "S5_figures/captions.md and captions.tex (explanation only); provenance: figs_manifest.json; plotted arrays: "
                  "figure_data.npz.")
        md.append("")
        md.append("| File | Rendered as | Vector |")
        md.append("|---|---|---|")
        figs = as_list(fm.get("figures"))
        for f in figs:
            md.append("| %s | %s | %s, %s |" % (f.get("file"), f.get("rendered"), f.get("vector"), f.get("eps")))
            d = f.get("data") or {}
            cap = d.get("caption_update_needed")
            if cap:
                S["decisions"].append("T-13 caption of %s (%s): %s (suggested caption text: S5_figures/captions.md)"
                                      % (f.get("rendered"), f.get("file"), cap))
            nout = as_list(d.get("n_outside_ylim"))
            if nout and any(nout):
                S["decisions"].append("%s range: %s eigenvalue(s) of panel (a) and %s of panel (b) lie outside the vertical range "
                                      "|Im Omega| <= %g that make_paper2_figs.m sets (the same in MATLAB); the caption says so. "
                                      "Widening the range (in both scripts) is a PI decision"
                                      % (f.get("rendered"), nout[0], nout[1] if len(nout) > 1 else 0, d.get("ylim")))
        md.append("")
        S["figures"] = [f.get("file") for f in figs]

    # ---- S6 ------------------------------------------------------------------------------------------
    st6 = load(os.path.join(run, "S6_phase2", "stage_status.json"))
    ps = load(os.path.join(run, "S6_phase2", "paper2_data", "phase2_summary.json"))
    if ps:
        c6 = (st6 or {}).get("checks") or {}
        md.append("## S6: phase-2 data products (run_phase2_all v0.1.0)")
        md.append("")
        md.append("Oracles v0.1.0 on the full-grid pass: %s/%s; CSV files %s of 10; figures %s of 6; stability map %s at n_tau = %s; "
                  "bundle: %s." % (ps.get("oracles_passed"), ps.get("oracles_total"), c6.get("csv_written"), c6.get("figures_written"),
                                   " x ".join(str(x) for x in (c6.get("stab_grid") or [])), c6.get("stab_n_tau"),
                                   ", ".join(c6.get("bundle") or [])))
        md.append("")
        md.append("| Quantity | Value |")
        md.append("|---|---|")
        for k in ("kappa_kinematic", "T_H", "thermality_slope", "kappa_fit", "flux_residual", "dS_tot", "eta_GSL", "E_N_max_nats",
                  "E_N_at_3kappa", "PT_nu_minus_max", "kappa_drive_spread", "kappa_cont_kg0", "greybody_range",
                  "bdg_gap_over_mu_d3_0", "bdg_maxIm_d3_0", "bdg_maxIm_d3_002", "radiation_rate_table"):
            md.append("| %s | %s |" % (k, fmt(ps.get(k))))
        md.append("")
        pv = load(os.path.join(run, "S6_phase2", "parity_vs_matlab.json"))
        if pv:
            md.append("Parity with the MATLAB R2025b products (reference/matlab_R2025b_U1): %d files compared; byte-identical: %s; "
                      "worst absolute deviation %s." % (pv.get("files_compared", 0), ", ".join(pv.get("byte_identical") or []) or "none",
                                                        fmt(pv.get("worst_abs_dev"))))
            md.append("")
            md.append("| File | Byte-identical | Worst column (max abs / max rel) |")
            md.append("|---|---|---|")
            for name, r in (pv.get("files") or {}).items():
                if not r.get("compared"):
                    continue
                if "eig_set_distance" in r:
                    ka = r.get("krein_counts_annulus") or {}
                    worst = "eigenvalue-set distance %s (annulus %s); Krein +/- in the annulus %s vs %s" % (
                        fmt(r["eig_set_distance"]), fmt(r.get("eig_set_distance_annulus")), ka.get("this"), ka.get("reference"))
                else:
                    cols = r.get("columns") or {}
                    if cols:
                        k, v = max(cols.items(), key=lambda kv: kv[1]["max_abs"])
                        worst = "%s: %s / %s" % (k, fmt(v["max_abs"]), fmt(v["max_rel"]))
                    else:
                        worst = "shape differs: rows %s" % r.get("rows")
                md.append("| %s | %s | %s |" % (name, r.get("byte_identical"), worst))
            md.append("")
        md.append("Caveat: the delta3 > 0 BdG rates of these products (radiation_rate_table.csv, bdg_spectrum_d3_020.csv, "
                  "stability_map.csv, fig2 (b), fig5) linearize about psi0 = N sech(N tau), which is not stationary for delta3 != 0 "
                  "(G3u audit, S1-1; oracles P1.4, P1.5 retired in v0.2.0). They are regression products of suite v0.1.0.")
        md.append("")
        S["phase2"] = {"oracles": [ps.get("oracles_passed"), ps.get("oracles_total")], "csv": c6.get("csv_written"),
                       "figures": c6.get("figures_written"), "parity_worst_abs_dev": (pv or {}).get("worst_abs_dev"),
                       "byte_identical": (pv or {}).get("byte_identical")}

    # ---- S7 ------------------------------------------------------------------------------------------
    sw = load(os.path.join(run, "S7_sweep", "robustness_sweep.json"))
    if sw:
        cl = sw.get("claims_of_ref25") or {}
        md.append("## S7: Ref. [25] robustness sweep (run_robustness_sweep)")
        md.append("")
        st7 = sw.get("settings") or {}
        md.append("Settings: Nt = %s, %s steps, xi_max = %s, %s saves, window %s%s; %d configurations in %.0f s." % (
            st7.get("Nt"), st7.get("n_steps"), st7.get("xi_max"), st7.get("n_save"), st7.get("tau_window"),
            " (FAST)" if st7.get("fast") else "", len(sw.get("rows") or []), sw.get("runtime_s") or 0))
        md.append("")
        md.append("| | This run | Ref. [25] |")
        md.append("|---|---|---|")
        md.append("| DeltaS_tot range (nats) | %s | %s |" % (fmt(sw.get("DeltaS_range")), fmt(cl.get("DeltaS_range"))))
        md.append("| eta_GSL range | %s | %s |" % (fmt(sw.get("eta_range")), fmt(cl.get("eta_range"))))
        md.append("| DeltaS > 0 and eta > 1/2 | %s of %d (%.1f%%) | all |" % (sw.get("pass_count"), len(sw.get("rows") or []),
                                                                            100 * (sw.get("pass_rate") or 0)))
        md.append("")
        md.append("| Subset | DeltaS range | eta range | Pass |")
        md.append("|---|---|---|---|")
        for grp, d in (sw.get("by") or {}).items():
            for name, v in d.items():
                md.append("| %s = %s | %s | %s | %d/%d |" % (grp, name, fmt(v["DeltaS_range"]), fmt(v["eta_range"]), v["pass"], v["n"]))
        md.append("")
        md.append("Note: %s." % sw.get("note"))
        md.append("")
        ck = load(os.path.join(run, "S7_sweep", "ref25_checks.json"))
        if ck:
            n1, pm, dm = ck.get("nominal") or {}, ck.get("phase_matching") or {}, ck.get("dimensional") or {}
            md.append("Other Ref. [25] functions (protocol of SOTHE_pkg validation/validate_paper_claims.py%s):" % (", FAST" if ck.get("fast") else ""))
            md.append("")
            md.append("| Check | Result | Ref. [25] / criterion | Pass |")
            md.append("|---|---|---|---|")
            md.append("| nominal run N = 3.5, delta3 = 0.02 | DeltaS_tot %s, eta %s, dS<0 fraction %s, photon drift %s | DeltaS > 0, eta > 1/2, "
                      "fraction in (0.30, 0.65), drift < 1e-3 | %s |" % (fmt(n1.get("Delta_S_tot")), fmt(n1.get("eta_GSL")),
                                                                         fmt(n1.get("neg_frac")), fmt(n1.get("drift")), n1.get("passed")))
            md.append("| phase matching, delta3 = 0.06 ... 0.10 | c_RR %s, R^2 %s, mean k_RR delta3 %s, spread %s | 1.01, 0.93, 1.02, < 0.05; "
                      "criterion 0.8 < mean < 1.2 and R^2 > 0.5 | %s |" % (fmt(pm.get("c_RR")), fmt(pm.get("R2")), fmt(pm.get("invariant_mean")),
                                                                            fmt(pm.get("invariant_spread")), pm.get("passed")))
            md.append("| gnlse_dimensional at its defaults | DeltaS_tot %s, eta %s, photon drift %s | (no reference) | ran |" % (
                fmt(dm.get("Delta_S_tot")), fmt(dm.get("eta_GSL")), fmt(dm.get("photon_drift"))))
            md.append("")
            S["ref25_checks"] = {"nominal": n1.get("passed"), "phase_matching": pm.get("passed"), "c_RR": pm.get("c_RR")}
        S["sweep"] = {"DeltaS_range": sw.get("DeltaS_range"), "eta_range": sw.get("eta_range"), "pass_count": sw.get("pass_count"),
                      "n": len(sw.get("rows") or []), "settings": st7}
        def outside(r, c):
            return r is not None and c is not None and (r[0] < c[0] - 0.005 or r[1] > c[1] + 0.005)

        n7 = len(sw.get("rows") or [])
        fails = (sw.get("pass_count") or 0) < n7
        off = outside(sw.get("DeltaS_range"), cl.get("DeltaS_range")) or outside(sw.get("eta_range"), cl.get("eta_range"))
        if not st7.get("fast") and (fails or off):
            S["decisions"].append("S7: the sweep differs from the Ref. [25] statement (DeltaS in [0.82, 2.42] nats, eta in "
                                  "[0.52, 0.62], every configuration passing): DeltaS %s, eta %s, %s of %d pass; "
                                  "see S7_sweep/robustness_sweep.json" % (fmt(sw.get("DeltaS_range")), fmt(sw.get("eta_range")),
                                                                          sw.get("pass_count"), n7))

    # ---- raw data kept for analysis --------------------------------------------------------------------
    data = []
    for sid, d, title in STAGES:
        p = os.path.join(run, d)
        if os.path.isdir(p):
            for dp, dn, fn in os.walk(p):
                dn.sort()
                for n in sorted(fn):
                    if n.endswith((".npz", ".csv", ".mat")) and n != "entropy_trajectory.csv":
                        q = os.path.join(dp, n)
                        data.append((os.path.relpath(q, run).replace(os.sep, "/"), os.path.getsize(q)))
    if data:
        md.append("## Raw data kept for analysis")
        md.append("")
        md.append("| File | Size (kB) |")
        md.append("|---|---|")
        for rel, size in data:
            md.append("| %s | %.0f |" % (rel, size / 1024))
        md.append("")
        S["data_files"] = [rel for rel, _ in data]

    # ---- tokens and decisions --------------------------------------------------------------------------
    tokens["[RELEASE-TAG]"] = {"value": "(PI, B3)", "source": "repository release"}
    tokens["[ZENODO-VERSION-DOI]"] = {"value": "(PI, B3)", "source": "Zenodo version DOI"}
    for v in tokens.values():
        v["status"] = "RECORDED (HR-3)"
    md.append("## Token values supplied by this run (RECORDED; no manuscript was edited)")
    md.append("")
    md.append("| Token | Value | Source |")
    md.append("|---|---|---|")
    for k in sorted(tokens):
        md.append("| `%s` | `%s` | %s |" % (k, tokens[k]["value"], tokens[k]["source"]))
    md.append("")
    md.append("Tokens listed: %d of 24." % len(tokens))
    md.append("")
    S["decisions"].insert(0, "HR-3: this run is RECORDED evidence only; the canonical question is open (PI)")
    md.append("## For the PI")
    md.append("")
    for d in S["decisions"]:
        md.append("- " + d)
    md.append("")
    with open(os.path.join(run, "SUMMARY.md"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(md) + "\n")
    with open(os.path.join(run, "summary.json"), "w", encoding="utf-8", newline="\n") as f:
        json.dump(S, f, indent=1)
    with open(os.path.join(run, "tokens_proposed.json"), "w", encoding="utf-8", newline="\n") as f:
        json.dump(tokens, f, indent=1, sort_keys=True)
    print("wrote %s (%d tokens; RECORDED)" % (os.path.join(run, "SUMMARY.md"), len(tokens)))


if __name__ == "__main__":
    main()
