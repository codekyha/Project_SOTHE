"""SUMMARY.md of a sweep run: provenance, oracle ledger, block-by-block findings mapped to the REVIEW_G5 items."""
import json
import os

from . import PACKAGE, __version__
from .grid import BLOCKS
from .util import RULE, fnum, read_json, utc

REVIEW_ITEMS = {
    "B1": "S1-1 option A (re-anchor), W11; kinematic layer and Eq. (7'') over the grid",
    "B2": "S1-1 (E_N on the Hawking window), R-02/R-02A, R-10 (loss-robust E_N)",
    "B3": "S1-1 A (loss at N = 1), T-10/T-12 inputs of the census, D-2",
    "B4": "S1-1 (escape threshold, window table), R-02, R-08 (census figure), W11",
    "B5": "A-2 zero-sector rank tests, BdG blocks at every N (S1-1 A unifies them at N = 1)",
    "B6": "S1-3 (Raman and self-steepening in silica), D-4 option b (Paper-3 scope)",
    "B7": "S1-2 / D-06 (scattering solve: protocol, numbers, checks), D-3; INV5: not a Paper-2 claim",
    "B8": "W2, FINDING S7, OD-P8 (Ref. [25] ranges)",
    "B9": "R-11 (the parity-table script), Supplement Tables 2 and 3",
    "B10": "R-02A (anchor at (N, k_op) = (1, 1.5)): paper2-style CSV files and the anchor numbers",
}


def _v(x):
    if isinstance(x, float):
        return fnum(x, "%.6g")
    if isinstance(x, (list, tuple)):
        return "[" + ", ".join(_v(y) for y in x) + "]"
    if isinstance(x, dict):
        return "{" + ", ".join("%s: %s" % (k, _v(v)) for k, v in x.items()) + "}"
    return str(x)


def _cell(x, n=90):
    s = _v(x).replace("|", "/")
    return s if len(s) <= n else s[:n - 3] + "..."


def write_summary(rundir, cfg, blocks, S, O, failed, errors, tm, wall_ab):
    prov = (read_json(os.path.join(rundir, "provenance.json")) or {}).get("current", {})
    L = []
    L.append("# %s results: run %s" % (PACKAGE, os.path.basename(rundir)))
    L.append("")
    L.append("**Status.** %s" % RULE)
    L.append("")
    L.append("## Run")
    L.append("")
    L.append("| Item | Value |")
    L.append("|---|---|")
    for k, v in (("package", "%s %s" % (PACKAGE, __version__)), ("written (UTC)", utc()), ("host", prov.get("host")),
                 ("SLURM job", prov.get("slurm_job_id") or "local"), ("workers", prov.get("workers")),
                 ("versions", ", ".join("%s %s" % kv for kv in (prov.get("versions") or {}).items())),
                 ("BLAS", ", ".join("%s: %s" % kv for kv in (prov.get("blas") or {}).items()) or "n/a"),
                 ("mode", "FAST (smoke test)" if cfg["fast"] else "full sweep"), ("blocks", " ".join(blocks)),
                 ("phase A+B wall time", "%.0f s" % wall_ab), ("failed tasks", len(failed))):
        L.append("| %s | %s |" % (k, _cell(v, 200)))
    L.append("")
    L.append("## Oracle ledger")
    L.append("")
    cnt = {c: [sum(1 for r in O.rows if r["class"] == c and r["pass"] is t) for t in (True, False, None)] for c in ("hard", "claim", "info")}
    L.append("| Class | pass | fail | no rule / not evaluated |")
    L.append("|---|---|---|---|")
    for c in ("hard", "claim", "info"):
        L.append("| %s | %d | %d | %d |" % (c, cnt[c][0], cnt[c][1], cnt[c][2]))
    L.append("")
    L.append("`hard`: identities, reproductions of recorded runs, exact checks (a failure means the pipeline is wrong). "
             "`claim`: a statement of the manuscript or of REVIEW_G5 put to the test (a failure is a finding). `info`: tracked, no rule.")
    L.append("")
    L.append("| Oracle | Block | Class | Pass | Value | Reference | Statement |")
    L.append("|---|---|---|---|---|---|---|")
    for r in O.rows:
        p = {True: "PASS", False: "**FAIL**", None: "-"}[r["pass"]]
        L.append("| %s | %s | %s | %s | %s | %s | %s |" % (r["id"], r["block"], r["class"], p, _cell(r["value"], 70), _cell(r["reference"], 40),
                                                         _cell(r["description"] + ((" (" + r["note"] + ")") if r["note"] else ""), 260)))
    L.append("")
    L.append("## Blocks")
    L.append("")
    for b in blocks:
        L.append("### %s: %s" % (b, BLOCKS[b]))
        L.append("")
        L.append("Review items: %s." % REVIEW_ITEMS[b])
        L.append("")
        L += _block_lines(b, S, cfg)
        t = tm.get(b)
        if t:
            L.append("- Tasks: %d, %.0f core-s in total, longest %.0f s." % (t["n"], t["wall_s"], t["max_s"]))
        L.append("")
    if failed:
        L.append("## Failed tasks")
        L.append("")
        for f in failed:
            L.append("- `%s` (traceback in `tasks/<block>/%s.error.json`; `resume` retries it)" % (f, f))
        L.append("")
    if errors:
        L.append("## Phase-C errors")
        L.append("")
        for k, v in errors.items():
            L.append("- %s: `%s`" % (k, v.strip().splitlines()[-1][:300]))
        L.append("")
    L.append("## Files")
    L.append("")
    L.append("- `ORACLES.json` (ledger), `summary.json` (all reduced numbers), `config.json` (the grids), `provenance.json`, `plan_A.json`, `plan_B.json`.")
    L.append("- `B1_kinematics/` ... `B10_anchor/`: CSV tables per block (floats at %.17g); `tasks/`: the raw result of every task (JSON).")
    L.append("- `figures/`: PDF + EPS (vector) + PNG (600 dpi), `captions.md`, `captions.tex`.")
    L.append("- `MANIFEST.sha256`: SHA-256 of every file of the run folder except the live logs.")
    L.append("")
    with open(os.path.join(rundir, "SUMMARY.md"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(L) + "\n")


def _block_lines(b, S, cfg):
    L = []
    s = S.get(b)
    if not s:
        return ["- No result (block not run or failed)."]
    if b == "B1":
        op = s.get("op") or {}
        L.append("- Grid rows: %d; kappa from %s to %s." % (s["n_rows"], fnum(s["kappa_range"][0]), fnum(s["kappa_range"][1])))
        if op:
            L.append("- Operating point: kappa_num = %s (closed form %s, deficit %s, Eq. (7'') %s); drive spread %s (band %s)." % (
                fnum(op["kappa_num"], "%.13g"), fnum(op["kappa_closed"], "%.10g"), fnum(op["deficit_rel"], "%.4e"), fnum(op["eq7pp_rel"], "%.4e"),
                fnum(op["drive_spread_abs"], "%.4e"), fnum(op["band_abs"], "%.4e")))
        L.append("- Largest |Eq. (7'') residual| / |deficit| over the grid: %s; largest spread / band: %s." % (fnum(s["eq7pp_rel_max"]),
                                                                                                         fnum(s["spread_over_band_max"])))
    elif b == "B9":
        L.append("- Table 2 max deviation from the published values: %s; Table 3 max relative deviation: %s." % (fnum(s["table2_max_dev"]),
                                                                                                          fnum(s["table3_max_rel"])))
    elif b == "B3":
        L.append("- Points: %d; loss over [0, 40] from %s to %s." % (s["n_points"], fnum(s["loss_0_40_range"][0]), fnum(s["loss_0_40_range"][1])))
        if s.get("op"):
            L.append("- Operating point (2, 0.05): loss %s, drift %s -> %s." % (fnum(s["op"]["loss_0_40"]), fnum(s["op"]["window_drift"][0]),
                                                                               fnum(s["op"]["window_drift"][4])))
        if s.get("anchor"):
            L.append("- (1, 0.05): loss %s, one Hawking period after xi = 10 (k_op = 1): %s, drift %s -> %s." % (
                fnum(s["anchor"]["loss_0_40"]), fnum(s["anchor"]["loss_first_hp"]), fnum(s["anchor"]["window_drift"][0]),
                fnum(s["anchor"]["window_drift"][4])))
        if s.get("reference"):
            L.append("- Job 530047 points reproduced: %d (max relative loss deviation %s, max drift deviation %s)." % (
                len(s["reference"]), fnum(max(max(c["loss_0_40_rel"], c["loss_rate_rel"], c["loss_first_hp_rel"]) for c in s["reference"])),
                fnum(max(c["window_drift_abs"] for c in s["reference"]))))
    elif b == "B4":
        L.append("- Operating point: pair channel closed until xi = %s; open below omega = %s in [30, 40]." % (s["op_pair_closed_until"],
                                                                                                          fnum(s["op_pair_open_below_late"])))
        L.append("- Largest grid N with a Hawking window: %s (bound sqrt(2 x 1.6) = 1.789)." % s["largest_N_with_window"])
        L.append("")
        L.append("| N | k_op | mu | pi N/S | n(mu) | E_N(mu) | lab edge | lab window | [30,40] V_d | [30,40] edge | [30,40] window |")
        L.append("|---|---|---|---|---|---|---|---|---|---|---|")
        for t in s["window_table"]:
            if t["k_op"] not in (1.0, 1.5, 2.0, 3.0):
                continue
            lw = "none" if t["lab_window"][0] is None else "[%s, %s] (%.0f%%)" % (fnum(t["lab_window"][0], "%.4g"), fnum(t["lab_window"][1], "%.4g"),
                                                                                 100 * t["lab_frac"])
            cw = "n/a" if t["late_window"] is None else ("none" if t["late_window"][0] is None else "[%s, %s] (%.0f%%)" % (
                fnum(t["late_window"][0], "%.4g"), fnum(t["late_window"][1], "%.4g"), 100 * t["late_frac"]))
            L.append("| %g | %g | %s | %s | %s | %s | %s | %s | %s | %s | %s |" % (
                t["N"], t["k_op"], fnum(t["mu"], "%.3g"), fnum(t["piN_over_S"], "%.4g"), fnum(t["n_at_mu"], "%.3g"), fnum(t["EN_at_mu"], "%.3g"),
                fnum(t["lab_edge"], "%.4g"), lw, fnum(t["late_Vd"], "%.4g"), fnum(t["late_edge"], "%.4g"), cw))
        L.append("")
    elif b == "B2":
        for k, v in s["curves"].items():
            L.append("- %s: E_N at omega = 0.05 %s nats, nu_- at 1.6 %s." % (k, fnum(v["EN_max"]), fnum(v["nu_max"])))
    elif b == "B5":
        L.append("")
        L.append("| N | delta3 | N delta3 | Newton | tail | annulus max Im | zero sector | dim ker M, M^2 | psi0 artefact max Im (/sqrt delta3) |")
        L.append("|---|---|---|---|---|---|---|---|---|")
        for r in s["rows"]:
            L.append("| %g | %g | %g | %s (%s) | %s | %s | %s | %s, %s | %s (%s) |" % (
                r["N"], r["d3"], r["N"] * r["d3"], "ok" if r["newton_converged"] else "NOT converged", fnum(r["newton_res"], "%.1e"),
                fnum(r["tail_rel"], "%.2e"), fnum(r["annulus_maxIm"], "%.1e"), fnum(r["zero_sector_max"], "%.1e"), r["dim_ker_M"], r["dim_ker_M2"],
                fnum(r.get("psi0_annulus_maxIm"), "%.5f"), fnum(r.get("psi0_annulus_over_sqrt_d3"), "%.4f")))
        L.append("")
        L.append("The manuscript's BdG scope is N delta3 <= 0.05; rows beyond it are reported, not claimed.")
        L.append("")
    elif b == "B6":
        L.append("- Raman models: Blow-Wood (f_R = 0.18, tau1 = 12.2 fs, tau2 = 32 fs; T_R = f_R M1 = %s fs), linear T_R = %s fs; "
                 "self-steepening s = lambda0 / (2 pi c T0) at 1550 nm." % (fnum(s["T_R_BW_fs"], "%.4g"), fnum(s["T_R_lin_fs"], "%.3g")))
        for g in s["gordon"]:
            L.append("- Gordon check, %s: rate %s against %s (ratio %s)." % (g["label"], fnum(g["rate_full"]), fnum(g["gordon"]), fnum(g["ratio"], "%.4f")))
        L.append("")
        L.append("| N | model | T0 (fs) | Gordon rate | initial rate, full field (core) | mean rate [0, 5] | dk over 1st Hawking period | total shift "
                 "(linear extrapolation) | x recoil | drift [30,40] |")
        L.append("|---|---|---|---|---|---|---|---|---|---|")
        for t in s["S1_3_table"]:
            L.append("| %g | %s | %g | %s | %s (%s) | %s | %s | %s (%s) | %s | %s |" % (
                t["N"], t["model"], t["T0_fs"], fnum(t["gordon_rate"], "%.4g"), fnum(t["initial_rate_kfull"], "%.4g"), fnum(t["initial_rate_ks"], "%.4g"),
                fnum(t["early_rate_ks"], "%.4g"),
                fnum(t["dk_first_HP"], "%.3g"), fnum(t["shift_total_ks"], "%.3g"), fnum(t["linear_extrapolation"], "%.3g"),
                fnum(t["times_recoil"], "%.3g"), fnum(t["late_drift"], "%.4g")))
        L.append("")
    elif b == "B7":
        f = s.get("findings")
        if f:
            L.append("- FINDINGS check (N = 2, k_op = 3): r_conf = %s at omega = %s (recorded %s); max relative deviation %s." % (
                ", ".join(fnum(x["r_conf"], "%.5f") for x in f["rows"]), ", ".join("%g" % x["omega"] for x in f["rows"]),
                ", ".join("%g" % x["findings"] for x in f["rows"]), fnum(f["max_rel_dev"], "%.2e")))
        a = s.get("anchor")
        if a:
            L.append("- (N, k_op) = (1, 1.5) at delta3 = 0: two-sided over %.0f%% of the band; n_H / Planck from %s to %s; "
                     "fit of ln n_H gives kappa = %s against kappa0 = %s." % (100 * a["frac_two_sided"], fnum(a["n_H_over_planck_min"], "%.3g"),
                                                                          fnum(a["n_H_over_planck_max"], "%.3g"),
                                                                          fnum((a["thermal_fit_n_H"] or {}).get("kappa_fit"), "%.4g"),
                                                                          fnum(a["kappa0"], "%.4g")))
        L.append("- Two-sided points (upstream u channel open somewhere in the band): %d of %d." % (len(s["two_sided_points"]), s["n_points"]))
        L.append("- INV5: none of these numbers is a Paper-2 claim; they document D-06.")
    elif b == "B8":
        L.append("")
        L.append("| Variant | configs | DeltaS range | eta range | all DeltaS>0, eta>1/2 | outside [0.82, 2.42] | outside [0.52, 0.62] |")
        L.append("|---|---|---|---|---|---|---|")
        for v in s["variants"]:
            L.append("| %s | %d | [%s, %s] | [%s, %s] | %d/%d | %d | %d |" % (v["variant"], v["n_configs"], fnum(v["DeltaS_range"][0], "%.4f"),
                                                                          fnum(v["DeltaS_range"][1], "%.4f"), fnum(v["eta_range"][0], "%.4f"),
                                                                          fnum(v["eta_range"][1], "%.4f"), v["pass_count"], v["n_configs"],
                                                                          v["n_DeltaS_below"] + v["n_DeltaS_above"], v["n_eta_below"] + v["n_eta_above"]))
        L.append("")
    elif b == "B10":
        L.append("- kappa = %s (closed form %s), T_H = %s, mu = %s, Hawking period %s; E_N(mu) = %s, E_N(1.6) = %s, E_N(0.05) = %s." % (
            fnum(s["kappa_kin"], "%.6g"), fnum(s["kappa_closed"], "%.6g"), fnum(s["T_H"], "%.5g"), fnum(s["mu"]), fnum(s["hawking_period"], "%.5g"),
            fnum(s["EN_at_mu"], "%.4g"), fnum(s["EN_at_1p6"], "%.4g"), fnum(s["EN_band_max"], "%.4g")))
        if s.get("loss"):
            L.append("- Loss at (1, 0.05): %s over [0, 40]; %s over one anchor Hawking period after xi = 10." % (
                fnum(s["loss"]["loss_0_40"], "%.3g"), fnum(s["loss"]["loss_one_anchor_HP_after_10"], "%.3g")))
        for c in s.get("census", []):
            L.append("- Census %s%s: V_d = %s, Phi_pair = %s, window %s." % (c["frame"], (" " + str(c["epoch"])) if c["epoch"] else "",
                                                                            fnum(c["Vd"], "%.4g"), fnum(c["Phi_pair"], "%.3f"), _v(c["window"])))
        L.append("- CSV files in the paper2 format: `B10_anchor/paper2_data_anchor/`.")
    return L
