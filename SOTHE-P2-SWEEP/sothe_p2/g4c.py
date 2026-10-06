"""Stage S2 (C-05b, gate G4c tasks T-10 and T-12): port of run_g4c_loss_scan.m 0.1.2.

Writes g4c_loss_report.json (points, scaling, fit_diagnostic, channels, tokens, protocol, wall_s,
runtime, fast, script_version) into the working folder.  The 27 propagations (9 points x base, dxi/2,
2n) and the 9 periodic Newton solutions are independent and run in parallel.  The manuscript tokens
are built with the same format strings as the MATLAB script.  RECORDED only (HR-3).

Sub-grid peak (introduced in P2_pypack 2.2.0): the drifts, the recoil check and every token that depends on them use it
(kit.splitstep_loss, peak="subgrid").  The same tokens with the grid-maximum peak of the MATLAB source
(0.1.1) are written as tokens_grid_peak, and tokens_changed_by_peak lists the ones whose text differs."""
import math
import os
import time

import numpy as np

from . import jsonio, kit
from . import compat as C
from .g4a import runtime_string
from .par import pmap

PTS = [(1, 0.02), (1, 0.03), (1, 0.04), (1, 0.05), (1, 0.06), (1, 0.08), (1, 0.10), (2, 0.025), (2, 0.05)]
KG = 0.3
KOP = 1.0
PROTOCOL = ("base xi_end=40 dxi=0.002 n=8192 L=400, launched from psi0; conv: dxi/2 and 2n; windows [0,5],[5,10],"
            "[10,20],[20,30],[30,40]; c = mean(rate*xi_mid) over the last three; kappa = N S(k_op=1) at kg=0.3; "
            "census band [0.05,1.6], 311 points; peak position tau_c: vertex of the parabola through ln|psi|^2 at the "
            "grid maximum and its two neighbours (SOTHE-P2; core window, Q and k_s anchored at the grid maximum); "
            "tokens_grid_peak: the same tokens with tau_c = tau[argmax|psi|^2] (p2_splitstep_loss.m 0.1.1)")
PEAK_ESTIMATOR = "subgrid"


def _task(kind, N, d3):
    if kind == "base":
        X = kit.splitstep_loss(N, d3, 40.0, 0.002, 8192, 400.0, 0.0)
        return X
    if kind == "half":
        X = kit.splitstep_loss(N, d3, 40.0, 0.001, 8192, 400.0, 0.0)
        X.update(Q0=float(X["Q"][0]), Qend=float(X["Q"][-1]))
        return X
    if kind == "n2":
        X = kit.splitstep_loss(N, d3, 40.0, 0.002, 16384, 400.0, 0.0)
        X.update(Q0=float(X["Q"][0]), Qend=float(X["Q"][-1]))
        return X
    if kind == "tod":
        try:
            Sk = kit.tod_soliton(N, d3, 500, 16.0)
            tail = float(np.max(np.abs(Sk["phi"][np.abs(Sk["tau"]) > 12.8])) / N)
            return {"tail": tail, "kmean": Sk["kmean"], "converged": Sk["converged"], "tau": Sk["tau"], "phi": Sk["phi"],
                    "newton_history": np.array(Sk["res"])}
        except Exception:                               # MATLAB: try ... catch -> NaN, false
            return {"tail": float("nan"), "kmean": float("nan"), "converged": False}
    raise ValueError(kind)


def _polyfit1(x, y):
    p = np.polyfit(x, y, 1)
    return float(p[0]), float(p[1])


def _census_and_tokens(rows, scaling, dk="window_drift", rk="recoil_rel_dev"):
    """T-12 channel census in the comoving frame and the manuscript tokens (MATLAB format strings), with the
    window drifts rows[*][dk] and the recoil deviation rows[*][rk]; dk/rk select the peak estimator
    ("window_drift", "recoil_rel_dev": sub-grid; "window_drift_grid", "recoil_rel_dev_grid": grid maximum)."""
    # ---- T-12: channel census in the frame comoving with the soliton (OD-P5 = A)
    op = [r for r in rows if r["N"] == 2 and abs(r["d3"] - 0.05) < 1e-12][-1]
    ch = {"lab_frame": kit.channel_census(2, KOP, 0.05, 0.0, KG)}
    ep = []
    closed_until = 40
    for iw in range(5):
        e = kit.channel_census(2, KOP, 0.05, float(op[dk][iw]), KG)
        e["window"] = [op["window_lo"][iw], op["window_hi"][iw]]
        ep.append(e)
        if e["Phi_pair"] > 0 and closed_until == 40:
            closed_until = op["window_lo"][iw]
    ch["epochs"] = ep
    ch["pair_closed_until"] = closed_until
    ch["pair_open_below_late"] = ep[4]["pair_open_below"]
    xs = np.array([r["d3"] for r in rows if r["N"] == 1])
    ys = np.array([r[dk][4] for r in rows if r["N"] == 1])
    o = np.argsort(xs, kind="stable")
    xs, ys = xs[o], ys[o]
    grid = {}
    for Ns in (0.6, 1.0, 1.4, 2.0):
        xg = 0.05 * Ns
        Vd = Ns * (-xg) if xg < xs[0] else Ns * C.interp1_extrap(xs, ys, xg)
        for kop in (1.0, 1.5, 2.0, 3.0):
            grid["N%02d_k%02d" % (int(C.mround(10 * Ns)), int(C.mround(10 * kop)))] = {
                "lab": kit.channel_census(Ns, kop, 0.05, 0.0, KG), "comoving": kit.channel_census(Ns, kop, 0.05, Vd, KG)}
    ch["grid"] = grid
    g15 = grid["N20_k15"]
    # ---- manuscript tokens (formatted as they enter the text)
    T = {}
    T["G4C_LOSS_C"] = kit.fmt_sci(op["c_over_xi"], 1)
    T["G4C_LOSS_0_40"] = C.mfmt(r"%.1f\%%", 100 * op["loss_0_40"])
    T["G4C_LOSS_FIRST_HP"] = kit.fmt_sci(op["loss_first_hp"], 1)
    T["G4C_DRIFT_EARLY"] = C.mfmt("%.2f", op[dk][0])
    T["G4C_DRIFT_LATE"] = C.mfmt("%.2f", op[dk][4])
    rd = [r["rel_dev_0_40"] for r in rows if not math.isnan(r["rel_dev_0_40"])]
    T["G4C_LOSS_CONV"] = "n/a (G4C_FAST run)" if not rd else kit.fmt_sci(max(rd), 1)
    T["G4C_OMEGA_LATE"] = C.mfmt("%.2f", ch["pair_open_below_late"])
    T["G4C_RECOIL_DEV"] = C.mfmt(r"%.0f\%%", 100 * op[rk])
    T["G4C_CLOSED_UNTIL"] = "%d" % int(C.mround(closed_until))
    tab = ("\\renewcommand{\\arraystretch}{1.15}\n\\begin{array}{c@{\\qquad}c@{\\qquad}l@{\\qquad}r}\n"
           "N\\dtres & N & 1-Q(40)/Q(0) & d\\tau_c/d\\xi\\\\[1pt]\n\\hline\n")
    order = sorted(range(len(PTS)), key=lambda i: (PTS[i][0] * PTS[i][1], PTS[i][0]))
    first = True
    for i in order:
        r = rows[i]
        rl = r"\rule{0pt}{2.4ex}" if first else ""
        first = False
        tab += C.mfmt("%.2f & %d & %s & %.3f%s\\\\\n", r["Nd3"], r["N"], kit.fmt_sci(r["loss_0_40"], 1), r[dk][4], rl)
    T["G4C_LOSS_TABLE"] = tab + "\\end{array}"
    sc = scaling["d3_050"]
    T["G4C_SCALING_SENTENCE"] = (
        r"The two amplitudes are related by the exact scaling $\psi(\tau,\xi)=N\Psi(N\tau,N^2\xi)$ with "
        r"$\tilde\dtres=N\dtres$, under which the $N=2$ run at $\dtres=0.05$ is the $N=1$ run at $\dtres=0.10$ continued to $\xi=160$: "
        r"the $N=1$ entry continued with its own late-time law, $c\ln4$ with $c=%s$, predicts $%s$ against the measured $%s$."
        % (kit.fmt_sci(sc["c_N1_at_2d3"], 1), kit.fmt_sci(sc["predicted_N2"], 2), kit.fmt_sci(sc["loss_0_40_N2"], 2)))
    cov = "the whole band" if g15["comoving"]["Phi_pair"] >= 0.995 else C.mfmt(r"$%.0f\%%$ of the band", 100 * g15["comoving"]["Phi_pair"])
    T["G4C_CHANNEL_MAP_SENTENCE"] = C.mfmt(
        r"at $N=2$ the pair channel covers %s at $k_{\rm op}=1.5$ once the drift is included, "
        r"against $%.0f\%%$ in the laboratory frame, so the emitting region begins at lower $k_{\rm op}$ than the laboratory-frame estimate",
        cov, 100 * g15["lab"]["Phi_pair"])
    it05 = [i for i, r in enumerate(rows) if r["N"] == 1 and abs(r["d3"] - 0.05) < 1e-12][-1]
    it10 = [i for i, r in enumerate(rows) if r["N"] == 1 and abs(r["d3"] - 0.10) < 1e-12][-1]
    T["G4C_TAIL_005"] = kit.fmt_sci(rows[it05]["newton_tail_rel"], 1)
    T["G4C_TAIL_010"] = kit.fmt_sci(rows[it10]["newton_tail_rel"], 1)
    return ch, T


def run(workdir, fast=False):
    t0 = time.time()
    jobs = []
    if not fast:
        jobs += [("half", N, d3) for N, d3 in PTS] + [("n2", N, d3) for N, d3 in PTS]
    jobs += [("base", N, d3) for N, d3 in PTS] + [("tod", N, d3) for N, d3 in PTS]
    res = dict(zip(jobs, pmap(_task, jobs, progress="G4c propagations and Newton solutions")))
    rows = []
    for N, d3 in PTS:
        X = res[("base", N, d3)]
        rec = {"N": N, "d3": d3, "Nd3": N * d3, "loss_rate": X["loss_rate"], "drift": X["drift"], "Q_ratio": X["Q_ratio"]}
        beta0 = 0.5 - d3
        S = abs(beta0) / math.sqrt(beta0 ** 2 + KG ** 2)
        rec["kappa"] = N * S
        HP = 2 * math.pi / rec["kappa"]
        rec["hawking_period"] = HP
        W = [(0, 5), (5, 10), (10, 20), (20, 30), (30, 40)]
        wr, wd, wg, wm = np.zeros(5), np.zeros(5), np.zeros(5), np.zeros(5)
        for iw, (lo, hi) in enumerate(W):
            sw = (X["xi"] >= lo) & (X["xi"] <= hi)
            wr[iw] = -_polyfit1(X["xi"][sw], np.log(X["Q"][sw]))[0]
            wd[iw] = _polyfit1(X["xi"][sw], X["tau_c"][sw])[0]
            wg[iw] = _polyfit1(X["xi"][sw], X["tau_c_grid"][sw])[0]
            wm[iw] = 0.5 * (lo + hi)
        rec["window_lo"] = [w[0] for w in W]
        rec["window_hi"] = [w[1] for w in W]
        rec["window_rate"] = wr
        rec["window_drift"] = wd
        rec["window_drift_grid"] = wg
        rec["window_rate_times_xi"] = wr * wm
        rec["c_over_xi"] = float(np.mean(wr[2:5] * wm[2:5]))
        rec["c_over_xi_spread"] = float((np.max(wr[2:5] * wm[2:5]) - np.min(wr[2:5] * wm[2:5])) / rec["c_over_xi"])
        rec["loss_0_40"] = float(1 - X["Q"][-1] / X["Q"][0])
        rec["loss_first_hp"] = float(1 - C.interp1(X["xi"], X["Q"], 10.0 + HP) / C.interp1(X["xi"], X["Q"], 10.0))
        rec["loss_hp_after_36"] = float(1 - C.interp1(X["xi"], X["Q"], min(36.0 + HP, 40.0)) / C.interp1(X["xi"], X["Q"], 36.0))
        ks = X["k_s"]
        Neff = X["Q"] / 2
        vpred = (ks - 3 * d3 * ks ** 2) - d3 * Neff ** 2
        late = (X["xi"] >= 10) & (X["xi"] <= 38)
        rec["k_s_early"] = float(np.mean(ks[X["xi"] <= 5]))
        rec["k_s_late"] = float(np.mean(ks[X["xi"] >= 30]))
        for key, suffix in (("tau_c", ""), ("tau_c_grid", "_grid")):
            dsm = C.conv_same(C.mgradient(X[key], X["xi"]), np.ones(9) / 9)
            rec["recoil_max_dev" + suffix] = float(np.max(np.abs(dsm[late] - vpred[late])))
            rec["recoil_rel_dev" + suffix] = rec["recoil_max_dev" + suffix] / float(np.max(np.abs(dsm[late])))
        rec["drift_leading_order"] = -d3 * N ** 2
        tk = res[("tod", N, d3)]
        rec["newton_tail_rel"] = tk["tail"]
        rec["newton_kmean"] = tk["kmean"]
        rec["newton_converged"] = tk["converged"]
        if not fast:
            Xh, Xn = res[("half", N, d3)], res[("n2", N, d3)]
            rec["loss_rate_dxi_half"] = Xh["loss_rate"]
            rec["loss_rate_n2"] = Xn["loss_rate"]
            rec["rel_dev"] = max(abs(Xh["loss_rate"] - X["loss_rate"]), abs(Xn["loss_rate"] - X["loss_rate"])) / abs(X["loss_rate"])
            rec["loss_0_40_dxi_half"] = 1 - Xh["Qend"] / Xh["Q0"]
            rec["loss_0_40_n2"] = 1 - Xn["Qend"] / Xn["Q0"]
            rec["rel_dev_0_40"] = max(abs(rec["loss_0_40_dxi_half"] - rec["loss_0_40"]),
                                      abs(rec["loss_0_40_n2"] - rec["loss_0_40"])) / rec["loss_0_40"]
            rec["converged"] = rec["rel_dev"] <= 0.10
        else:
            rec["loss_rate_dxi_half"] = float("nan")
            rec["loss_rate_n2"] = float("nan")
            rec["rel_dev"] = float("nan")
            rec["rel_dev_0_40"] = float("nan")
            rec["converged"] = False
        rec["loss_per_hawking_period_windowfit"] = X["loss_rate"] * HP
        rows.append(rec)
        print("N=%d d3=%.3f  loss[0,40] %.4e  c/xi %.3e (spread %.0f%%)  first-HP loss %.3e  drift %.4f -> %.4f "
              "(grid peak %.4f)  k_s %.3f -> %.3f (recoil dev %.3f; grid peak %.3f)  Newton tail %.1e  rel_dev %s conv %d" % (
                  N, d3, rec["loss_0_40"], rec["c_over_xi"], 100 * rec["c_over_xi_spread"], rec["loss_first_hp"],
                  wd[0], wd[4], wg[4], rec["k_s_early"], rec["k_s_late"], rec["recoil_rel_dev"], rec["recoil_rel_dev_grid"],
                  rec["newton_tail_rel"], "%.4g" % rec["rel_dev"], rec["converged"]))
    R = {}
    # exact-scaling check on the integrated loss: loss(2, d3; [0,40]) = loss(1, 2 d3; [0,40]) + c(1, 2 d3) ln 4
    scaling = {}
    for r in rows:
        if r["N"] == 2:
            for s in rows:
                if s["N"] == 1 and abs(s["d3"] - 2 * r["d3"]) < 1e-12:
                    pred = s["loss_0_40"] + s["c_over_xi"] * math.log(4)
                    scaling["d3_%03d" % int(C.mround(1000 * r["d3"]))] = {
                        "loss_0_40_N2": r["loss_0_40"], "loss_0_40_N1_at_2d3": s["loss_0_40"],
                        "c_N1_at_2d3": s["c_over_xi"], "predicted_N2": pred, "ratio": r["loss_0_40"] / pred,
                        "drift_N2_window_5_10": float(r["window_drift"][1]),
                        "two_drift_N1_at_2d3_window_30_40": float(2 * s["window_drift"][4]),
                        "windowfit_ratio_0_1_0": r["loss_rate"] / (4 * s["loss_rate"])}
    R["points"] = None                                  # placeholder keeps the MATLAB field order
    R["scaling"] = scaling
    x, y = [], []
    for r in rows:
        if r["N"] == 1 and r["loss_rate"] > 1e-7 and (fast or r["converged"]):
            x.append(1 / r["Nd3"])
            y.append(math.log(r["loss_rate"]))
    if len(x) >= 3:
        x, y = np.array(x), np.array(y)
        pf = np.polyfit(x, y, 1)
        yhat = np.polyval(pf, x)
        R["fit_diagnostic"] = {"form": "ln(loss_rate) = a + b/(N d3)  [diagnostic only, not a law]", "a": float(pf[1]),
                               "b": float(pf[0]), "n_points": int(x.size),
                               "R2": float(1 - np.sum((y - yhat) ** 2) / np.sum((y - np.mean(y)) ** 2))}
    else:
        R["fit_diagnostic"] = {"form": "ln(loss_rate) = a + b/(N d3)", "n_points": len(x), "note": "fewer than 3 usable points"}
    R["points"] = rows
    ch, T = _census_and_tokens(rows, scaling, "window_drift", "recoil_rel_dev")
    chg, Tg = _census_and_tokens(rows, scaling, "window_drift_grid", "recoil_rel_dev_grid")
    R["channels"] = ch
    R["tokens"] = T
    R["tokens_grid_peak"] = Tg
    R["channels_grid_peak"] = {"pair_closed_until": chg["pair_closed_until"], "pair_open_below_late": chg["pair_open_below_late"]}
    R["tokens_changed_by_peak"] = sorted(k for k in T if T[k] != Tg.get(k))
    sc = scaling["d3_050"]
    closed_until, ep = ch["pair_closed_until"], ch["epochs"]
    R["protocol"] = PROTOCOL
    R["wall_s"] = time.time() - t0
    R["runtime"] = runtime_string()
    R["fast"] = bool(fast)
    R["script_version"] = "0.1.2 (Python port, SOTHE-P2 1.0.0) + sub-grid peak"
    print("\n==== T-10 / T-12 (Python port of run_g4c_loss_scan 0.1.2; RECORDED) ====")
    print("scaling check (2,0.05) vs (1,0.10)+c ln4: measured %.4e predicted %.4e ratio %.4f" % (
        sc["loss_0_40_N2"], sc["predicted_N2"], sc["ratio"]))
    print("pair channel (comoving, TOD included): closed until xi = %d, open below omega = %.3f in [30,40]; "
          "upstream low-k open fraction %.3f" % (int(C.mround(closed_until)), ep[4]["pair_open_below"], ep[4]["Phi_up_low"]))
    for k, v in T.items():
        if k != "G4C_LOSS_TABLE":
            print("  %-26s %s" % (k, v))
    print(T["G4C_LOSS_TABLE"])
    print("peak estimator: %s; tokens whose text differs with the grid-maximum peak (MATLAB 0.1.1): %s" % (
        PEAK_ESTIMATOR, ", ".join("%s %s -> %s" % (k, Tg[k] if k != "G4C_LOSS_TABLE" else "(table)", T[k] if k != "G4C_LOSS_TABLE" else "(table)")
                                  for k in R["tokens_changed_by_peak"]) or "none"))
    fd = R["fit_diagnostic"]
    if "b" in fd:
        print("diagnostic fit: ln(loss) = %.3f %+.4f/(N d3)   (%d points, R^2 = %.4f)" % (fd["a"], fd["b"], fd["n_points"], fd["R2"]))
    print("wall %.0f s, %s" % (R["wall_s"], R["runtime"]))
    rep = {"points": R["points"], "scaling": R["scaling"], "fit_diagnostic": R["fit_diagnostic"], "channels": R["channels"],
           "tokens": R["tokens"], "protocol": R["protocol"], "wall_s": R["wall_s"], "runtime": R["runtime"],
           "fast": R["fast"], "script_version": R["script_version"],
           "peak_estimator": PEAK_ESTIMATOR, "tokens_grid_peak": R["tokens_grid_peak"],
           "channels_grid_peak": R["channels_grid_peak"], "tokens_changed_by_peak": R["tokens_changed_by_peak"]}
    jsonio.write(os.path.join(workdir, "g4c_loss_report.json"), rep, pretty=False)
    # raw data for analysis (not part of the MATLAB report): every propagation and every Newton solution
    from .g4a import write_series_csv
    sd = os.path.join(workdir, "series")
    os.makedirs(sd, exist_ok=True)
    for (kind, N, d3), X in res.items():
        if kind in ("base", "half", "n2"):
            write_series_csv(os.path.join(sd, "%s_N%d_d%.3f.csv" % (kind, N, d3)), X)
    npz = {}
    for (kind, N, d3), X in res.items():
        if kind == "tod" and "phi" in X:
            tag = "N%d_d%.3f" % (N, d3)
            npz["tau"] = X["tau"]
            npz["phi_" + tag] = X["phi"]
            npz["newton_history_" + tag] = X["newton_history"]
    if npz:
        np.savez_compressed(os.path.join(workdir, "newton_solutions.npz"), **npz)
    return rep


def op_series():
    """The operating-point propagation (N = 2, delta3 = 0.05, base protocol) used by Fig. 4."""
    return kit.splitstep_loss(2.0, 0.05, 40.0, 0.002, 8192, 400.0, 0.0)
