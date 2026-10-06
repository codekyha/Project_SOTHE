"""Phase C: per-block reductions, tables (CSV), the oracle ledger (ORACLES.json), summary.json and SUMMARY.md.

Oracle classes:
  hard   identities, reproductions of recorded runs and exact checks: a failure means the pipeline is wrong (exit 2)
  claim  a statement of the manuscript or of REVIEW_G5 put to the test by the sweep: pass/fail is a finding
  info   a number worth tracking, with no pass rule
Every value is RECORDED Python output (HR-3)."""
import math
import os

import numpy as np

from sothe_p2 import phase2, suite

from . import PACK_ROOT, b_census, b_loss, b_raman, b_ref25, grid
from .tasks import _f
from .util import RULE, fnum, read_json, utc, write_csv, write_json

REF = os.path.join(PACK_ROOT, "reference", "job530047")
FULL_LOSS = {"xi_end": 40.0, "dxi": 0.002, "n": 8192, "L": 400.0}


# ============================================================ ledger
class Ledger(object):
    def __init__(self):
        self.rows = []

    def add(self, oid, block, cls, desc, value, ref, ok, tol=None, note=""):
        self.rows.append({"id": oid, "block": block, "class": cls, "description": desc, "value": value, "reference": ref,
                          "tolerance": tol, "pass": (None if ok is None else bool(ok)), "note": note})

    def missing(self, oid, block, cls, desc, why):
        self.rows.append({"id": oid, "block": block, "class": cls, "description": desc, "value": None, "reference": None,
                          "tolerance": None, "pass": False if cls == "hard" else None, "note": "not evaluated: " + why})

    def hard_failures(self):
        return [r for r in self.rows if r["class"] == "hard" and r["pass"] is False]


def _mk(rundir, name):
    d = os.path.join(rundir, name)
    os.makedirs(d, exist_ok=True)
    return d


def _near(a, b, tol=1e-9):
    return abs(float(a) - float(b)) <= tol


def _find(rows, **kv):
    for r in rows:
        if all(_near(r[k], v) if isinstance(v, float) else r[k] == v for k, v in kv.items()):
            return r
    return None


def _same_protocol(P):
    return all(_near(P[k], v) for k, v in FULL_LOSS.items())


# ============================================================ B1, B9
def red_B1(rundir, cfg, res, O, S):
    d = _mk(rundir, "B1_kinematics")
    rows, flows, ident = [], [], None
    for N in cfg["N"]:
        r = res.get("B1_kin_N%s" % _f(N))
        if r is None:
            continue
        rows += r["rows"]
        flows += r["flows"]
        ident = ident or r["identity"]
    S["_B1_rows"] = rows
    if not rows:
        O.missing("O-B1", "B1", "hard", "kinematic layer", "no B1 result")
        return
    cols = ["N", "k_op", "d3", "kg", "S", "c_match", "lambda", "kappa_closed", "kappa_num", "deficit_rel", "eq7pp_rel", "eq7pp_dev",
            "drive_spread_abs", "drive_spread_rel", "band_rel", "band_abs", "eq7pp_drive_max_dev", "mirror_max", "T_H", "HP", "mu"]
    write_csv(os.path.join(d, "kin_grid.csv"), cols, [[r[c] for c in cols] for r in rows])
    fr = []
    for f in flows:
        for g, kn, kc in zip(f["kg"], f["kappa_num"], f["kappa_closed"]):
            fr.append([f["N"], f["k_op"], f["d3"], g, kn, kc, kn / (2 * math.pi)])
    write_csv(os.path.join(d, "kin_flows.csv"), ["N", "k_op", "d3", "kappa_g", "kappa_num", "kappa_closed", "T_H"], fr)
    T = suite.phase2_params({})["T"]
    op = _find(rows, N=cfg["op"]["N"], k_op=cfg["op"]["k_op"], d3=cfg["op"]["d3"], kg=cfg["op"]["kg"])
    if ident is not None:
        O.add("O-B1-identity", "B1", "hard", "kappa_extract(n = 4096) is suite.kinematic_kappa at the operating point (bitwise)",
              ident["here"], ident["suite"], ident["equal"])
    if op is not None:
        O.add("O-B1-op-kappa", "B1", "hard", "kappa_num at (N, k_op, delta3, kappa_g) = (2, 1, 0.05, 0.3) against the v0.1.0 pin",
              op["kappa_num"], T["kappa_kinematic"], _near(op["kappa_num"], T["kappa_kinematic"], 1e-12), 1e-12)
        O.add("O-B1-op-spread", "B1", "hard", "drive spread of kappa over drive in [0.7, 1.3] at the operating point (pin 9.10e-5, 3 figures)",
              op["drive_spread_abs"], T["kappa_drive_spread"], _near(op["drive_spread_abs"], T["kappa_drive_spread"], 5e-7), 5e-7)
    mm = max(r["mirror_max"] for r in rows)
    O.add("O-B1-mirror", "B1", "hard", "kappa(drive) = kappa(2 - drive) over the whole grid (symmetric nodes)", mm, 0.0, mm <= 1e-10, 1e-10)
    mono = all(f["monotone"] for f in flows)
    O.add("O-B1-flow", "B1", "hard", "kappa decreases strictly with kappa_g on every flow (N, k_op, delta3)", mono, True, mono)
    rel = max(abs(r["eq7pp_dev"]) / abs(r["deficit_rel"]) for r in rows if r["deficit_rel"] != 0)
    O.add("O-B1-eq7pp", "B1", "hard", "Eq. (7''): max |extracted deficit - Eq. (7'')| / |deficit| over the grid (drive = 1)",
          rel, 0.0, rel <= 0.02, 0.02)
    worst = max(r["drive_spread_rel"] / r["band_rel"] for r in rows)
    O.add("O-B1-band", "B1", "claim", "drive spread of kappa within the Eq. (7'') band lambda^2 dtau^2 / 4 at every grid point",
          worst, 1.0, worst <= 1.0 + 1e-6, note="ratio spread / band, maximum over the grid")
    S["B1"] = {"n_rows": len(rows), "op": op, "identity": ident, "mirror_max": mm, "eq7pp_rel_max": rel, "spread_over_band_max": worst,
               "kappa_range": [min(r["kappa_num"] for r in rows), max(r["kappa_num"] for r in rows)]}


def red_B9(rundir, cfg, res, O, S):
    d = _mk(rundir, "B9_parity")
    r = res.get("B9_parity")
    if r is None:
        O.missing("O-B9", "B9", "hard", "parity tables", "no B9 result")
        return
    t2, t3 = r["table2"], r["table3"]
    write_csv(os.path.join(d, "table2_parity.csv"), ["parity", "n_tau", "deficit_over_lamdt2", "eq7pp_over_lamdt2", "s_over_dt", "published",
                                                     "dev_published"], [[x[k] for k in ("parity", "n_tau", "deficit_over_lamdt2",
                                                                                         "eq7pp_over_lamdt2", "s_over_dt", "published", "dev_published")] for x in t2])
    write_csv(os.path.join(d, "table3_continuum.csv"), ["N", "extracted", "eq7pp", "relative", "published_extracted", "published_eq7pp"],
              [[x[k] for k in ("N", "extracted", "eq7pp", "relative", "published_extracted", "published_eq7pp")] for x in t3])
    devs = [abs(x["dev_published"]) for x in t2 if x["dev_published"] is not None]
    m2 = max(devs) if devs else None
    O.add("O-B9-table2", "B9", "hard", "Supplement Table 2 (deficit / (lambda dtau)^2, odd and even n_tau) against the published 6 decimals",
          m2, 0.0, m2 is not None and m2 <= 5e-7, 5e-7)
    worst = 0.0
    for x in t3:
        for key, pk in (("extracted", "published_extracted"), ("eq7pp", "published_eq7pp")):
            if x[pk]:
                worst = max(worst, abs(x[key] / x[pk] - 1))
    O.add("O-B9-table3", "B9", "hard", "Supplement Table 3 (N = 1, 2, 3 at kappa_g = 0) against the published 5 figures", worst, 0.0,
          worst <= 1e-4, 1e-4, note="relative deviation")
    S["B9"] = {"table2_max_dev": m2, "table3_max_rel": worst, "table2": t2, "table3": t3}


# ============================================================ B3
def loss_reductions(cfg, res):
    P = cfg["loss"]
    out = {}
    for N, d3 in P["points"]:
        tid = "B3_loss_N%s_d%s_%%s" % (_f(N), _f(d3))
        base = res.get(tid % "base")
        if base is None:
            continue
        out["%s|%s" % (_f(N), _f(d3))] = b_loss.reduce_point(base, res.get(tid % "half"), res.get(tid % "n2"))
    return out


def red_B3(rundir, cfg, res, O, S):
    d = _mk(rundir, "B3_loss")
    red = loss_reductions(cfg, res)
    S["_B3_red"] = red
    if not red:
        O.missing("O-B3", "B3", "hard", "loss and drift", "no B3 result")
        return
    pts = sorted(red.values(), key=lambda r: (r["d3"], r["N"]))
    cols = ["N", "d3", "Nd3", "loss_rate", "loss_0_40", "loss_first_hp", "c_over_xi", "kappa_hp", "hawking_period", "k_s_early", "k_s_late",
            "recoil_rel_dev", "drift_leading_order", "rel_dev", "rel_dev_0_40", "drift_conv_max_abs", "converged"]
    wc = ["window_drift", "window_rate", "window_ks_rate"]
    head = cols + ["%s_%d" % (c, i + 1) for c in wc for i in range(5)]
    write_csv(os.path.join(d, "loss_points.csv"), head, [[r[c] for c in cols] + [r[c][i] for c in wc for i in range(5)] for r in pts])
    sd = _mk(rundir, os.path.join("B3_loss", "series"))
    for N, d3 in cfg["loss"]["points"]:
        b = res.get("B3_loss_N%s_d%s_base" % (_f(N), _f(d3)))
        if b:
            write_csv(os.path.join(sd, "loss_N%s_d%s.csv" % (_f(N), _f(d3))), ["xi", "Q", "tau_c", "k_s", "tau_c_grid"],
                      list(zip(b["xi"], b["Q"], b["tau_c"], b["k_s"], b["tau_c_grid"])))
    ref = read_json(os.path.join(REF, "g4c_loss_report.json"))
    cmp = []
    if ref and _same_protocol(cfg["loss"]):
        for rp in ref["points"]:
            k = "%s|%s" % (_f(float(rp["N"])), _f(float(rp["d3"])))
            if k not in red:
                continue
            r = red[k]
            dv = {"N": rp["N"], "d3": rp["d3"],
                  "loss_0_40_rel": abs(r["loss_0_40"] / rp["loss_0_40"] - 1), "loss_rate_rel": abs(r["loss_rate"] / rp["loss_rate"] - 1),
                  "loss_first_hp_rel": abs(r["loss_first_hp"] / rp["loss_first_hp"] - 1),
                  "window_drift_abs": max(abs(a - b) for a, b in zip(r["window_drift"], rp["window_drift"])),
                  "window_rate_rel": max(abs(a / b - 1) for a, b in zip(r["window_rate"], rp["window_rate"]) if b != 0)}
            if cfg["loss"]["conv"] and rp.get("rel_dev") is not None:
                dv["rel_dev_abs"] = abs(r["rel_dev"] - rp["rel_dev"])
            cmp.append(dv)
        worst_rel = max(max(c["loss_0_40_rel"], c["loss_rate_rel"], c["loss_first_hp_rel"]) for c in cmp) if cmp else None
        worst_dr = max(c["window_drift_abs"] for c in cmp) if cmp else None
        ok = bool(cmp) and len(cmp) == 9 and worst_rel <= 1e-6 and worst_dr <= 1e-7
        O.add("O-B3-530047", "B3", "hard", "the 9 recorded points of job 530047 (S2): loss over [0, 40], loss rate, loss over one Hawking "
              "period, window drifts", {"points": len(cmp), "loss_rel_max": worst_rel, "drift_abs_max": worst_dr}, "g4c_loss_report.json",
              ok, "rel 1e-6 / abs 1e-7")
    else:
        O.add("O-B3-530047", "B3", "info", "comparison with job 530047", None, None, None, note="protocol differs (fast run): not compared")
    write_json(os.path.join(d, "loss_report.json"), {"points": pts, "reference_comparison": cmp, "protocol": cfg["loss"], "status": "RECORDED"})
    nconv = [r for r in pts if r["converged"] is not None]
    if nconv:
        O.add("O-B3-conv", "B3", "info", "points whose loss rate changes by <= 10% under dxi/2 and 2n", "%d of %d" % (sum(r["converged"] for r in nconv),
              len(nconv)), None, None, note="non-converged points carry an exponentially small loss rate (noise-dominated fit)")
    a = red.get("1|0.05")
    if a:
        O.add("O-B3-anchor-loss", "B3", "claim", "REVIEW_G5 S1-1 A: loss at (N, delta3) = (1, 0.05) over xi <= 40 is 9.3e-4, over one "
              "Hawking period after xi = 10 is 1.2e-4 (k_op = 1 period)", [a["loss_0_40"], a["loss_first_hp"]], [9.3e-4, 1.2e-4],
              abs(a["loss_0_40"] - 9.3e-4) <= 0.05e-4 and abs(a["loss_first_hp"] - 1.2e-4) <= 0.05e-4, "half a unit of the last digit")
    o = red.get("2|0.05")
    if o:
        O.add("O-B3-op-drift", "B3", "claim", "REVIEW_G5 B.6: drift at the operating point, windows [0,5] and [30,40]: -0.3133, -0.5256",
              [o["window_drift"][0], o["window_drift"][4]], [-0.3133, -0.5256],
              abs(o["window_drift"][0] + 0.3133) <= 5e-5 and abs(o["window_drift"][4] + 0.5256) <= 5e-5, 5e-5)
    S["B3"] = {"n_points": len(pts), "reference": cmp, "op": o, "anchor": a,
               "loss_0_40_range": [min(r["loss_0_40"] for r in pts), max(r["loss_0_40"] for r in pts)]}


# ============================================================ B4, B2
def _census_rows(cfg, res):
    rows = []
    for N in cfg["N"]:
        for d3 in cfg["d3_kin"]:
            r = res.get("B4_census_N%s_d%s" % (_f(N), _f(d3)))
            if r:
                rows += r["rows"]
    return rows


def red_B4(rundir, cfg, res, O, S):
    d = _mk(rundir, "B4_census")
    rows = _census_rows(cfg, res)
    S["_B4_rows"] = rows
    if not rows:
        O.missing("O-B4", "B4", "hard", "channel census", "no B4 result")
        return
    cols = ["N", "d3", "k_op", "kg", "frame", "Vd", "S", "c_match", "two_c", "omega_max_pair_used", "omega_max_pair", "omega_max_pair_no_tod",
            "Phi_pair", "Phi_up_low", "pair_open_below", "k_star", "mu", "kappa", "window_lo", "window_hi", "window_frac", "piN_over_S",
            "n_at_mu", "EN_at_mu"]
    write_csv(os.path.join(d, "census.csv"), cols + ["epoch_lo", "epoch_hi"],
              [[r.get(c) for c in cols] + ((r["epoch_window"] or [None, None]) if r.get("epoch_window") else [None, None]) for r in rows])
    op = cfg["op"]
    opr = [r for r in rows if _near(r["N"], op["N"]) and _near(r["d3"], op["d3"]) and _near(r["k_op"], op["k_op"]) and _near(r["kg"], op["kg"])]
    ep = [r for r in opr if r["frame"] != "lab"]
    closed_until, late_edge = 40, None
    for e in ep:
        if e["Phi_pair"] > 0 and closed_until == 40:
            closed_until = e["epoch_window"][0]
    if ep:
        late_edge = ep[-1]["pair_open_below"]
    ref = read_json(os.path.join(REF, "g4c_loss_report.json"))
    full_census = cfg["census_nw"] == 311
    if ref and ep and _same_protocol(cfg["loss"]) and full_census:
        rc = ref["channels"]
        ok = (closed_until == rc["pair_closed_until"]) and _near(late_edge, rc["pair_open_below_late"], 1e-12)
        O.add("O-B4-op-epochs", "B4", "hard", "comoving census at the operating point from the B3 drifts: pair channel closed until xi = 20, "
              "open below omega = 0.205 in [30, 40] (job 530047)", [closed_until, late_edge], [rc["pair_closed_until"], rc["pair_open_below_late"]], ok)
    elif ep:
        O.add("O-B4-op-epochs", "B4", "info", "comoving census at the operating point", [closed_until, late_edge], None, None, note="fast run")
    if ref:
        worst, n = 0.0, 0
        for key, g in ref["channels"]["grid"].items():
            lab = g["lab"]
            r = _find([x for x in rows if x["frame"] == "lab"], N=float(lab["N"]), k_op=float(lab["k_op"]), d3=float(lab["d3"]), kg=float(lab["kg"]))
            if r is None:
                continue
            n += 1
            for k in ("omega_max_pair", "omega_max_pair_no_tod", "Phi_pair", "Phi_up_low", "pair_open_below"):
                worst = max(worst, abs(r[k] - lab[k]))
        if n:
            O.add("O-B4-lab-grid", "B4", "hard" if full_census else "info", "laboratory-frame census at the %d (N, k_op) points of job 530047 channels.grid (same function)" % n,
                  worst, 0.0, (worst <= 1e-12) if full_census else None, 1e-12)
    # Hawking-window table (REVIEW_G5 S1-1) at delta3 = 0.05, kappa_g = 0.3: lab frame and the last drift epoch
    tab = []
    for N in cfg["N"]:
        for kop in cfg["k_op"]:
            sel = [r for r in rows if _near(r["N"], N) and _near(r["k_op"], kop) and _near(r["d3"], 0.05) and _near(r["kg"], 0.3)]
            if not sel:
                continue
            lab = [r for r in sel if r["frame"] == "lab"][0]
            late = [r for r in sel if r["frame"] != "lab"]
            late = late[-1] if late else None
            tab.append({"N": N, "k_op": kop, "mu": lab["mu"], "kappa": lab["kappa"], "piN_over_S": lab["piN_over_S"], "n_at_mu": lab["n_at_mu"],
                        "EN_at_mu": lab["EN_at_mu"], "lab_edge": lab["omega_max_pair_used"], "lab_window": [lab["window_lo"], lab["window_hi"]],
                        "lab_frac": lab["window_frac"], "late_Vd": late["Vd"] if late else None,
                        "late_edge": late["omega_max_pair_used"] if late else None,
                        "late_window": [late["window_lo"], late["window_hi"]] if late else None, "late_frac": late["window_frac"] if late else None,
                        "late_Phi_pair": late["Phi_pair"] if late else None})
    write_csv(os.path.join(d, "hawking_window_d3_005.csv"),
              ["N", "k_op", "mu", "kappa", "piN_over_S", "n_at_mu", "EN_at_mu", "lab_edge", "lab_lo", "lab_hi", "lab_frac", "late_Vd",
               "late_edge", "late_lo", "late_hi", "late_frac", "late_Phi_pair"],
              [[t["N"], t["k_op"], t["mu"], t["kappa"], t["piN_over_S"], t["n_at_mu"], t["EN_at_mu"], t["lab_edge"], t["lab_window"][0],
                t["lab_window"][1], t["lab_frac"], t["late_Vd"], t["late_edge"], (t["late_window"] or [None, None])[0],
                (t["late_window"] or [None, None])[1], t["late_frac"], t["late_Phi_pair"]] for t in tab])
    review = {(0.6, 1.0): 0.748, (1.0, 1.0): 0.4975, (1.4, 1.0): 0.158, (2.0, 1.0): -0.152,
              (0.6, 1.5): 1.965, (1.0, 1.5): 1.760, (1.4, 1.5): 1.511, (2.0, 1.5): 1.611}
    devs = []
    for (N, kop), val in review.items():
        t = [x for x in tab if _near(x["N"], N) and _near(x["k_op"], kop)]
        if t and t[0]["late_edge"] is not None:
            devs.append({"N": N, "k_op": kop, "review": val, "sweep": t[0]["late_edge"], "dev": t[0]["late_edge"] - val})
    if devs:
        O.add("O-B4-review-table", "B4", "info", "REVIEW_G5 S1-1 table: comoving edge omega_max in [30, 40] (review: drift of N = 1 scaled; "
              "sweep: drift measured at each N)", max(abs(x["dev"]) for x in devs), 0.0, None, note="max |sweep - review|; rows in summary.json")
    wmax = max((x["N"] for x in tab if x["lab_window"][0] is not None or (x["late_window"] and x["late_window"][0] is not None)), default=None)
    O.add("O-B4-Nmax", "B4", "claim", "REVIEW_G5 S1-1: escape inside the band requires N < sqrt(2 omega_max) = 1.789 (largest grid N with a window)",
          wmax, math.sqrt(3.2), wmax is None or wmax < math.sqrt(3.2))
    S["B4"] = {"n_rows": len(rows), "op_pair_closed_until": closed_until, "op_pair_open_below_late": late_edge, "window_table": tab,
               "review_table_devs": devs, "largest_N_with_window": wmax}


def red_B2(rundir, cfg, res, O, S):
    d = _mk(rundir, "B2_entanglement")
    rows = S.get("_B1_rows") or []
    cen = S.get("_B4_rows") or []
    if not rows:
        O.missing("O-B2", "B2", "hard", "entanglement layer", "no B1 rows")
        return
    om = cfg["omega"]
    out = []
    for r in rows:
        k = r["kappa_num"]
        L = b_census.entanglement_layer(k, om)
        sel = [c for c in cen if _near(c["N"], r["N"]) and _near(c["k_op"], r["k_op"]) and _near(c["d3"], r["d3"]) and _near(c["kg"], r["kg"])]
        lab = [c for c in sel if c["frame"] == "lab"]
        late = [c for c in sel if c["frame"] != "lab"]
        lab = lab[0] if lab else None
        late = late[-1] if late else None

        def en(w):
            return b_census.EN_thermal(w, k) if (w is not None and k > 0) else None

        mu = 0.5 * r["N"] ** 2
        out.append([r["N"], r["k_op"], r["d3"], r["kg"], k, k / (2 * math.pi), L["EN_max"], L["nu_max"], L["IR_asymptote"],
                    en(mu) if mu <= om[-1] else None,
                    lab["window_lo"] if lab else None, lab["window_hi"] if lab else None, lab["window_frac"] if lab else None,
                    en(lab["window_lo"]) if lab else None, en(lab["window_hi"]) if lab else None,
                    late["window_lo"] if late else None, late["window_hi"] if late else None, late["window_frac"] if late else None,
                    en(late["window_lo"]) if late else None, en(late["window_hi"]) if late else None])
    write_csv(os.path.join(d, "entanglement_grid.csv"),
              ["N", "k_op", "d3", "kg", "kappa", "T_H", "EN_max_band", "nu_max_band", "IR_asymptote", "EN_at_mu", "lab_lo", "lab_hi", "lab_frac",
               "EN_lab_lo", "EN_lab_hi", "late_lo", "late_hi", "late_frac", "EN_late_lo", "EN_late_hi"], out)
    T = suite.phase2_params({})["T"]
    op = _find(rows, N=cfg["op"]["N"], k_op=cfg["op"]["k_op"], d3=cfg["op"]["d3"], kg=cfg["op"]["kg"])
    curves = {}
    for name, pt in (("op", cfg["op"]), ("anchor", cfg["anchor"])):
        r = _find(rows, N=pt["N"], k_op=pt["k_op"], d3=pt["d3"], kg=pt["kg"])
        if r is None:
            continue
        L = b_census.entanglement_layer(r["kappa_num"], om)
        curves[name] = L
        write_csv(os.path.join(d, "EN_curves_%s.csv" % name), ["omega", "E_N", "nu_minus", "nbar_max", "E_N_eta0.9", "E_N_eta0.5", "E_N_eta0.1"],
                  list(zip(L["omega"], L["EN"], L["nu"], L["nbar_max"], L["EN_loss_eta0.9"], L["EN_loss_eta0.5"], L["EN_loss_eta0.1"])))
    if op is not None and "op" in curves:
        L = curves["op"]
        O.add("O-B2-op-ENmax", "B2", "hard", "E_N at omega = 0.05 at the operating point against the v0.1.0 pin", L["EN_max"], T["E_N_max_nats"],
              _near(L["EN_max"], T["E_N_max_nats"], 1e-9), 1e-9)
        O.add("O-B2-EN3k", "B2", "hard", "E_N at omega = 3 kappa (kappa-independent)", L["EN_at_3kappa"], T["E_N_at_3kappa"],
              abs(L["EN_at_3kappa"] / T["E_N_at_3kappa"] - 1) <= 1e-9, 1e-9)
    S["B2"] = {"n_rows": len(out), "curves": {k: {"EN_max": v["EN_max"], "nu_max": v["nu_max"]} for k, v in curves.items()}}


# ============================================================ B5
def red_B5(rundir, cfg, res, O, S):
    d = _mk(rundir, "B5_bdg")
    sd = _mk(rundir, os.path.join("B5_bdg", "spectra"))
    rows = []
    for N, d3 in cfg["bdg"]["points"]:
        r = res.get("B5_bdg_N%s_d%s" % (_f(N), _f(d3)))
        if r is None:
            continue
        s = r["summary"]
        rows.append(s)
        re, im, kr = np.array(r["spectrum"]["re"]), np.array(r["spectrum"]["im"]), np.array(r["spectrum"]["krein"])
        win = np.abs(re) <= 3 * s["mu"]
        write_csv(os.path.join(sd, "bdg_N%s_d%s.csv" % (_f(N), _f(d3))), ["Re_omega", "Im_omega", "krein"], list(zip(re[win], im[win], kr[win])))
    S["_B5_rows"] = rows
    if not rows:
        O.missing("O-B5", "B5", "hard", "BdG", "no B5 result")
        return
    cols = ["N", "d3", "n_tau", "tau_max", "mu", "newton_converged", "newton_res", "newton_iters", "kmean", "tail_rel", "annulus_maxIm",
            "gap_over_mu", "n_annulus", "krein_pos", "krein_neg", "zero_sector_max", "rho", "sqrt_eps_rho", "dim_ker_M", "dim_ker_M2",
            "sv_M_gap", "sv_M2_gap", "psi0_annulus_maxIm", "psi0_annulus_over_sqrt_d3", "psi0_quartet_maxIm", "psi0_quartet_over_sqrt_d3"]
    write_csv(os.path.join(d, "bdg_points.csv"), cols, [[r.get(c) for c in cols] for r in rows])
    full = cfg["bdg"]["n_tau"] == 500 and _near(cfg["bdg"]["tau_max"], 16.0)
    for d3v in (0.02, 0.05):
        r = _find(rows, N=1.0, d3=d3v)
        if r is None:
            continue
        cls = "hard" if full else "info"
        O.add("O-B5-T02-%03d" % round(d3v * 1000), "B5", cls, "G4a T-02 at (N, delta3) = (1, %g): Newton converged, residual <= 1e-12" % d3v,
              r["newton_res"], 1e-12, r["newton_converged"] and r["newton_res"] <= 1e-12)
        O.add("O-B5-T03-%03d" % round(d3v * 1000), "B5", cls, "G4a T-03 at (1, %g): annulus max|Im| <= 1e-10, zero sector <= 1e-6, "
              "dim ker M = 2, dim ker M^2 = 4" % d3v, [r["annulus_maxIm"], r["zero_sector_max"], r["dim_ker_M"], r["dim_ker_M2"]],
              [1e-10, 1e-6, 2, 4], r["annulus_maxIm"] <= 1e-10 and r["zero_sector_max"] <= 1e-6 and r["dim_ker_M"] == 2 and r["dim_ker_M2"] == 4)
    r = _find(rows, N=1.0, d3=0.02)
    if r is not None and full:
        T = suite.phase2_params({})["T"]
        g4a = read_json(os.path.join(REF, "g4a_report.json")) or {}
        rec = g4a.get("T01_maxIm_frozen")
        O.add("O-B5-psi0-020", "B5", "hard", "psi0-linearization rate at (1, 0.02) (the artifact pinned by v0.1.0, retired in v0.2.0) against "
              "the pin and job 530047 (G4a T-01)", r["psi0_annulus_maxIm"], [T["bdg_maxIm_d3_002"], rec],
              _near(r["psi0_annulus_maxIm"], T["bdg_maxIm_d3_002"], 1e-9) and (rec is None or _near(r["psi0_annulus_maxIm"], rec, 1e-9)), 1e-9)
    r = _find(rows, N=1.0, d3=0.05)
    if r is not None and full:
        O.add("O-B5-tail-005", "B5", "hard", "tail of the TOD soliton at (1, 0.05), |phi| / N at |tau| > 0.8 tau_max (G4C_TAIL_005 3.6e-5)",
              r["tail_rel"], 3.5981933823064204e-05, abs(r["tail_rel"] / 3.5981933823064204e-05 - 1) <= 1e-6, 1e-6)
    inside = [r for r in rows if r["d3"] > 0 and r["N"] * r["d3"] <= 0.05 + 1e-12]
    outside = [r for r in rows if r["N"] * r["d3"] > 0.05 + 1e-12]
    if inside:
        w = max(inside, key=lambda r: r["annulus_maxIm"])
        O.add("O-B5-stable-scope", "B5", "claim", "manuscript scope N delta3 <= 0.05: the stationary TOD soliton is spectrally stable at every "
              "grid point (annulus max|Im| <= 1e-10, Newton converged)", {"max_annulus_maxIm": w["annulus_maxIm"], "at": [w["N"], w["d3"]],
                                                                          "all_converged": all(r["newton_converged"] for r in inside)}, 1e-10,
              all(r["newton_converged"] and r["annulus_maxIm"] <= 1e-10 for r in inside))
    if outside:
        O.add("O-B5-beyond-scope", "B5", "info", "N delta3 > 0.05 (outside the manuscript's BdG scope): annulus max|Im| and soliton tail",
              [[r["N"], r["d3"], r["annulus_maxIm"], r["tail_rel"]] for r in outside], None, None,
              note="[N, delta3, max|Im|, tail]; a large tail means the periodic stationary state carries radiation")
    S["B5"] = {"n_points": len(rows), "rows": [{k: r.get(k) for k in ("N", "d3", "newton_converged", "newton_res", "tail_rel", "annulus_maxIm",
                                                                      "zero_sector_max", "dim_ker_M", "dim_ker_M2", "psi0_annulus_maxIm",
                                                                      "psi0_annulus_over_sqrt_d3", "gap_over_mu")} for r in rows]}


# ============================================================ B6
def red_B6(rundir, cfg, res, O, S):
    d = _mk(rundir, "B6_raman")
    sd = _mk(rundir, os.path.join("B6_raman", "series"))
    R = cfg["raman"]
    runs, seen = [], set()
    for N, d3 in R["points"]:
        for model in R["models"]:
            for T0 in R["T0_fs"]:
                tid = "B6_raman_N%s_d%s_%s_T%s" % (_f(N), _f(d3), model, _f(T0) if model != "none" else "any")
                if tid in seen or tid not in res:
                    continue
                seen.add(tid)
                runs.append(res[tid])
    S["_B6_runs"] = runs
    gord = [res[k] for k in sorted(res) if k.startswith("B6_gordon_")]
    if not runs and not gord:
        O.missing("O-B6", "B6", "hard", "Raman extension", "no B6 result")
        return
    nw = max((len(r["windows"]) for r in runs), default=0)
    cols = ["N", "d3", "T0_fs", "model", "f_R", "s", "tauR_eff", "gordon_rate", "initial_ks_rate", "initial_kfull_rate", "early_ks_rate",
            "early_kfull_rate", "shift_total_ks", "loss_end", "k_s_end", "k_full_end", "tau_c_end", "n_roll"]
    head = cols + ["%s_%d" % (c, i + 1) for c in ("window_drift", "window_ks_rate", "window_kfull_rate") for i in range(nw)]
    rows = []
    for r in runs:
        row = [r[c] if c != "T0_fs" or r["model"] != "none" else None for c in cols]
        for c in ("window_drift", "window_ks_rate", "window_kfull_rate"):
            row += [r[c][i] if i < len(r[c]) else None for i in range(nw)]
        rows.append(row)
        s = r["series"]
        write_csv(os.path.join(sd, "raman_N%s_d%s_%s_T%s.csv" % (_f(r["N"]), _f(r["d3"]), r["model"], _f(r["T0_fs"]) if r["model"] != "none" else "any")),
                  ["xi", "Q", "tau_c", "k_s", "k_full", "P"], list(zip(s["xi"], s["Q"], s["tau_c"], s["k_s"], s["k_full"], s["P"])))
    write_csv(os.path.join(d, "raman_runs.csv"), head, rows)
    write_csv(os.path.join(d, "raman_gordon.csv"), ["label", "tauR_eff", "rate_full", "rate_windowed", "gordon", "ratio", "photon_drift"],
              [[g["label"], g["tauR_eff"], g["rate_full"], g["rate_windowed"], g["gordon"], g["ratio"], g["photon_drift"]] for g in gord])
    if gord:
        worst = max(abs(g["ratio"] - 1) for g in gord)
        O.add("O-B6-gordon", "B6", "hard", "soliton self-frequency shift at N = 1, delta3 = 0: full-field centroid rate over xi in [1, 5] against "
              "Gordon's -(8/15) tau_R (linear model and Blow-Wood at T0 = 1000 fs)", [[g["label"], g["ratio"]] for g in gord], 1.0,
              worst <= 0.05, 0.05, note="ratio rate / Gordon")
    # f_R = s = 0 against the split step of B3 (same point, same protocol)
    red3 = S.get("_B3_red") or {}
    for r in [x for x in runs if x["model"] == "none"]:
        b = red3.get("%s|%s" % (_f(r["N"]), _f(r["d3"])))
        if b is None or not _same_protocol(dict(xi_end=r["xi_end"], dxi=r["dxi"], n=r["n"], L=r["L"])) or not _same_protocol(cfg["loss"]):
            continue
        dl = abs(r["loss_end"] / b["loss_0_40"] - 1)
        dd = max(abs(a - c) for a, c in zip(r.get("window_drift_05") or r["window_drift"], b["window_drift"]))
        O.add("O-B6-none-N%s" % _f(r["N"]), "B6", "hard", "RK4IP with f_R = s = 0 against the split step of B3 at (N, delta3) = (%g, %g): "
              "loss over [0, 40] and window drifts" % (r["N"], r["d3"]), [dl, dd], [0, 0], dl <= 1e-3 and dd <= 1e-3, "rel 1e-3 / abs 1e-3")
    # the S1-3 table, simulated
    tab = []
    for N, d3 in R["points"]:
        base = [x for x in runs if x["model"] == "none" and _near(x["N"], N) and _near(x["d3"], d3)]
        recoil = None
        if base:
            b = base[0]
            mids = np.array([0.5 * (lo + hi) for lo, hi in b["windows"]])
            wd = np.array(b.get("window_drift_05") or b["window_drift"])
            recoil = float(np.polyfit(mids, wd, 1)[0]) if len(mids) >= 2 else None
        for model in ("lin3", "bw", "bw_ss", "ss"):
            for T0 in R["T0_fs"]:
                x = [y for y in runs if y["model"] == model and _near(y["N"], N) and _near(y["d3"], d3) and _near(y["T0_fs"], T0)]
                if not x:
                    continue
                x = x[0]
                kop = 1.0 if _near(N, 2.0) else 1.5
                kap = b_census.kappa_closed(N, kop, d3, 0.3)
                HP = 2 * math.pi / kap
                xs, ks = np.asarray(x["series"]["xi"]), np.asarray(x["series"]["k_s"])
                dk_hp = float(np.interp(HP, xs, ks) - ks[0]) if xs[-1] >= HP else None
                tab.append({"N": N, "d3": d3, "model": model, "T0_fs": T0, "tauR_eff": x["tauR_eff"], "gordon_rate": x["gordon_rate"],
                            "initial_rate_ks": x["initial_ks_rate"], "initial_rate_kfull": x["initial_kfull_rate"],
                            "early_rate_ks": x["early_ks_rate"], "dk_first_HP": dk_hp,
                            "HP": HP, "k_op_for_HP": kop, "shift_total_ks": x["shift_total_ks"], "linear_extrapolation": x["gordon_rate"] * x["xi_end"],
                            "recoil_dVd_dxi": recoil, "times_recoil": (abs(x["initial_kfull_rate"]) / abs(recoil)) if recoil else None,
                            "late_drift": x["window_drift"][-1] if x["window_drift"] else None, "loss_end": x["loss_end"]})
    keys = ["N", "d3", "model", "T0_fs", "tauR_eff", "gordon_rate", "initial_rate_kfull", "initial_rate_ks", "early_rate_ks", "dk_first_HP", "HP", "k_op_for_HP",
            "shift_total_ks", "linear_extrapolation", "recoil_dVd_dxi", "times_recoil", "late_drift", "loss_end"]
    write_csv(os.path.join(d, "raman_S1-3_table.csv"), keys, [[t[k] for k in keys] for t in tab])
    c = [t for t in tab if t["model"] == "lin3" and _near(t["N"], 2.0) and _near(t["T0_fs"], 50.0)]
    if c:
        c = c[0]
        O.add("O-B6-S1-3-rate", "B6", "claim", "REVIEW_G5 S1-3: at N = 2, T0 = 50 fs (T_R = 3 fs) the self-frequency shift starts at Gordon's "
              "0.512 per xi, about 78 times the recoil-driven drift change", [abs(c["initial_rate_kfull"]), abs(c["initial_rate_ks"]),
                                                                               c["times_recoil"]], [0.512, 0.512, 78],
              abs(abs(c["initial_rate_kfull"]) / 0.512 - 1) <= 0.10 and c["times_recoil"] is not None and c["times_recoil"] > 30,
              note="[full-field centroid, core centroid, ratio] over xi in [0, 0.5]; 10% band on the full-field rate (the core window "
                   "biases the rate low by about 4%, see O-B6-gordon)")
        if c["shift_total_ks"] is not None and _near(cfg["raman"]["xi_end"], 40.0):
            O.add("O-B6-S1-3-total", "B6", "claim", "REVIEW_G5 S1-3: the shift accumulates linearly, about 20.5 over xi = 40 (and 1.93 per "
                  "Hawking period)", [c["shift_total_ks"], c["dk_first_HP"]], [c["linear_extrapolation"], c["gordon_rate"] * c["HP"]],
                  abs(c["shift_total_ks"] / c["linear_extrapolation"] - 1) <= 0.2,
                  note="simulated k_s(40) - k_s(0) and k_s(HP) - k_s(0); the rate falls as the red-shifted soliton meets stronger dispersion")
    # comoving census along the Raman runs (k_op = 1 and 1.5)
    cen = []
    for r in runs:
        for kop in (1.0, 1.5):
            for (lo, hi), vd in zip(r["windows"], r["window_drift"]):
                c = b_census.census_one(r["N"], kop, r["d3"], vd, 0.3, 311)
                cen.append([r["N"], r["d3"], r["model"], r["T0_fs"] if r["model"] != "none" else None, kop, lo, hi, vd, c["omega_max_pair_used"],
                            c["Phi_pair"], c["Phi_up_low"], c["window_lo"], c["window_hi"], c["window_frac"]])
    write_csv(os.path.join(d, "raman_census.csv"), ["N", "d3", "model", "T0_fs", "k_op", "xi_lo", "xi_hi", "Vd", "omega_max_pair", "Phi_pair",
                                                   "Phi_up_low", "window_lo", "window_hi", "window_frac"], cen)
    S["B6"] = {"n_runs": len(runs), "gordon": gord, "S1_3_table": tab, "T_R_BW_fs": b_raman.F_R * b_raman.M1_FS, "T_R_lin_fs": b_raman.T_R_LIN_FS}


# ============================================================ B7
def red_B7(rundir, cfg, res, O, S):
    d = _mk(rundir, "B7_scatter")
    Sc = cfg["scatter"]
    pts, orow = [], []
    nerr = 0
    for N, kop in Sc["points"]:
        r = res.get("B7_scatter_N%s_k%s" % (_f(N), _f(kop)))
        if r is None:
            continue
        pts.append(r["summary"])
        nerr += r["summary"]["n_err"]
        for x in r["omega_records"]:
            orow.append([N, kop, x["omega"], " ".join(x["ins"]), " ".join(x["outs"]), x["two_sided"], x["r_conf"], x["n_H"], x["planck"],
                         x["boltzmann"], x["n_H_over_planck"], x["r_conf_over_boltzmann"], x["fluxbal_max"], x["resid_max"]])
    S["_B7_pts"] = pts
    if orow:
        write_csv(os.path.join(d, "scatter_omega.csv"), ["N", "k_op", "omega", "ins", "outs", "two_sided", "r_conf", "n_H", "planck", "boltzmann",
                                                         "n_H_over_planck", "r_conf_over_boltzmann", "fluxbal", "resid"], orow)
    if pts:
        cols = ["N", "k_op", "kg", "kappa0", "c", "S", "mu", "L", "h", "n_omega", "n_ok", "n_err", "frac_two_sided", "fluxbal_max", "resid_max",
                "r_conf_min", "r_conf_max", "r_conf_spread_rel", "n_H_over_planck_min", "n_H_over_planck_max"]
        write_csv(os.path.join(d, "scatter_points.csv"), cols + ["kappa_fit_n_H"],
                  [[p[c] for c in cols] + [(p["thermal_fit_n_H"] or {}).get("kappa_fit")] for p in pts])
        fb = max(p["fluxbal_max"] for p in pts)
        O.add("O-B7-flux", "B7", "hard", "Krein-signed flux balance sum_out s|S|^2 - s_in over every solve of the grid", fb, 0.0, fb <= 1e-8, 1e-8)
        O.add("O-B7-errors", "B7", "hard", "solves that failed (non-square exterior-mode system)", nerr, 0, nerr == 0)
    kn = res.get("B7_kaup_null")
    if kn:
        write_json(os.path.join(d, "kaup_null.json"), kn)
        R2 = max(x["R2"] for x in kn["kaup"])
        T2 = max(abs(x["T2"] - 1) for x in kn["kaup"])
        O.add("O-B7-kaup", "B7", "hard", "N = 1, V = 0 (static NLS soliton, reflectionless): max |R|^2 and max ||T|^2 - 1|", [R2, T2], [0, 0],
              R2 <= 1e-6 and T2 <= 1e-6, 1e-6)
        mix = max(x["uv_mix_max"] for x in kn["null"])
        O.add("O-B7-null", "B7", "hard", "flow without soliton: no u <-> v conversion (exact zero)", mix, 0.0, mix == 0.0)
    fd = res.get("B7_findings")
    if fd:
        write_csv(os.path.join(d, "findings_check.csv"), ["omega", "r_conf", "findings", "rel_dev", "boltzmann", "fluxbal"],
                  [[x["omega"], x["r_conf"], x["findings"], x["rel_dev"], x["boltzmann"], x["fluxbal_max"]] for x in fd["rows"]])
        O.add("O-B7-findings", "B7", "info", "FINDINGS.md table (N = 2, k_op = 3, delta3 = 0): conversion ratio |S[Rv<-Ru]|^2/|S[Ru<-Ru]|^2 against "
              "the recorded 0.47966 ... 0.50422", fd["max_rel_dev"], 0.0, None, note="max relative deviation; FINDINGS unitarity ~0.4%")
        flat = [x["r_conf"] for x in fd["rows"] if x["r_conf"] is not None]
        if flat:
            O.add("O-B7-flat", "B7", "claim", "manuscript Sec. IV: the conversion ratio at k_op = 3, delta3 = 0 is nearly frequency-independent "
                  "(spread < 10% over omega in [0.2, 1.6], against a Boltzmann factor falling 80-fold)", (max(flat) - min(flat)) / np.mean(flat), 0.10,
                  (max(flat) - min(flat)) / np.mean(flat) < 0.10)
    conv = []
    for N, kop in Sc["conv"]:
        r = res.get("B7_conv_N%s_k%s" % (_f(N), _f(kop)))
        if r:
            conv.append(r)
    if conv:
        write_csv(os.path.join(d, "scatter_convergence.csv"), ["N", "k_op", "h", "conv_h", "L", "conv_L", "S2_dev_h_half", "nH_rel_h_half",
                                                               "S2_dev_L_wide", "nH_rel_L_wide", "fluxbal_max"],
                  [[c["N"], c["k_op"], c["h"], c["conv_h"], c["L"], c["conv_L"], c["h_half"]["S2_max_abs_dev"], c["h_half"]["n_H_max_rel_dev"],
                    c["L_wide"]["S2_max_abs_dev"], c["L_wide"]["n_H_max_rel_dev"], c["fluxbal_max"]] for c in conv])
        worst = max(max(c["h_half"]["S2_max_abs_dev"], c["L_wide"]["S2_max_abs_dev"]) for c in conv)
        O.add("O-B7-conv", "B7", "info", "|S|^2 change under h -> h/2 and L -> %g at %d points" % (Sc["conv_L"], len(conv)), worst, None, None)
    anchor = [p for p in pts if _near(p["N"], cfg["anchor"]["N"]) and _near(p["k_op"], cfg["anchor"]["k_op"])]
    S["B7"] = {"n_points": len(pts), "findings": fd, "kaup_null": kn, "convergence": conv, "anchor": anchor[0] if anchor else None,
               "two_sided_points": [[p["N"], p["k_op"], p["frac_two_sided"]] for p in pts if p["frac_two_sided"] > 0]}


# ============================================================ B8
def red_B8(rundir, cfg, res, O, S):
    d = _mk(rundir, "B8_ref25")
    idt = res.get("B8_ext_identity")
    if idt:
        write_json(os.path.join(d, "ext_identity.json"), idt)
        O.add("O-B8-ext-identity", "B8", "hard", "gnlse_ext = src.gnlse_dimensionless for the four pulses of the script (bitwise)", idt["all_equal"],
              True, idt["all_equal"])
    summ = []
    for name in cfg["ref25"]["variants"]:
        v = grid.ref25_variant(name)
        rows = []
        for key in sorted(k for k in res if k.startswith("B8_%s_" % name)):
            rows += res[key]["rows"]
        if not rows:
            continue
        s, rows = b_ref25.summarize(v, rows)
        write_csv(os.path.join(d, "ref25_%s.csv" % name), ["N_sol", "delta3", "pulse_shape", "mask_width", "mask_order", "DeltaS_tot",
                                                            "eta_GSL_final", "pass", "photon_drift_max", "wall_s"],
                  [[r["N_sol"], r["delta3"], r["pulse_shape"], r["mask_width"], r["mask_order"], r["DeltaS_tot"], r["eta_GSL_final"], r["passed"],
                    r["photon_drift_max"], r["wall_s"]] for r in rows])
        expected = len(v["N"]) * len(v["d3"]) * len(v["shapes"]) * len(v["masks"])
        s["complete"] = len(rows) == expected
        if name == "V1_shipped":
            cmp = b_ref25.compare_reference(rows, os.path.join(REF, "robustness_table.csv"))
            s["vs_530047"] = cmp
            O.add("O-B8-V1-530047", "B8", "hard", "shipped script (V1) against job 530047 S7, 300 configurations row by row",
                  [cmp["n_matched"], cmp["DeltaS_max_abs_dev"], cmp["eta_max_abs_dev"], cmp["n_within_1e-6"]], [300, 0, 0, 300],
                  cmp["n_matched"] == 300 and cmp["DeltaS_max_abs_dev"] <= 1e-3 and cmp["eta_max_abs_dev"] <= 1e-4, "1e-3 / 1e-4",
                  note="[matched, max |dDeltaS|, max |deta|, rows within 1e-6]; round-off of another FFT/BLAS build grows in the N = 4.5 "
                       "fission runs; on the Altay module the rows should agree to round-off")
        allpass = s["pass_count"] == s["n_configs"]
        O.add("O-B8-%s-positive" % name, "B8", "claim" if not name.endswith("fast") else "info",
              "Ref. [25]: every configuration has DeltaS_tot > 0 and eta_GSL > 1/2 (%s)" % name, "%d of %d" % (s["pass_count"], s["n_configs"]),
              "all", allpass if not name.endswith("fast") else None)
        inr = s["n_DeltaS_below"] + s["n_DeltaS_above"] + s["n_eta_below"] + s["n_eta_above"] == 0
        O.add("O-B8-%s-ranges" % name, "B8", "claim" if not name.endswith("fast") else "info",
              "Ref. [25] Eq. (14): DeltaS_tot in [0.82, 2.42], eta_GSL in [0.52, 0.62] (%s)" % name,
              {"DeltaS": s["DeltaS_range"], "eta": s["eta_range"]}, {"DeltaS": [0.82, 2.42], "eta": [0.52, 0.62]},
              inr if not name.endswith("fast") else None)
        summ.append(s)
    write_json(os.path.join(d, "ref25_variants.json"), {"variants": summ, "paper1_ranges": b_ref25.RANGES_PAPER1, "status": "RECORDED"})
    if not summ and not idt:
        O.missing("O-B8", "B8", "hard", "Ref. [25] variants", "no B8 result")
    S["B8"] = {"variants": summ, "ext_identity": idt["all_equal"] if idt else None}


# ============================================================ B10
def red_B10(rundir, cfg, res, O, S):
    d = _mk(rundir, "B10_anchor")
    dd = _mk(rundir, os.path.join("B10_anchor", "paper2_data_anchor"))
    A = cfg["anchor"]
    N, d3, kg, kop = A["N"], A["d3"], A["kg"], A["k_op"]
    kap, aux = suite.kinematic_kappa(N, d3, kg, 1.0, kop)
    om = np.asarray(cfg["omega"], dtype=float)
    a2, b2, ratio = suite.bogoliubov_thermal(om, kap)
    phase2.write_csv(os.path.join(dd, "thermality_ratio.csv"), ["omega", "ln_ratio", "ratio", "alpha2", "beta2"],
                     np.column_stack([om, np.log(ratio), ratio, a2, b2]))
    drive = np.asarray(cfg["drive"], dtype=float)
    kd = np.array([suite.kinematic_kappa(N, d3, kg, x, kop)[0] for x in drive])
    phase2.write_csv(os.path.join(dd, "kappa_drive_independence.csv"), ["drive", "kappa_kin"], np.column_stack([drive, kd]))
    kgf = np.asarray(cfg["kg_flow"], dtype=float)
    kf = np.array([suite.kinematic_kappa(N, d3, g, 1.0, kop)[0] for g in kgf])
    phase2.write_csv(os.path.join(dd, "kinematic_kappa_flow.csv"), ["kappa_g", "kappa_kin", "T_H"], np.column_stack([kgf, kf, kf / (2 * np.pi)]))
    EN, ENx, num = suite.partner_logneg(om, kap)
    phase2.write_csv(os.path.join(dd, "partner_log_negativity.csv"), ["omega", "E_N", "E_N_cross", "nu_minus"], np.column_stack([om, EN, ENx, num]))
    slope, icpt, R2, kfit = suite.thermality_fit(om, ratio)
    mu = 0.5 * N ** 2
    closed = N * float(suite.slow_light(kop, d3, kg))
    cen = [r for r in (S.get("_B4_rows") or []) if _near(r["N"], N) and _near(r["d3"], d3) and _near(r["k_op"], kop) and _near(r["kg"], kg)]
    red3 = (S.get("_B3_red") or {}).get("%s|%s" % (_f(N), _f(d3)))
    bdg = _find(S.get("_B5_rows") or [], N=N, d3=d3)
    sc = [p for p in (S.get("_B7_pts") or []) if _near(p["N"], N) and _near(p["k_op"], kop)]
    HP = 2 * math.pi / kap
    loss_hp = None
    b = res.get("B3_loss_N%s_d%s_base" % (_f(N), _f(d3)))
    if b:
        from sothe_p2 import compat as C
        xi, Q = np.asarray(b["xi"]), np.asarray(b["Q"])
        loss_hp = float(1 - C.interp1(xi, Q, 10.0 + HP) / C.interp1(xi, Q, 10.0))
    summ = {"anchor": A, "kappa_kin": kap, "kappa_closed": closed, "deficit_rel": kap / closed - 1, "T_H": kap / (2 * math.pi), "mu": mu,
            "hawking_period": HP, "thermality_slope": slope, "thermality_R2": R2, "kappa_fit": kfit,
            "flux_residual": suite.flux_residual(a2, b2), "drive_spread": float(kd.max() - kd.min()),
            "EN_band_max": float(EN[0]), "EN_at_mu": b_census.EN_thermal(mu, kap), "EN_at_1p6": b_census.EN_thermal(1.6, kap),
            "n_at_mu": b_census.beta2_thermal(mu, kap),
            "census": [{"frame": c["frame"], "epoch": c["epoch_window"], "Vd": c["Vd"], "Phi_pair": c["Phi_pair"], "Phi_up_low": c["Phi_up_low"],
                        "omega_max_pair": c["omega_max_pair_used"], "window": [c["window_lo"], c["window_hi"]], "window_frac": c["window_frac"]}
                       for c in cen],
            "loss": ({"loss_0_40": red3["loss_0_40"], "loss_first_hp_kop1": red3["loss_first_hp"], "loss_one_anchor_HP_after_10": loss_hp,
                      "window_drift": red3["window_drift"]} if red3 else None),
            "bdg": ({k: bdg[k] for k in ("annulus_maxIm", "zero_sector_max", "tail_rel", "dim_ker_M", "dim_ker_M2", "gap_over_mu")} if bdg else None),
            "scatter_delta3_0": sc[0] if sc else None, "status": "RECORDED (HR-3); closed forms and recorded-protocol numerics only"}
    write_json(os.path.join(d, "anchor_summary.json"), summ)
    O.add("O-B10-kappa", "B10", "claim", "REVIEW_G5 S1-1 A: kappa = 0.9541, T_H = 0.1519, Hawking period 6.585 at (N, k_op) = (1, 1.5)",
          [kap, kap / (2 * math.pi), HP], [0.9541, 0.1519, 6.585],
          abs(kap - 0.9541) <= 5e-5 and abs(kap / (2 * math.pi) - 0.1519) <= 5e-5 and abs(HP - 6.585) <= 5e-4)
    O.add("O-B10-EN", "B10", "claim", "REVIEW_G5 S1-1 A: E_N = 0.390 at omega = mu, 0.0103 at omega = 1.6", [summ["EN_at_mu"], summ["EN_at_1p6"]],
          [0.390, 0.0103], abs(summ["EN_at_mu"] - 0.390) <= 5e-4 and abs(summ["EN_at_1p6"] - 0.0103) <= 5e-5)
    eps = [c for c in cen if c["frame"] != "lab"]
    if eps:
        ok = all(_near(c["window_lo"] or -1, 0.5) and _near(c["window_hi"] or -1, 1.6) for c in eps)
        ok2 = all(c["Phi_pair"] >= 0.999 for c in eps)
        lab = [c for c in cen if c["frame"] == "lab"]
        O.add("O-B10-window", "B10", "claim", "REVIEW_G5 S1-1 A: in the comoving frame the pair channel is open over the whole band and the "
              "Hawking window is [0.5, 1.6] (71% of the band), at every drift epoch",
              {"windows": [[c["frame"], c["window_lo"], c["window_hi"]] for c in eps], "Phi_pair_min": min(c["Phi_pair"] for c in eps),
               "lab": [lab[0]["window_lo"], lab[0]["window_hi"], lab[0]["Phi_pair"]] if lab else None},
              {"window": [0.5, 1.6], "Phi_pair": 1.0}, ok and ok2, note="the laboratory frame (V_d = 0) is listed for comparison")
    S["B10"] = summ


# ============================================================ timing, summary, driver
def timing(rundir, A, B):
    from .tasks import task_paths
    by = {}
    for t in list(A) + list(B):
        p = task_paths(rundir, t)[0]
        j = read_json(p) if os.path.exists(p) else None
        if j is None:
            continue
        b = by.setdefault(t["block"], {"n": 0, "wall_s": 0.0, "max_s": 0.0, "est_s": 0.0})
        b["n"] += 1
        b["wall_s"] += j["wall_s"]
        b["max_s"] = max(b["max_s"], j["wall_s"])
        b["est_s"] += t["cost"]
    return by


def reduce_all(rundir, cfg, blocks, res, A, B, failed, wall_ab=0.0):
    O = Ledger()
    S = {}
    order = [("B1", red_B1), ("B9", red_B9), ("B3", red_B3), ("B4", red_B4), ("B2", red_B2), ("B5", red_B5), ("B6", red_B6), ("B7", red_B7),
             ("B8", red_B8), ("B10", red_B10)]
    errors = {}
    for b, fn in order:
        if b not in blocks:
            continue
        try:
            fn(rundir, cfg, res, O, S)
        except Exception:
            import traceback
            errors[b] = traceback.format_exc()
            O.missing("O-%s-reduce" % b, b, "hard", "reduction of block %s" % b, "exception in phase C (see summary.json reduce_errors)")
    try:
        from . import figs
        S["figures"] = figs.make_all(rundir, cfg, S, res)
    except Exception:
        import traceback
        errors["figures"] = traceback.format_exc()
        S["figures"] = []
    tm = timing(rundir, A, B)
    hard_fail = O.hard_failures()
    code = 2 if (hard_fail or errors) else 0
    public = {k: v for k, v in S.items() if not k.startswith("_")}
    write_json(os.path.join(rundir, "ORACLES.json"), {"rule": RULE, "written": utc(), "rows": O.rows,
                                                      "counts": {c: {"pass": sum(1 for r in O.rows if r["class"] == c and r["pass"] is True),
                                                                     "fail": sum(1 for r in O.rows if r["class"] == c and r["pass"] is False),
                                                                     "n/a": sum(1 for r in O.rows if r["class"] == c and r["pass"] is None)}
                                                                 for c in ("hard", "claim", "info")}})
    write_json(os.path.join(rundir, "summary.json"), {"rule": RULE, "written": utc(), "fast": cfg["fast"], "blocks": blocks, "failed_tasks": failed,
                                                      "reduce_errors": errors, "timing": tm, "wall_phase_AB_s": wall_ab, "results": public,
                                                      "exit_code": 3 if failed else code})
    from .summary_md import write_summary
    write_summary(rundir, cfg, blocks, S, O, failed, errors, tm, wall_ab)
    return code
