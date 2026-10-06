"""The task list of the sweep (phases A and B) and the task functions that the worker processes run.

A task is a plain dict: id, block, fn (a name in FUNCS), kw (keyword arguments, JSON-able), cost (estimated
single-core seconds, used for longest-first scheduling).  Its result is written atomically to
tasks/<block>/<id>.json (or tasks/<block>/<id>.error.json with the traceback); a resumed run skips every task whose
result file exists."""
import os
import time
import traceback

from . import b_bdg, b_census, b_kin, b_loss, b_raman, b_ref25, b_scatter, grid
from .util import write_json


# ============================================================ task functions (module level: picklable under spawn)
def t_kin(N, KOP, D3, KG, drive, kg_flow, n):
    return b_kin.kin_chunk(N, KOP, D3, KG, drive, kg_flow, n)


def t_parity(cfg_par, op):
    return b_kin.parity(cfg_par, op)


def t_loss(N, d3, kind, P):
    return b_loss.propagate(N, d3, kind, P)


def t_bdg(N, d3, n_tau, tau_max):
    summ, spec = b_bdg.bdg_point(N, d3, n_tau, tau_max)
    w, kr = spec["w"], spec["krein"]
    return {"summary": summ, "spectrum": {"re": w.real.tolist(), "im": w.imag.tolist(), "krein": kr.tolist()}}


def t_raman(N, d3, T0_fs, model, P):
    return b_raman.run_model(N, d3, T0_fs, model, P)


def t_gordon(tauR_lin, T0_bw, xi_end, dxi, n, L):
    return b_raman.gordon_check(tauR_lin=tauR_lin, T0_bw=T0_bw, xi_end=xi_end, dxi=dxi, n=n, L=L)


def t_scatter(N, k_op, kg, omegas, L, h):
    return b_scatter.scatter_point(N, k_op, kg, omegas, L, h)


def t_scatter_conv(N, k_op, kg, omegas, L, h, conv_h, conv_L):
    return b_scatter.convergence_point(N, k_op, kg, omegas, L, h, conv_h, conv_L)


def t_kaup_null(omegas, hs, L):
    return b_scatter.kaup_null(omegas, hs, L)


def t_findings(N, k_op, kg, omegas, ratios, L, h):
    return b_scatter.findings_check(N, k_op, kg, omegas, ratios, L, h)


def t_ref25(variant, iN, iS, iM):
    return {"variant": variant, "iN": iN, "iS": iS, "iM": iM, "rows": b_ref25.chunk(grid.ref25_variant(variant), iN, iS, iM)}


def t_ref25_identity():
    return b_ref25.ext_identity()


def t_census(N, d3, KOP, KG, epochs, nw):
    return {"N": N, "d3": d3, "rows": b_census.census_rows(N, d3, KOP, KG, [tuple(e) for e in epochs], nw)}


FUNCS = {f.__name__: f for f in (t_kin, t_parity, t_loss, t_bdg, t_raman, t_gordon, t_scatter, t_scatter_conv, t_kaup_null,
                                  t_findings, t_ref25, t_ref25_identity, t_census)}


# ============================================================ the plan
def _f(x):
    return ("%.4f" % x).rstrip("0").rstrip(".") if isinstance(x, float) else str(x)


def plan_phase_a(cfg, blocks):
    """Phase A: every task that needs nothing but the configuration."""
    T = []

    def add(block, tid, fn, kw, cost):
        T.append({"id": tid, "block": block, "fn": fn, "kw": kw, "cost": float(cost)})

    fast = cfg["fast"]
    if "B1" in blocks:
        for N in cfg["N"]:
            add("B1", "B1_kin_N%s" % _f(N), "t_kin", dict(N=N, KOP=cfg["k_op"], D3=cfg["d3_kin"], KG=cfg["kg"], drive=cfg["drive"],
                                                          kg_flow=cfg["kg_flow"], n=cfg["tau_kin_n"]), 1)
    if "B9" in blocks:
        add("B9", "B9_parity", "t_parity", dict(cfg_par=cfg["parity"], op=cfg["op"]), 1)
    if "B3" in blocks:
        P = cfg["loss"]
        scale = (P["xi_end"] / 40.0) * (P["n"] / 8192.0) * (0.002 / P["dxi"])
        kinds = ["base"] + (["half", "n2"] if P["conv"] else [])
        for N, d3 in P["points"]:
            for kind in kinds:
                c = {"base": 10, "half": 20, "n2": 32}[kind] * scale
                add("B3", "B3_loss_N%s_d%s_%s" % (_f(N), _f(d3), kind), "t_loss", dict(N=N, d3=d3, kind=kind, P=P), c)
    if "B5" in blocks:
        B = cfg["bdg"]
        for N, d3 in B["points"]:
            add("B5", "B5_bdg_N%s_d%s" % (_f(N), _f(d3)), "t_bdg", dict(N=N, d3=d3, n_tau=B["n_tau"], tau_max=B["tau_max"]),
                11 * (B["n_tau"] / 500.0) ** 3)
    if "B6" in blocks:
        R = cfg["raman"]
        scale = (R["xi_end"] / 40.0) * (R["n"] / 8192.0) * (0.002 / R["dxi"])
        for N, d3 in R["points"]:
            done_none = False
            for model in R["models"]:
                for T0 in R["T0_fs"]:
                    if model == "none":
                        if done_none:
                            continue
                        done_none = True
                    c = {"none": 22, "bw_ss": 56, "bw": 38, "lin3": 36, "ss": 41}[model] * scale
                    add("B6", "B6_raman_N%s_d%s_%s_T%s" % (_f(N), _f(d3), model, _f(T0) if model != "none" else "any"), "t_raman",
                        dict(N=N, d3=d3, T0_fs=T0, model=model, P=R), c)
        G = R["gordon"]
        gP = dict(xi_end=G["xi_end"], dxi=R["dxi"], n=R["n"], L=R["L"])
        for tl in G["tauR_lin"]:
            add("B6", "B6_gordon_lin_%s" % _f(tl), "t_gordon", dict(tauR_lin=tl, T0_bw=None, **gP), 9 * G["xi_end"] / 10.0)
        for T0 in G["T0_bw"]:
            add("B6", "B6_gordon_bw_T%s" % _f(T0), "t_gordon", dict(tauR_lin=None, T0_bw=T0, **gP), 12 * G["xi_end"] / 10.0)
    if "B7" in blocks:
        Sc = cfg["scatter"]
        n_scale = (Sc["L"] / 25.0) * (0.01 / Sc["h"])
        for N, kop in Sc["points"]:
            add("B7", "B7_scatter_N%s_k%s" % (_f(N), _f(kop)), "t_scatter",
                dict(N=N, k_op=kop, kg=Sc["kg"], omegas=cfg["omega"], L=Sc["L"], h=Sc["h"]), 0.02 * len(cfg["omega"]) * n_scale)
        for N, kop in Sc["conv"]:
            add("B7", "B7_conv_N%s_k%s" % (_f(N), _f(kop)), "t_scatter_conv",
                dict(N=N, k_op=kop, kg=Sc["kg"], omegas=cfg["omega"], L=Sc["L"], h=Sc["h"], conv_h=Sc["conv_h"], conv_L=Sc["conv_L"]),
                0.1 * len(cfg["omega"]))
        add("B7", "B7_kaup_null", "t_kaup_null", dict(omegas=Sc["kaup_omega"], hs=Sc["kaup_h"], L=25.0), 1)
        F = cfg["findings"]
        add("B7", "B7_findings", "t_findings", dict(N=F["N"], k_op=F["k_op"], kg=Sc["kg"], omegas=F["omega"], ratios=F["ratio"],
                                                   L=25.0, h=0.01), 1)
    if "B8" in blocks:
        add("B8", "B8_ext_identity", "t_ref25_identity", {}, 1)
        for name in cfg["ref25"]["variants"]:
            v = grid.ref25_variant(name)
            per = 1.8 * (v["Nt"] / 8192.0) * (v["n_steps"] / 4000.0)
            for iM in range(len(v["masks"])):
                for iS in range(len(v["shapes"])):
                    for iN in range(len(v["N"])):
                        add("B8", "B8_%s_N%d_S%d_M%d" % (name, iN, iS, iM), "t_ref25", dict(variant=name, iN=iN, iS=iS, iM=iM),
                            per * len(v["d3"]))
    if fast:
        for t in T:
            t["cost"] = max(t["cost"], 0.01)
    return T


def epochs_from_loss(red):
    """Drift epochs of one B3 point: [(label, Vd, [lo, hi]), ...] with the laboratory frame first."""
    ep = [("lab", 0.0, None)]
    for lo, hi, vd in zip(red["window_lo"], red["window_hi"], red["window_drift"]):
        ep.append(("xi%g-%g" % (lo, hi), float(vd), [lo, hi]))
    return ep


def plan_phase_b(cfg, blocks, loss_red):
    """Phase B: the comoving census (B4) with the drift epochs measured in B3; delta3 = 0 in the laboratory frame."""
    T = []
    if "B4" not in blocks:
        return T
    for N in cfg["N"]:
        for d3 in cfg["d3_kin"]:
            key = "%s|%s" % (_f(N), _f(d3))
            if d3 == 0.0:
                ep = [("lab", 0.0, None)]
            elif key in loss_red:
                ep = epochs_from_loss(loss_red[key])
            else:
                ep = [("lab", 0.0, None)]
            T.append({"id": "B4_census_N%s_d%s" % (_f(N), _f(d3)), "block": "B4", "fn": "t_census",
                      "kw": dict(N=N, d3=d3, KOP=cfg["k_op"], KG=cfg["kg"], epochs=ep, nw=cfg["census_nw"]),
                      "cost": 0.0007 * len(cfg["k_op"]) * len(cfg["kg"]) * len(ep) * cfg["census_nw"] / 10.0})
    return T


# ============================================================ one task in a worker
def task_paths(rundir, t):
    d = os.path.join(rundir, "tasks", t["block"])
    return os.path.join(d, t["id"] + ".json"), os.path.join(d, t["id"] + ".error.json")


def run_task(rundir, t):
    """Run one task; write its result (or its traceback).  Returns (id, ok, wall_s, error_text)."""
    out, err = task_paths(rundir, t)
    t0 = time.time()
    try:
        res = FUNCS[t["fn"]](**t["kw"])
        wall = time.time() - t0
        write_json(out, {"id": t["id"], "block": t["block"], "fn": t["fn"], "kw": t["kw"], "ok": True, "wall_s": wall,
                         "pid": os.getpid(), "result": res}, pretty=False)
        if os.path.exists(err):
            os.remove(err)
        return t["id"], True, wall, None
    except Exception:
        wall = time.time() - t0
        tb = traceback.format_exc()
        write_json(err, {"id": t["id"], "block": t["block"], "fn": t["fn"], "kw": t["kw"], "ok": False, "wall_s": wall, "error": tb})
        return t["id"], False, wall, tb
