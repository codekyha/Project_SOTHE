"""Runtime probe: versions, BLAS, the process pool (spawn), a numerical smoke test of every block, and a speed
calibration that turns the task-cost model into a wall-time estimate for the full sweep on this node."""
import math
import os
import time

import numpy as np

from . import b_bdg, b_kin, b_raman, b_ref25, b_scatter, grid, tasks
from .util import provenance, write_json


def _check(rows, name, ok, value, ref):
    rows.append({"check": name, "pass": bool(ok), "value": value, "reference": ref})
    print("  check %-44s %s   value %s   reference %s" % (name, "PASS" if ok else "FAIL", value, ref), flush=True)
    return ok


def _pool_ok(workers):
    import multiprocessing as mp
    from concurrent.futures import ProcessPoolExecutor
    kw = dict(N=1.0, KOP=[1.0], D3=[0.05], KG=[0.3], drive=[1.0], kg_flow=[0.0, 0.3], n=4096)
    with ProcessPoolExecutor(max_workers=max(2, min(workers, 4)), mp_context=mp.get_context("spawn")) as ex:
        r = [f.result() for f in [ex.submit(tasks.t_kin, **kw) for _ in range(4)]]
    return all(x["rows"][0]["kappa_num"] == r[0]["rows"][0]["kappa_num"] for x in r)


def run(rundir, workers=2):
    P = provenance(False, workers)
    print("--- versions %s" % ", ".join("%s %s" % kv for kv in P["versions"].items()))
    print("--- blas     %s" % (", ".join("%s: %s" % kv for kv in P["blas"].items()) or "n/a"))
    print("--- threads  %s   cpus %s" % (P["threads"], P["cpu_count"]))
    rows = []
    ok = True
    t0 = time.time()
    k, _ = b_kin.kappa_extract(2.0, 0.05, 0.3, 1.0, 1.0)
    ok &= _check(rows, "B1 kappa at the operating point", abs(k - 1.6635880450803207) < 1e-12, "%.16g" % k, "1.6635880450803207")
    p = b_kin.parity({"odd": [4095], "even": [4096], "cont_N": [1.0]}, grid.OP)
    ok &= _check(rows, "B9 Table 2 (n = 4095, 4096)", all(abs(x["dev_published"]) < 5e-7 for x in p["table2"]),
                 ["%.6f" % x["deficit_over_lamdt2"] for x in p["table2"]], ["-0.333263", "-0.583065"])
    kappa0, c, _ = b_scatter.kinematics(2.0, 3.0, 0.3)
    S, info = b_scatter.solve(2.0, c, kappa0, 0.6, L=25.0, h=0.01)
    r = b_scatter.reduce_omega(S, info, kappa0)
    ok &= _check(rows, "B7 conversion ratio (2, 3), omega = 0.6", abs(r["r_conf"] / 0.48300 - 1) < 0.02 and r["fluxbal_max"] < 1e-9,
                 "%.5f (flux %.1e)" % (r["r_conf"], r["fluxbal_max"]), "0.48300 (FINDINGS), flux 0")
    kn = b_scatter.kaup_null([1.2], [0.02])
    ok &= _check(rows, "B7 Kaup reflectionless, no-soliton null", kn["kaup"][0]["R2"] < 1e-6 and kn["null"][0]["uv_mix_max"] == 0.0,
                 "|R|^2 %.1e, mix %.1e" % (kn["kaup"][0]["R2"], kn["null"][0]["uv_mix_max"]), "< 1e-6, 0")
    idt = b_ref25.ext_identity(Nt=2 ** 9, n_steps=100, xi_max=1.0)
    ok &= _check(rows, "B8 extended solver = Ref. [25] port (bitwise)", idt["all_equal"], idt["all_equal"], True)
    g = b_raman.gordon_check(tauR_lin=0.02, xi_end=5.0, dxi=0.004, n=2048, L=100.0)
    ok &= _check(rows, "B6 Raman SSFS against Gordon (tau_R = 0.02)", abs(g["ratio"] - 1) < 0.08, "%.4f" % g["ratio"], "1 (8%)")
    s, _ = b_bdg.bdg_point(1.0, 0.02, 160, 16.0)
    ok &= _check(rows, "B5 TOD soliton, BdG (n_tau = 160)", s["newton_converged"] and s["annulus_maxIm"] < 1e-8 and s["dim_ker_M"] == 2,
                 "res %.1e, maxIm %.1e, ker %d" % (s["newton_res"], s["annulus_maxIm"], s["dim_ker_M"]), "converged, < 1e-8, 2")
    try:
        pool = _pool_ok(workers)
    except Exception as e:
        pool = False
        print("  process pool error: %s" % e)
    ok &= _check(rows, "process pool (spawn), identical results", pool, pool, True)
    # speed calibration: one short split-step propagation at the production grid
    from sothe_p2 import kit
    t1 = time.time()
    kit.splitstep_loss(2.0, 0.05, 12.0, 0.002, 8192, 400.0)
    t_cal = (time.time() - t1) * 40.0 / 12.0
    speed = t_cal / 10.0                                   # ratio to the cost model (10 s for a base run of B3)
    cfg = grid.config(False)
    A = tasks.plan_phase_a(cfg, grid.expand_blocks(list(grid.ALL_BLOCKS)))
    tot = sum(t["cost"] for t in A) * speed
    longest = max(t["cost"] for t in A) * speed
    est = {}
    for w in (32, 64, 128):
        est[w] = max(tot / w, longest) + 60.0
    print("  speed factor %.2f (cost model = 1.00); full sweep ~%.0f core-s, %d tasks in phase A" % (speed, tot, len(A)))
    print("  estimate: wall ~%s" % ", ".join("%d min on %d workers" % (math.ceil(est[w] / 60.0), w) for w in (32, 64, 128)))
    write_json(os.path.join(rundir, "probe.json"), {"provenance": P, "checks": rows, "pass": bool(ok), "speed_factor": speed,
                                                    "full_sweep_core_s": tot, "longest_task_s": longest, "wall_estimate_s": est,
                                                    "probe_wall_s": time.time() - t0})
    return bool(ok)
