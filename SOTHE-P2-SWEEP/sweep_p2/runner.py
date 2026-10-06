"""The sweep driver: phase A (independent tasks) and phase B (census) in one process pool, then phase C (reduce).

Scheduling: longest task first (estimated cost), one task per worker process, single-threaded BLAS/FFT in every
process, so the numbers do not depend on the worker count.  Every result is written atomically by the worker; a
resumed run (resume RUN) recomputes only the tasks without a result file (failed tasks are retried) and then redoes
phase C.  Exit codes: 0 all tasks ran and every hard oracle passed; 2 a hard oracle failed; 3 a task failed."""
import json
import multiprocessing as mp
import os
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool

from . import grid, tasks
from .util import provenance, read_json, utc, write_json


def log(msg):
    print("[%s] %s" % (time.strftime("%H:%M:%S", time.gmtime()), msg), flush=True)


def _fmt_s(s):
    s = int(max(0, s))
    return "%dh%02dm" % (s // 3600, (s % 3600) // 60) if s >= 3600 else ("%dm%02ds" % (s // 60, s % 60) if s >= 60 else "%ds" % s)


def execute(rundir, tlist, workers, label):
    """Run the tasks without a result file; returns (n_ok, n_failed, failed ids)."""
    todo = [t for t in tlist if not os.path.exists(tasks.task_paths(rundir, t)[0])]
    skipped = len(tlist) - len(todo)
    if skipped:
        log("%s: %d of %d tasks already done (resume), %d to run" % (label, skipped, len(tlist), len(todo)))
    if not todo:
        return 0, 0, []
    todo.sort(key=lambda t: -t["cost"])
    total_cost = sum(t["cost"] for t in todo)
    log("%s: %d tasks, estimated %.0f core-s (%s on %d workers; longest task ~%s)" % (
        label, len(todo), total_cost, _fmt_s(max(total_cost / max(workers, 1), todo[0]["cost"])), workers, _fmt_s(todo[0]["cost"])))
    t0 = time.time()
    done_cost, n_ok, failed = 0.0, 0, []
    last, last_i = 0.0, -1
    step = max(1, len(todo) // 20)

    def report(i, force=False):
        nonlocal last, last_i
        el = time.time() - t0
        if i == last_i:
            return
        if force or (i % step == 0) or (el - last > 120):
            last, last_i = el, i
            frac = done_cost / total_cost if total_cost > 0 else 1.0
            eta = el / frac * (1 - frac) if frac > 0 else float("nan")
            log("  %s: %d/%d done, %d failed, %s elapsed, ~%s left" % (label, i, len(todo), len(failed), _fmt_s(el),
                                                                      _fmt_s(eta) if eta == eta else "?"))

    if workers <= 1 or len(todo) == 1:
        for i, t in enumerate(todo, 1):
            tid, ok, wall, err = tasks.run_task(rundir, t)
            done_cost += t["cost"]
            n_ok += ok
            if not ok:
                failed.append(tid)
                log("  FAILED %s: %s" % (tid, (err or "").strip().splitlines()[-1][:200]))
            report(i)
        report(len(todo), True)
        return n_ok, len(failed), failed
    ctx = mp.get_context("spawn")
    remaining = {t["id"]: t for t in todo}
    try:
        with ProcessPoolExecutor(max_workers=min(workers, len(todo)), mp_context=ctx) as ex:
            fut = {ex.submit(tasks.run_task, rundir, t): t for t in todo}
            for i, f in enumerate(as_completed(fut), 1):
                t = fut[f]
                remaining.pop(t["id"], None)
                try:
                    tid, ok, wall, err = f.result()
                except BrokenProcessPool:
                    raise
                except Exception as e:                       # e.g. a result that could not be pickled
                    tid, ok, err = t["id"], False, "%s: %s" % (type(e).__name__, e)
                done_cost += t["cost"]
                n_ok += ok
                if not ok:
                    failed.append(tid)
                    log("  FAILED %s: %s" % (tid, (err or "").strip().splitlines()[-1][:200]))
                report(i)
    except BrokenProcessPool:
        log("WARNING: the worker pool broke (a worker died: memory?); running the %d unfinished tasks serially" % len(remaining))
        for t in list(remaining.values()):
            if os.path.exists(tasks.task_paths(rundir, t)[0]):
                continue
            tid, ok, wall, err = tasks.run_task(rundir, t)
            n_ok += ok
            if not ok:
                failed.append(tid)
    report(len(todo), True)
    return n_ok, len(failed), failed


def load_results(rundir, tlist):
    """id -> result payload for every task with a result file."""
    out = {}
    for t in tlist:
        p = tasks.task_paths(rundir, t)[0]
        if os.path.exists(p):
            with open(p, encoding="utf-8") as f:
                out[t["id"]] = json.load(f)["result"]
    return out


def run_sweep(rundir, fast=False, blocks=None, workers=1, resume=False):
    """The whole sweep in RUNDIR.  Returns the exit code."""
    from . import reduce as RD
    os.makedirs(rundir, exist_ok=True)
    cfg_path = os.path.join(rundir, "config.json")
    if resume and os.path.exists(cfg_path):
        stored = read_json(cfg_path)
        cfg, blocks, fast = stored["cfg"], stored["blocks"], stored["cfg"]["fast"]
        log("resume: configuration of %s (fast=%s, blocks %s)" % (os.path.basename(rundir), fast, " ".join(blocks)))
    else:
        blocks = grid.expand_blocks(blocks or list(grid.ALL_BLOCKS))
        cfg = grid.config(fast)
        write_json(cfg_path, {"cfg": cfg, "blocks": blocks, "written": utc()})
    prov = provenance(fast, workers)
    prov["blocks"] = blocks
    prov["resume"] = bool(resume)
    hist = read_json(os.path.join(rundir, "provenance.json"), {}) or {}
    runs = hist.get("sessions", []) + [prov]
    write_json(os.path.join(rundir, "provenance.json"), {"sessions": runs, "current": prov})
    log("%s %s: blocks %s, fast=%s, %d workers, %s" % (prov["package"], prov["version"], " ".join(blocks), fast, workers,
                                                        ", ".join("%s %s" % kv for kv in prov["versions"].items())))
    t0 = time.time()
    A = tasks.plan_phase_a(cfg, blocks)
    write_json(os.path.join(rundir, "plan_A.json"), A, pretty=False)
    okA, badA, failA = execute(rundir, A, workers, "phase A")
    resA = load_results(rundir, A)
    loss_red = RD.loss_reductions(cfg, resA) if "B3" in blocks else {}
    B = tasks.plan_phase_b(cfg, blocks, loss_red)
    write_json(os.path.join(rundir, "plan_B.json"), B, pretty=False)
    okB, badB, failB = execute(rundir, B, workers, "phase B")
    resB = load_results(rundir, B)
    failed = failA + failB
    log("phase C: reductions, tables, figures, oracles, summary")
    res = dict(resA)
    res.update(resB)
    code_c = RD.reduce_all(rundir, cfg, blocks, res, A, B, failed, wall_ab=time.time() - t0)
    code = 3 if failed else code_c
    log("sweep finished in %s: %d tasks ran or were present, %d failed; exit code %d" % (
        _fmt_s(time.time() - t0), len(A) + len(B) - len(failed), len(failed), code))
    sys.stdout.flush()
    return code
