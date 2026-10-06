"""Run the catalogue (or a part of it) into a results folder, with provenance and a manifest."""
import concurrent.futures as cf
import datetime
import hashlib
import json
import os
import platform
import sys
import time

import numpy as np

from . import __version__, experiments, solver


def sha256_file(p):
    h = hashlib.sha256()
    with open(p, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def versions():
    import numpy
    out = {"python": platform.python_version(), "numpy": numpy.__version__, "platform": platform.platform(),
           "machine": platform.machine(), "processor_count": os.cpu_count(), "sothe_rad": __version__}
    try:
        import scipy
        out["scipy"] = scipy.__version__
    except ImportError:
        out["scipy"] = None
    try:
        import matplotlib
        out["matplotlib"] = matplotlib.__version__
    except ImportError:
        out["matplotlib"] = None
    try:
        cfg = numpy.show_config(mode="dicts")
        bl = cfg.get("Build Dependencies", {}).get("blas", {})
        out["blas"] = "%s %s" % (bl.get("name", "?"), bl.get("version", "?"))
    except Exception:
        out["blas"] = "unknown"
    return out


def _one(args):
    run, folder = args
    cfg = {k: v for k, v in run.items() if k not in ("id", "group", "note")}
    series, info = solver.integrate(cfg)
    d = os.path.join(folder, "runs", run["id"])
    os.makedirs(d, exist_ok=True)
    np.savez_compressed(os.path.join(d, "series.npz"), **series)
    meta = {"id": run["id"], "group": run["group"], "cfg": info["cfg"], "wall_s": round(info["wall_s"], 1), "steps": info["steps"],
            "h": info["h"], "kmax": info["kmax"], "absorb_start": info["absorb_start"], "samples": int(series["xi"].size)}
    with open(os.path.join(d, "meta.json"), "w", encoding="utf-8") as f:
        json.dump(meta, f, indent=1, sort_keys=True)
    return run["id"], info["wall_s"]


def write_manifest(folder):
    rows = []
    for dp, dn, fn in os.walk(folder):
        dn.sort()
        for n in sorted(fn):
            if n == "MANIFEST.sha256":
                continue
            p = os.path.join(dp, n)
            rows.append("%s  %s" % (sha256_file(p), os.path.relpath(p, folder).replace(os.sep, "/")))
    rows.sort(key=lambda r: r.split("  ", 1)[1])
    with open(os.path.join(folder, "MANIFEST.sha256"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(rows) + "\n")
    return len(rows)


def check_manifest(folder):
    p = os.path.join(folder, "MANIFEST.sha256")
    bad, n = [], 0
    listed = set()
    with open(p, encoding="utf-8") as f:
        for ln in f:
            if not ln.strip():
                continue
            h, rel = ln.rstrip("\n").split("  ", 1)
            listed.add(rel)
            n += 1
            q = os.path.join(folder, *rel.split("/"))
            if not os.path.isfile(q) or sha256_file(q) != h:
                bad.append(rel)
    extra = []
    for dp, dn, fn in os.walk(folder):
        for nm in fn:
            rel = os.path.relpath(os.path.join(dp, nm), folder).replace(os.sep, "/")
            if rel != "MANIFEST.sha256" and rel not in listed:
                extra.append(rel)
    return n, bad, extra


def run_catalogue(out_root, sets=("all",), workers=1, stamp=None, command=""):
    runs = experiments.select(experiments.catalogue(), sets)
    if not runs:
        raise SystemExit("no run selected")
    stamp = stamp or datetime.datetime.now(datetime.timezone.utc).strftime("%Y%m%dT%H%M%SZ")
    folder = os.path.join(out_root, stamp)
    if os.path.exists(folder):
        raise SystemExit("%s exists" % folder)
    os.makedirs(os.path.join(folder, "runs"))
    t0 = time.time()
    prov = {"package": "SOTHE-P2-RAD", "version": __version__, "stamp": stamp, "started_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
            "command": command or " ".join(sys.argv), "sets": list(sets), "workers": workers, "versions": versions(),
            "host": platform.node(), "runs": [r["id"] for r in runs]}
    with open(os.path.join(folder, "catalogue.json"), "w", encoding="utf-8") as f:
        json.dump(runs, f, indent=1, sort_keys=True)
    order = sorted(runs, key=experiments.cost_estimate, reverse=True)
    walls = {}
    if workers <= 1:
        for r in order:
            rid, w = _one((r, folder))
            walls[rid] = w
            print("[%6.0f s] %-40s %6.1f s" % (time.time() - t0, rid, w), flush=True)
    else:
        with cf.ProcessPoolExecutor(max_workers=workers) as ex:
            futs = {ex.submit(_one, (r, folder)): r["id"] for r in order}
            for fu in cf.as_completed(futs):
                rid, w = fu.result()
                walls[rid] = w
                print("[%6.0f s] %-40s %6.1f s" % (time.time() - t0, rid, w), flush=True)
    prov["finished_utc"] = datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")
    prov["wall_total_s"] = round(time.time() - t0, 1)
    prov["core_s"] = round(sum(walls.values()), 1)
    with open(os.path.join(folder, "RUN.json"), "w", encoding="utf-8") as f:
        json.dump(prov, f, indent=1, sort_keys=True)
    return folder
