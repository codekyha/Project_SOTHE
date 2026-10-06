"""Provenance, HR-3 status and stage records.

Every stage writes provenance.json and stage_status.json.  Under HR-3 only MATLAB R2025b output can be
canonical, so canonical_eligible is always false here: every number of this package is RECORDED.  The
canonical question is left open by the PI (2026-09-26); this module does not decide it."""
import os
import platform
import socket
import sys
import time
from datetime import datetime, timezone

import numpy as np

from . import PACK_ROOT, PACKAGE, PORT_OF, __version__
from . import par

RULE = ("HR-3: only MATLAB R2025b output can be canonical. This package is a Python port: every number it writes "
        "is RECORDED, whatever the settings. The canonical question is open (PI decision).")


def utc(compact=False):
    t = datetime.now(timezone.utc)
    return t.strftime("%Y%m%dT%H%M%SZ" if compact else "%Y-%m-%dT%H:%M:%SZ")


def _versions():
    v = {"python": platform.python_version(), "numpy": np.__version__}
    for mod in ("scipy", "matplotlib"):
        try:
            v[mod] = __import__(mod).__version__
        except Exception:
            v[mod] = None
    return v


def blas_info():
    out = {"blas": "", "lapack": ""}
    try:
        d = np.show_config(mode="dicts").get("Build Dependencies", {})
        for k in ("blas", "lapack"):
            if k in d:
                out[k] = "%s %s" % (d[k].get("name", ""), d[k].get("version", ""))
    except Exception:
        try:
            cfg = np.__config__
            out["blas"] = ",".join(getattr(cfg, "blas_opt_info", {}).get("libraries", []))
            out["lapack"] = ",".join(getattr(cfg, "lapack_opt_info", {}).get("libraries", []))
        except Exception:
            pass
    return out


def provenance(fast=False):
    v = _versions()
    b = blas_info()
    S = {"engine": "Python", "version": v["python"], "release": "NumPy %s" % v["numpy"], "versions": v,
         "blas": b["blas"], "lapack": b["lapack"],
         "threads": {"OMP_NUM_THREADS": os.environ.get("OMP_NUM_THREADS"), "MKL_NUM_THREADS": os.environ.get("MKL_NUM_THREADS"),
                     "OPENBLAS_NUM_THREADS": os.environ.get("OPENBLAS_NUM_THREADS"), "workers": par.workers(),
                     "cpu_count": os.cpu_count(), "slurm_cpus_per_task": os.environ.get("SLURM_CPUS_PER_TASK")},
         "computer": platform.platform(), "host": socket.gethostname(), "utc": utc(),
         "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""), "package": PACKAGE, "pack_version": __version__,
         "runner": os.environ.get("P2_RUNNER", ""), "fast_mode": bool(fast),
         "canonical_release_required": "MATLAB R2025b", "canonical_eligible": False, "rule": RULE, "port_of": PORT_OF,
         "runtime_short": "Python %s / NumPy %s" % (v["python"], v["numpy"])}
    return S


class Status(dict):
    """Stage record (ok = the stage ran; pass = its acceptance criteria held)."""

    def __init__(self, stage, title, fast):
        super().__init__(stage=stage, title=title, ok=False, **{"pass": False}, error="", checks={}, outputs=[], notes=[],
                         provenance=provenance(fast))
        self._t0 = time.time()

    def finish(self, outdir):
        from . import jsonio
        self["wall_s"] = time.time() - self._t0
        jsonio.write(os.path.join(outdir, "provenance.json"), self["provenance"])
        jsonio.write(os.path.join(outdir, "stage_status.json"), dict(self))
        print("\n[%s] %s: %s, %s (%.0f s)" % (self["stage"], self["title"], "ran" if self["ok"] else "FAILED",
                                              "PASS" if self["pass"] else "not passed", self["wall_s"]))
        if self["error"]:
            print("[%s] error: %s" % (self["stage"], self["error"].strip().splitlines()[-1]))
        return self


class Tee:
    """Write to the console and to a stage's console.log at the same time."""

    def __init__(self, f, stream):
        self.f = f                      # an open text file shared by the stdout and stderr tees
        self.stream = stream

    def write(self, s):
        self.stream.write(s)
        self.f.write(s)
        return len(s)

    def flush(self):
        self.stream.flush()
        self.f.flush()

    def isatty(self):
        return False


def pack_root():
    return PACK_ROOT


def python_desc():
    return "%s (%s)" % (sys.executable, platform.python_version())
