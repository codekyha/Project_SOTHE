"""Small shared helpers: JSON/CSV writers, atomic files, hashing, provenance."""
import hashlib
import json
import math
import os
import platform
import socket
import sys
import tempfile
from datetime import datetime, timezone

import numpy as np

from . import PACKAGE, __version__

RULE = ("HR-3: only MATLAB R2025b output can be canonical. Every number of this sweep is RECORDED Python output. "
        "Scattering results (block B7) are not Paper-2 claims (INV5); they document D-06 and serve Paper 3.")


def utc(compact=False):
    t = datetime.now(timezone.utc)
    return t.strftime("%Y%m%dT%H%M%SZ" if compact else "%Y-%m-%dT%H:%M:%SZ")


def jsonable(x):
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        if np.iscomplexobj(x):
            return {"re": jsonable(x.real.tolist()), "im": jsonable(x.imag.tolist())}
        return jsonable(x.tolist())
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, np.integer):
        return int(x)
    if isinstance(x, (float, np.floating)):
        v = float(x)
        return v if math.isfinite(v) else None
    if isinstance(x, complex):
        return {"re": jsonable(x.real), "im": jsonable(x.imag)}
    return x


def write_json(path, obj, pretty=True):
    """Atomic JSON write (temporary file in the same folder, then rename)."""
    d = os.path.dirname(os.path.abspath(path))
    os.makedirs(d, exist_ok=True)
    fd, tmp = tempfile.mkstemp(prefix=".tmp_", suffix=".json", dir=d)
    with os.fdopen(fd, "w", encoding="utf-8", newline="\n") as f:
        f.write(json.dumps(jsonable(obj), indent=1 if pretty else None, allow_nan=False, ensure_ascii=True) + "\n")
    os.replace(tmp, path)


def read_json(path, default=None):
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError):
        return default


def write_csv(path, header, rows):
    """rows: list of sequences; floats at %.17g, None/NaN as empty-free 'nan'."""
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)

    def fmt(v):
        if v is None:
            return "nan"
        if isinstance(v, (bool, np.bool_)):
            return "1" if v else "0"
        if isinstance(v, (int, np.integer)):
            return "%d" % v
        if isinstance(v, (float, np.floating)):
            return "%.17g" % float(v) if math.isfinite(float(v)) else ("nan" if math.isnan(float(v)) else ("inf" if v > 0 else "-inf"))
        return str(v).replace(",", ";")

    with open(path, "w", encoding="ascii", errors="replace", newline="\n") as f:
        f.write(",".join(header) + "\n")
        for r in rows:
            f.write(",".join(fmt(v) for v in r) + "\n")


def read_csv(path):
    """Header and rows (strings) of a comma-separated file."""
    with open(path, encoding="utf-8") as f:
        lines = [l.rstrip("\r\n") for l in f if l.strip()]
    head = lines[0].split(",")
    return head, [l.split(",") for l in lines[1:]]


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def versions():
    v = {"python": platform.python_version(), "numpy": np.__version__}
    for mod in ("scipy", "matplotlib"):
        try:
            v[mod] = __import__(mod).__version__
        except Exception:
            v[mod] = None
    return v


def blas():
    try:
        d = np.show_config(mode="dicts").get("Build Dependencies", {})
        return {k: "%s %s" % (d[k].get("name", ""), d[k].get("version", "")) for k in ("blas", "lapack") if k in d}
    except Exception:
        return {}


def provenance(fast=False, workers=1):
    return {"package": PACKAGE, "version": __version__, "utc": utc(), "host": socket.gethostname(),
            "platform": platform.platform(), "python": sys.executable, "versions": versions(), "blas": blas(),
            "slurm_job_id": os.environ.get("SLURM_JOB_ID", ""), "workers": workers, "cpu_count": os.cpu_count(),
            "threads": {k: os.environ.get(k) for k in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS")},
            "fast_mode": bool(fast), "status": "RECORDED (HR-3)", "rule": RULE}


def fnum(x, fmt="%.6g"):
    try:
        v = float(x) + 0.0                     # + 0.0 turns -0.0 into 0.0
    except (TypeError, ValueError):
        return str(x)
    return "n/a" if not math.isfinite(v) else fmt % v
