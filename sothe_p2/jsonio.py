"""JSON output with the conventions of MATLAB jsonencode: NaN and Inf -> null, logical -> true/false,
vectors -> arrays, matrices -> arrays of rows, structs -> objects (field order kept)."""
import json
import math

import numpy as np


def jsonable(x):
    if isinstance(x, dict):
        return {str(k): jsonable(v) for k, v in x.items()}
    if isinstance(x, (list, tuple)):
        return [jsonable(v) for v in x]
    if isinstance(x, np.ndarray):
        if np.iscomplexobj(x):
            raise TypeError("complex arrays are not written to JSON")
        return jsonable(x.tolist())
    if isinstance(x, (bool, np.bool_)):
        return bool(x)
    if isinstance(x, (np.integer,)):
        return int(x)
    if isinstance(x, (float, np.floating)):
        v = float(x)
        return v if math.isfinite(v) else None
    return x


def dumps(obj, pretty=True):
    return json.dumps(jsonable(obj), indent=1 if pretty else None, allow_nan=False, ensure_ascii=True)


def write(path, obj, pretty=True):
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write(dumps(obj, pretty) + "\n")


def read(path, default=None):
    try:
        with open(path, encoding="utf-8") as f:
            return json.load(f)
    except (OSError, ValueError):
        return default
