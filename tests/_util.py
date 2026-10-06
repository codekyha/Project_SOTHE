"""Shared set-up of the tests: the pack root on sys.path, one BLAS thread, the reference folders."""
import os
import sys

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("MPLBACKEND", "Agg")
sys.dont_write_bytecode = True

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
if ROOT not in sys.path:
    sys.path.insert(0, ROOT)

REF = os.path.join(ROOT, "reference", "matlab_R2025b_U1")
REC = os.path.join(ROOT, "reference", "recorded_octave840")
TRAJ = os.path.join(ROOT, "data", "entropy_trajectory.csv")


def read_csv(path):
    import numpy as np
    with open(path) as f:
        head = f.readline().strip().split(",")
    M = np.genfromtxt(path, delimiter=",", skip_header=1)
    return head, M.reshape(-1, len(head))


def text(path):
    with open(path, encoding="utf-8", newline="") as f:
        return f.read()
