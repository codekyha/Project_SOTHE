"""Tests of SOTHE-P2-SWEEP (about a minute on one core).  Run:  python tests/test_sweep.py   (or: python -m pytest tests)

They check the identities the sweep relies on, the new numerics against exact results, and the command line."""
import json
import math
import os
import shutil
import subprocess
import sys
import tempfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
os.environ.setdefault("MPLBACKEND", "Agg")

import numpy as np  # noqa: E402

from sothe_p2 import kit, suite  # noqa: E402
from sweep_p2 import b_bdg, b_census, b_kin, b_loss, b_raman, b_ref25, b_scatter, cli, grid, tasks  # noqa: E402


def test_blocks_consistent():
    assert [b[0] for b in cli.BLOCK_INFO] == list(grid.ALL_BLOCKS)
    assert grid.expand_blocks(["B10"]) == ["B1", "B2", "B3", "B4", "B5", "B7", "B10"]
    assert grid.expand_blocks(["B6"]) == ["B6"]
    assert set(cli.DEFAULTS["BLOCKS"].split()) == set(grid.ALL_BLOCKS)


def test_plan_ids_unique():
    for fast in (False, True):
        cfg = grid.config(fast)
        A = tasks.plan_phase_a(cfg, list(grid.ALL_BLOCKS))
        ids = [t["id"] for t in A]
        assert len(ids) == len(set(ids)), "duplicate task ids"
        assert all(t["fn"] in tasks.FUNCS for t in A)
        json.dumps(A)                                                      # the plan is plain data
    cfg = grid.config(False)
    A = tasks.plan_phase_a(cfg, list(grid.ALL_BLOCKS))
    n8 = sum(1 for t in A if t["block"] == "B8" and t["fn"] == "t_ref25")
    assert n8 == 60 + 60 + 15 + 15 + 60 + 60                               # (N x shape x mask) chunks of 5 delta3 values


def test_kinematic_identity():
    a = suite.kinematic_kappa(2.0, 0.05, 0.3, 1.0, 1.0)[0]
    b = b_kin.kappa_extract(2.0, 0.05, 0.3, 1.0, 1.0, 4096)[0]
    assert a == b and abs(a - 1.6635880450803207) < 1e-12
    pred = b_kin.eq7pp(2.0, 0.05, 0.3, 1.0, 1.0)[0]
    assert abs((b / (2.0 * float(suite.slow_light(1.0, 0.05, 0.3))) - 1) - pred) < 5e-7


def test_parity_published():
    p = b_kin.parity({"odd": [4095, 8191], "even": [4096, 8192], "cont_N": [1.0, 2.0, 3.0]}, grid.OP)
    assert all(abs(x["dev_published"]) < 5e-7 for x in p["table2"])
    for x in p["table3"]:
        assert abs(x["extracted"] / x["published_extracted"] - 1) < 1e-4


def test_census_window_closed_form():
    w = b_census.window(1.0, 1.5, 0.05, 0.3, 1.76)
    assert abs(w["window_lo"] - 0.5) < 1e-15 and abs(w["window_hi"] - 1.6) < 1e-15
    assert abs(w["EN_at_mu"] - 0.390) < 5e-4
    assert b_census.window(2.0, 1.5, 0.05, 0.3, 1.611)["window_lo"] is None      # mu = 2 > 1.6
    c = b_census.census_one(2.0, 1.0, 0.05, 0.0, 0.3, 311)
    r = kit.channel_census(2.0, 1.0, 0.05, 0.0, 0.3)
    assert c["Phi_pair"] == r["Phi_pair"] and c["omega_max_pair"] == r["omega_max_pair"]


def test_scatter_kaup_null_flux():
    kn = b_scatter.kaup_null([0.8, 2.0], [0.02])
    assert max(x["R2"] for x in kn["kaup"]) < 1e-6
    assert max(abs(x["T2"] - 1) for x in kn["kaup"]) < 1e-6
    assert max(x["uv_mix_max"] for x in kn["null"]) == 0.0
    k0, c, _ = b_scatter.kinematics(1.0, 1.5, 0.3)
    S, info = b_scatter.solve(1.0, c, k0, 1.0, L=20.0, h=0.02)
    assert info["ins"] == ["Lu", "Ru", "Rv"] and info["outs"] == ["Lu", "Ru", "Rv"]
    assert max(abs(f) for f in info["fluxbal"]) < 1e-9


def test_scatter_findings():
    r = b_scatter.findings_check(2.0, 3.0, 0.3, [0.2, 1.6], [0.47966, 0.50422], 25.0, 0.01)
    assert r["max_rel_dev"] < 0.02


def test_ref25_extended_solver_identity():
    assert b_ref25.ext_identity(Nt=2 ** 9, n_steps=80, xi_max=1.0)["all_equal"]
    v = grid.ref25_variant("V3_sg8_shipped")
    assert v["shapes"] == ["super_gaussian_m4"] and v["Nt"] == 2 ** 13
    tau = np.linspace(-2, 2, 5)
    assert np.allclose(b_ref25._initial("super_gaussian_m4", tau, 3.0), 3.0 * np.exp(-tau ** 8 / 2))
    c = b_ref25.configs(grid.ref25_variant("V1_shipped"))
    assert len(c) == 300 and (c[0]["N_sol"], c[1]["N_sol"], c[5]["delta3"]) == (2.5, 3.0, c[5]["delta3"])
    assert c[5]["iD"] == 1 and c[25]["iS"] == 1 and c[100]["iM"] == 1           # MATLAB ndgrid order


def test_raman_gordon_and_zero_limit():
    g = b_raman.gordon_check(tauR_lin=0.02, xi_end=5.0, dxi=0.004, n=2048, L=100.0)
    assert abs(g["ratio"] - 1) < 0.05
    assert abs(b_raman.M1_FS - 8.122) < 1e-3 and abs(b_raman.F_R * b_raman.M1_FS - 1.462) < 1e-3
    X = b_raman.rk4ip(2.0, 0.05, 12.0, 0.002, 1024, 100.0)                      # f_R = s = 0
    Y = kit.splitstep_loss(2.0, 0.05, 12.0, 0.002, 1024, 100.0)
    assert abs((1 - X["Q"][-1] / X["Q"][0]) / (1 - Y["Q"][-1] / Y["Q"][0]) - 1) < 1e-3
    assert np.max(np.abs(X["tau_c"] - Y["tau_c"])) < 1e-3


def test_raman_tracking_exact_roll():
    X = b_raman.rk4ip(1.0, 0.0, 2.0, 0.004, 1024, 100.0, tauR_lin=0.02, track=True)
    Y = b_raman.rk4ip(1.0, 0.0, 2.0, 0.004, 1024, 100.0, tauR_lin=0.02, track=False)
    assert np.allclose(X["tau_c"], Y["tau_c"]) and X["n_roll"] == 0


def test_bdg_small():
    s, _ = b_bdg.bdg_point(1.0, 0.02, 160, 16.0)
    assert s["newton_converged"] and s["annulus_maxIm"] < 1e-8 and s["dim_ker_M"] == 2 and s["dim_ker_M2"] == 4


def test_loss_reduction_runs():
    P = {"xi_end": 40.0, "dxi": 0.004, "n": 1024, "L": 100.0}
    base = b_loss.propagate(1.0, 0.05, "base", P)
    r = b_loss.reduce_point(base)
    assert len(r["window_drift"]) == 5 and 0 < r["loss_0_40"] < 1e-2 and math.isnan(r["rel_dev"])


def test_cli_dry_run_and_help():
    tmp = tempfile.mkdtemp(prefix="sweep_test_")
    try:
        env = dict(os.environ, SWEEP_RUNS=os.path.join(tmp, "runs"), SWEEP_RESULTS_DIR=tmp)
        out = subprocess.run([sys.executable, os.path.join(ROOT, "sweep.py"), "help"], stdout=subprocess.PIPE, universal_newlines=True, env=env)
        assert "submit" in out.stdout and out.returncode == 0
        out = subprocess.run([sys.executable, os.path.join(ROOT, "sweep.py"), "submit", "--dry-run", "--blocks", "B7"], stdout=subprocess.PIPE,
                             stderr=subprocess.STDOUT, universal_newlines=True, env=env)
        assert out.returncode == 0, out.stdout
        runs = os.listdir(os.path.join(tmp, "runs"))
        assert len(runs) == 1 and runs[0].startswith("dryrun_")
        js = open(os.path.join(tmp, "runs", runs[0], "job.sbatch")).read()
        for s in ("#SBATCH -A modfkt", "#SBATCH -p cpu2dq", "#SBATCH -c 128", "module load ANACONDA/Anaconda3-2024.06-1-python-3.12",
                  "OMP_NUM_THREADS=1", "sweep.py _job", "--blocks B7"):
            assert s in js, s
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def _all_tests():
    return [(k, v) for k, v in sorted(globals().items()) if k.startswith("test_") and callable(v)]


if __name__ == "__main__":
    import time
    bad = 0
    for name, fn in _all_tests():
        t0 = time.time()
        try:
            fn()
            print("PASS  %-40s %5.1f s" % (name, time.time() - t0), flush=True)
        except Exception as e:
            bad += 1
            print("FAIL  %-40s %5.1f s  %s: %s" % (name, time.time() - t0, type(e).__name__, e), flush=True)
    print("%d tests, %d failed" % (len(_all_tests()), bad))
    sys.exit(1 if bad else 0)
