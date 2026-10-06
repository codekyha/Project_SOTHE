"""Suite v0.1.0 numerics against the MATLAB R2025b products, and the oracle ledgers.

The ledger tests take about a minute (seven 1000 x 1000 eigenproblems, two Newton solves, one short propagation);
SOTHE_P2_QUICK=1 skips them."""
import os
import unittest

import _util  # noqa: F401
import numpy as np

from sothe_p2 import compat as C
from sothe_p2 import suite

QUICK = os.environ.get("SOTHE_P2_QUICK") == "1"
KAPPA_MATLAB = 1.6635880450803306          # phase2_summary.json of the MATLAB R2025b runs (Updates 1 and 4)


class TestAgainstMatlab(unittest.TestCase):
    def test_kinematic_kappa(self):
        k, _ = suite.kinematic_kappa(2.0, 0.05, 0.3, 1.0, 1.0)
        self.assertLess(abs(k / KAPPA_MATLAB - 1), 1e-13)

    def test_thermality_columns(self):
        head, M = _util.read_csv(os.path.join(_util.REF, "thermality_ratio.csv"))
        om = C.mlinspace(0.05, 1.60, 60)
        a2, b2, ratio = suite.bogoliubov_thermal(om, KAPPA_MATLAB)
        self.assertTrue(np.array_equal(om, M[:, 0]))
        for j, col in ((2, ratio), (3, a2), (4, b2), (1, np.log(ratio))):
            self.assertLess(float(np.max(np.abs(col / M[:, j] - 1))), 1e-13, head[j])

    def test_bdg_rate(self):
        head, M = _util.read_csv(os.path.join(_util.REF, "radiation_rate_table.csv"))
        b = suite.bdg_spectrum(1.0, 0.02, 500, 16.0)
        self.assertLess(abs(b["maxIm"] - M[2, 1]), 1e-9)

    def test_params_env_is_explicit(self):
        P = suite.phase2_params(suite.SMALL_STAB)
        self.assertEqual((P["stab_n_tau"], len(P["stab_N"]), len(P["stab_d3"])), (48, 2, 2))
        self.assertEqual(len(suite.phase2_params({})["stab_N"]), 13)


@unittest.skipIf(QUICK, "SOTHE_P2_QUICK=1")
class TestLedgers(unittest.TestCase):
    def test_ledgers(self):
        from sothe_p2 import oracles
        P = suite.phase2_params(suite.SMALL_STAB)
        R = suite.compute_all(P, _util.TRAJ)
        n, m, _ = oracles.run_oracles(P, R, verbose=False)
        self.assertEqual((n, m), (29, 29))
        P2 = suite.phase2_params_v020(suite.SMALL_STAB)
        n2, m2, L2, _ = oracles.run_oracles_v020(P2, R, _util.TRAJ, fast=True, verbose=False)
        self.assertEqual((n2, m2), (35, 35), [r["id"] for r in L2 if not r["pass"]])


if __name__ == "__main__":
    unittest.main()
