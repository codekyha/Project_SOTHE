"""MATLAB semantics of sothe_p2.compat."""
import math
import os
import unittest

import _util  # noqa: F401  (path and thread set-up)
import numpy as np

from sothe_p2 import compat as C


class TestLinspace(unittest.TestCase):
    def test_matches_matlab_grids(self):
        """Every grid column of the MATLAB R2025b CSV products is mlinspace, bit for bit."""
        cases = [("thermality_ratio.csv", "omega", 0.05, 1.6, 60, False), ("kappa_drive_independence.csv", "drive", 0.7, 1.3, 9, False),
                 ("kinematic_kappa_flow.csv", "kappa_g", 0.0, 1.0, 21, False), ("false_positive_map.csv", "kappa_g", 0.0, 1.0, 61, True),
                 ("false_positive_map.csv", "drive", 0.0, 1.5, 61, True), ("stability_map.csv", "N", 0.6, 2.4, 13, True),
                 ("stability_map.csv", "delta3", 0.0, 0.08, 13, True)]
        for name, col, a, b, n, uniq in cases:
            head, M = _util.read_csv(os.path.join(_util.REF, name))
            x = M[:, head.index(col)]
            if uniq:
                x = np.unique(x)
            self.assertTrue(np.array_equal(x, C.mlinspace(a, b, n)), "%s %s" % (name, col))

    def test_numpy_linspace_differs(self):
        self.assertFalse(np.array_equal(np.linspace(0.05, 1.6, 60), C.mlinspace(0.05, 1.6, 60)))

    def test_edges(self):
        self.assertEqual(list(C.mlinspace(1, 5, 1)), [5.0])
        self.assertEqual(C.mlinspace(1, 5, 0).size, 0)
        y = C.mlinspace(-20, 20, 4096)
        self.assertEqual((y[0], y[-1]), (-20.0, 20.0))


class TestRoundAndFormats(unittest.TestCase):
    def test_mround(self):
        self.assertEqual([C.mround(v) for v in (0.5, 1.5, 2.5, -0.5, -2.5, 2.4999)], [1, 2, 3, -1, -3, 2])
        self.assertTrue(np.array_equal(C.mround_array([0.5, 1.5, 2.5, -2.5]), [1, 2, 3, -3]))

    def test_mfmt(self):
        self.assertEqual(C.mfmt("%.2f %g %g", float("nan"), float("inf"), -float("inf")), "NaN Inf -Inf")

    def test_trapz(self):
        x = np.array([0.0, 1.0, 3.0])
        y = np.array([1.0, 2.0, 4.0])
        self.assertEqual(C.mtrapz(x, y), 7.5)
        self.assertTrue(np.array_equal(C.mcumtrapz(x, y), [0.0, 1.5, 7.5]))

    def test_interp_and_conv(self):
        self.assertTrue(math.isnan(C.interp1([0, 1], [0, 1], 1.5)))
        self.assertEqual(C.interp1_extrap([0, 1], [0, 2], 1.5), 3.0)
        self.assertTrue(np.array_equal(C.conv_same([1, 2, 3], [1, 1, 1]), [3, 6, 5]))
        self.assertTrue(np.array_equal(C.mgradient([0, 1, 4, 9], [0, 1, 2, 3]), [1, 2, 4, 5]))


if __name__ == "__main__":
    unittest.main()
