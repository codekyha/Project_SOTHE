"""Ports of the Ref. [25] MATLAB functions (sothe_p2/ref25/src.py)."""
import math
import unittest

import _util  # noqa: F401
import numpy as np

from sothe_p2.ref25 import src


class TestRef25(unittest.TestCase):
    def test_spectral_entropy_uniform(self):
        k = np.linspace(-3.0, 5.0, 801)
        self.assertAlmostEqual(src.spectral_entropy(k, np.ones_like(k), np.ones_like(k)), math.log(8.0), places=12)

    def test_spectral_entropy_checks(self):
        k = np.linspace(0, 1, 5)
        with self.assertRaises(src.MatlabError):
            src.spectral_entropy(k, np.array([1, 1, -1, 1, 1.0]), np.ones(5))
        with self.assertRaises(src.MatlabError):
            src.spectral_entropy(k, np.ones(4), np.ones(5))
        self.assertEqual(src.spectral_entropy(k, np.ones(5), np.zeros(5)), 0.0)

    def test_eta_gsl(self):
        xi = np.linspace(0, 10, 101)
        self.assertEqual(src.eta_gsl(xi, xi ** 2)[-1], 1.0)
        e = src.eta_gsl(xi, np.sin(8 * xi))
        self.assertTrue(0.35 < e[-1] < 0.65)
        self.assertEqual(e[0], 0.0)
        with self.assertRaises(src.MatlabError):
            src.eta_gsl(xi[::-1], xi)

    def test_soliton_mask(self):
        k = np.linspace(-20, 20, 4001)
        P = np.exp(-(k - 2.5) ** 2)
        mS, mR, kc = src.soliton_mask(k, P, 3.0)
        self.assertAlmostEqual(kc, 2.5, places=10)
        self.assertTrue(np.allclose(mS + mR, 1.0))
        self.assertAlmostEqual(float(mS[np.argmin(np.abs(k - 2.5))]), 1.0, places=6)

    def test_find_rr_lobe(self):
        k = np.linspace(-50, 50, 4096)
        P = np.exp(-(k + 15) ** 2 / 2) + 0.05 * np.exp(-(k - 25) ** 2 / 8)
        self.assertAlmostEqual(src.find_rr_lobe(k, P, -15.0, 4.0), 25.0, places=2)
        self.assertTrue(math.isnan(src.find_rr_lobe(k, P, 40.0, 4.0)))

    def test_rr_phase_match_fit(self):
        d3 = np.array([0.06, 0.07, 0.08, 0.09, 0.10])
        out = src.rr_phase_match_fit(d3, 1.01 / d3)
        self.assertAlmostEqual(out["c_RR"], 1.01, places=12)
        self.assertAlmostEqual(out["R2"], 1.0, places=12)
        self.assertAlmostEqual(out["invariant_spread"], 0.0, places=12)

    def test_solver_matches_validation_port(self):
        """src.gnlse_dimensionless (MATLAB statements) against the unchanged validation port on a small grid."""
        from sothe_p2.ref25 import gnlse_dimensionless as v
        a = src.gnlse_dimensionless(N_sol=3.5, delta3=0.02, xi_max=3.0, n_steps=600, Nt=2048, n_save=40, verbose=False)
        b = v.gnlse_dimensionless(N_sol=3.5, delta3=0.02, xi_max=3.0, n_steps=600, Nt=2048, n_save=40, verbose=False)
        self.assertLess(float(np.max(np.abs(a["S_tot"] - b.S_tot))), 1e-10)
        self.assertLess(abs(a["photon_number"][-1] / a["photon_number"][0] - 1), 1e-3)   # validate_paper_claims.py bound

    def test_dimensional_runs(self):
        r = src.gnlse_dimensional(Nt=2 ** 11, n_steps=300, n_save=20, verbose=False)
        self.assertEqual(r["S_tot"].size, 20)
        self.assertTrue(np.isfinite(r["Delta_S_tot"]))


if __name__ == "__main__":
    unittest.main()
