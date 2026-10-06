"""kit.py against the recorded Octave 8.4.0 reports (reference/recorded_octave840): formats, census, layer."""
import json
import math
import os
import unittest

import _util  # noqa: F401
import numpy as np

from sothe_p2 import kit


def _rec(name):
    with open(os.path.join(_util.REC, name)) as f:
        return json.load(f)


class TestTokens(unittest.TestCase):
    def setUp(self):
        self.g4c = _rec("g4c_loss_report_octave840_recorded_v012.json")
        self.pts = self.g4c["points"]
        self.op = [p for p in self.pts if p["N"] == 2 and abs(p["d3"] - 0.05) < 1e-12][0]

    def test_fmt_sci_on_recorded_numbers(self):
        T = self.g4c["tokens"]
        self.assertEqual(kit.fmt_sci(self.op["c_over_xi"], 1), T["G4C_LOSS_C"])
        self.assertEqual(kit.fmt_sci(self.op["loss_first_hp"], 1), T["G4C_LOSS_FIRST_HP"])
        rd = [p["rel_dev_0_40"] for p in self.pts]
        self.assertEqual(kit.fmt_sci(max(rd), 1), T["G4C_LOSS_CONV"])
        t05 = [p for p in self.pts if p["N"] == 1 and abs(p["d3"] - 0.05) < 1e-12][0]
        t10 = [p for p in self.pts if p["N"] == 1 and abs(p["d3"] - 0.10) < 1e-12][0]
        self.assertEqual(kit.fmt_sci(t05["newton_tail_rel"], 1), T["G4C_TAIL_005"])
        self.assertEqual(kit.fmt_sci(t10["newton_tail_rel"], 1), T["G4C_TAIL_010"])

    def test_fmt_sci_rules(self):
        self.assertEqual(kit.fmt_sci(0, 1), "0")
        self.assertEqual(kit.fmt_sci(9.96e-3, 1), "1.0\\times10^{-2}")        # mantissa rounding carries into the exponent
        self.assertEqual(kit.fmt_sci(0.2346, 2), "0.235")                    # plain decimal for 0.1 <= |x| < 1000
        self.assertEqual(kit.fmt_sci(float("nan"), 1), "NaN\\times10^{NaN}")

    def test_census_on_recorded_drifts(self):
        """The census with the recorded (grid-peak) drifts reproduces the recorded epochs."""
        for e in self.g4c["channels"]["epochs"]:
            c = kit.channel_census(2, 1.0, 0.05, e["Vd"], 0.3)
            for k in ("S", "c_match", "two_c", "omega_max_pair", "omega_max_pair_no_tod", "Phi_pair", "Phi_up_low", "k_star"):
                self.assertTrue(math.isclose(c[k], e[k], rel_tol=1e-12, abs_tol=1e-15), (e["window"], k, c[k], e[k]))
            self.assertLess(abs(c["pair_open_below"] - e["pair_open_below"]), 1e-15)

    def test_gauss_layer_tokens(self):
        g4a = _rec("g4a_report_octave840_recorded_v011.json")
        self.assertEqual(kit.gauss_layer(1.664100588676)["tokens"], g4a["T06_tokens"])

    def test_subgrid_peak_exact_for_gaussian(self):
        h = 0.05
        tau = -5 + h * np.arange(201)
        for t0 in (0.013, -0.021, 0.0249):
            I2 = np.exp(-(tau - t0) ** 2 / (2 * 0.3 ** 2))
            ip = int(np.argmax(I2))
            self.assertAlmostEqual(kit.subgrid_peak(I2, ip, tau, h), t0, places=12)

    def test_kaup_orders(self):
        K = kit.kaup_residual(0.7, 8.0, (0.2, 0.1, 0.05, 0.025))
        self.assertTrue(np.allclose(K["res"], K["recorded"], rtol=2e-3))
        self.assertGreaterEqual(float(np.min(K["order"])), 3.7)


if __name__ == "__main__":
    unittest.main()
