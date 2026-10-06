#!/usr/bin/env python3
"""Unit tests of SOTHE-P2-RAD (python tests/test_rad.py, or python rad.py test).  About 15 seconds on one core.
Prints 'RESULT: <n> tests, <m> failed' as its last line."""
import math
import os
import sys
import unittest

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
sys.dont_write_bytecode = True

from sothe_rad import experiments, model, ode, solver  # noqa: E402

RECORDED_LOSS_N1_D005 = 9.328024929e-04     # SOTHE-P2 1.0.0, UHeM Altay job 530047, g4c_loss_report.json, (N, delta3) = (1, 0.05)


def k_res_lab_direct(N, ks, d3):
    """Largest real root of the laboratory-frame resonance condition beta0(k) + theta - (k - ks) V = 0 of the local soliton,
    theta = -beta0(ks) + N^2/2 and V its velocity, solved as a cubic in k (independent of the scaled root of model.q_res)."""
    V = model.soliton_velocity(N, ks, d3)
    theta = -model.beta0(ks, d3) + N ** 2 / 2
    r = np.roots([d3, -0.5, V, -(theta + ks * V)])
    return float(max(x.real for x in r if abs(x.imag) < 1e-9 * max(1.0, abs(x))))



class TestModel(unittest.TestCase):
    def test_roots_solve_the_cubic(self):
        for dt in (0.03, 0.05, 0.08, 0.12):
            q = model.q_res(dt, "drift")
            self.assertLess(abs(dt * q ** 3 - 0.5 * q ** 2 - dt * q - 0.5), 1e-9 * q ** 3)
            q2 = model.q_res(dt, "nodrift")
            self.assertLess(abs(dt * q2 ** 3 - 0.5 * q2 ** 2 - 0.5), 1e-9 * q2 ** 3)
            self.assertGreater(q, q2)

    def test_small_dt_limit(self):
        for dt in (0.005, 0.01):
            self.assertLess(abs(model.q_res(dt) * 2 * dt - 1), 3 * dt)

    def test_lab_root_maps_to_scaled_root(self):
        N, ks, d3 = 1.1, -0.05, 0.07
        k = model.k_res_lab(N, ks, d3)
        D = 1 - 6 * d3 * ks
        V = model.soliton_velocity(N, ks, d3)
        theta = -model.beta0(ks, d3) + N ** 2 / 2
        self.assertLess(abs(model.beta0(k, d3) + theta - (k - ks) * V), 1e-9 * k ** 3)
        self.assertAlmostEqual(model.delta_tilde(N, ks, d3), N * d3 / D ** 1.5, places=14)

    def test_measured_root_equals_model_root_for_model_inputs(self):
        N, d3 = 1.0, 0.05
        V = model.soliton_velocity(N, 0.0, d3)
        self.assertAlmostEqual(model.k_res_measured(N ** 2 / 2, V, d3), model.k_res_lab(N, 0.0, d3), places=9)

    def test_born_rate_against_quadrature(self):
        """Closed form against the golden rule with a numerical Fourier transform of the source -d3 psi0'''."""
        N, d3 = 1.0, 0.05
        L, n = 80.0, 2 ** 16
        tau = -L / 2 + (L / n) * np.arange(n)
        psi0 = N / np.cosh(N * tau)
        k = 2 * np.pi * np.fft.fftfreq(n, L / n)
        src = np.fft.ifft(-d3 * (1j * k) ** 3 * np.fft.fft(psi0))
        q = model.q_res(N * d3, "drift")
        S = abs(np.sum(src * np.exp(-1j * q * tau)) * (L / n))
        rate = S ** 2 / abs(model.v_rel_scaled(q, N * d3, "drift")) / 2.0
        self.assertLess(abs(rate / model.rate_born_scaled(N * d3, "drift") - 1), 1e-6)


class TestSolver(unittest.TestCase):
    small = {"L": 60.0, "n": 1024, "absorb_strength": 0.0, "sample": 0.5}

    def test_nls_soliton_is_stationary(self):
        # absorb_frac 0.51: every grid point lies in the box |tau| < a, so Qbox is the norm of the whole periodic grid,
        # which the scheme conserves to round-off when there is no absorber
        c = dict(self.small, N=1.0, d3=0.0, xi_end=10.0, absorb_frac=0.51)
        s, info = solver.integrate(c, keep_field=True)
        tau, psi = info["tau"], info["psi_end"]
        self.assertLess(np.max(np.abs(np.abs(psi) - 1 / np.cosh(tau))), 2e-6)
        # exact in exact arithmetic; in floating point the norm drifts by about one unit in the last place per step
        # (2.7e-16 per step here, proportional to the number of steps and independent of the step size)
        self.assertLess(abs(s["Qbox"][-1] / s["Qbox"][0] - 1), 1e-15 * info["steps"])
        self.assertLess(abs(s["Q"][-1] / s["Q"][0] - 1), 1e-9)
        ph = np.unwrap(s["phase"])
        self.assertLess(abs((ph[-1] - ph[4]) / (s["xi"][-1] - s["xi"][4]) - 0.5), 1e-5)

    def test_invariants_with_tod(self):
        c = dict(self.small, N=1.0, d3=0.05, xi_end=5.0, absorb_frac=0.51)      # whole periodic grid, no absorber
        s, info = solver.integrate(c)
        self.assertLess(abs(s["Qbox"][-1] / s["Qbox"][0] - 1), 1e-15 * info["steps"])     # round-off only (see above)
        self.assertLess(abs(s["Mbox"][-1] - s["Mbox"][0]), 1e-9)
        self.assertLess(abs(s["E"][-1] - s["E"][0]) / abs(s["E"][0]), 1e-5)

    def test_moving_frame_is_a_translation(self):
        c0 = dict(self.small, N=1.0, d3=0.05, xi_end=4.0)
        s0, i0 = solver.integrate(dict(c0, V=0.0), keep_field=True)
        s1, i1 = solver.integrate(dict(c0, V=-0.05), keep_field=True)
        k = 2 * np.pi * np.fft.fftfreq(i0["psi_end"].size, i0["h"])
        f0, f1 = np.fft.fft(i0["psi_end"]), np.fft.fft(i1["psi_end"])
        # psi_V(tau) = psi_0(tau + V xi): a translation by V xi = -0.2, i.e. the spectrum times exp(i k V xi)
        self.assertLess(np.max(np.abs(f1 - f0 * np.exp(1j * k * (-0.05) * 4.0))), 1e-9 * np.max(np.abs(f0)))
        # tc of the moving frame plus V xi is tc of the laboratory frame; the parabolic vertex of ln|psi|^2 through three
        # grid points (h = 0.059) carries a sub-grid bias of order 1e-6 that depends on the offset of the peak from the grid
        self.assertLess(abs((s1["tc"][-1] - 0.05 * 4.0) - s0["tc"][-1]), 1e-5)

    def test_amplitude_scaling_is_exact_without_absorber(self):
        a = solver.integrate({"N": 1.0, "d3": 0.05, "L": 120.0, "n": 2048, "dxi": 0.004, "xi_end": 8.0, "sample": 0.5, "absorb_strength": 0.0})[0]
        b = solver.integrate({"N": 2.0, "d3": 0.025, "L": 60.0, "n": 2048, "dxi": 0.001, "xi_end": 2.0, "sample": 0.125, "absorb_strength": 0.0})[0]
        self.assertLess(np.max(np.abs(a["Q"] / a["Q"][0] - b["Q"] / b["Q"][0])), 1e-10)

    def test_reproduces_recorded_loss(self):
        s, _ = solver.integrate({"N": 1.0, "d3": 0.05, "xi_end": 40.0, "sample": 0.5})
        loss = 1 - s["Q"][-1] / s["Q"][0]
        self.assertLess(abs(loss / RECORDED_LOSS_N1_D005 - 1), 1e-6)

    def test_flux_estimator_on_a_synthetic_tail(self):
        c = solver.cfg_full({"N": 1.0, "d3": 0.05})
        tau, k, h, sig, a = solver.grid(c)
        kr = model.k_res_lab(1.0, 0.0, 0.05)
        amp = 3e-5
        ramp = 0.5 * (1 + np.tanh((-tau - 10.0) / 2.0)) * 0.5 * (1 + np.tanh((tau + 135.0) / 2.0))
        psi = 1 / np.cosh(tau) + amp * ramp * np.exp(1j * kr * tau)
        m = solver.measure(psi.astype(complex), c, tau, k, h, a, 0.05)
        self.assertLess(abs(m["Iband"] / amp ** 2 - 1), 2e-3)
        self.assertLess(abs(m["kband"] - kr), 1e-3)


class TestModelODE(unittest.TestCase):
    def test_momentum_balance_along_solution(self):
        """d(k_s Q)/dxi = k_r dQ/dxi at states along a solution, with k_r the root of the laboratory-frame resonance
        condition of the local soliton (k_res_lab_direct), computed independently of the scaled root used by ode.rhs."""
        d3 = 0.10
        rs = lambda dt: 50.0 * dt ** -3 * math.exp(-math.pi * model.q_res(dt))
        T = ode.integrate(d3, 1.98, -0.02, 10.0, 60.0, rs)
        for j in range(0, T["xi"].size, 20):
            Q, ks = T["Q"][j], T["ks"][j]
            dlnQ, dks = ode.rhs(T["xi"][j], [math.log(Q), ks], d3, rs)
            dQ = Q * dlnQ
            dP = ks * dQ + Q * dks
            kr = k_res_lab_direct(Q / (2 * math.sqrt(1 - 6 * d3 * ks)), ks, d3)
            self.assertLess(abs(dP - kr * dQ), 1e-9 * abs(kr * dQ))
        # and the integrated balance: Delta(k_s Q) = integral of k_r dQ (trapezoid on the 0.25 output grid)
        kr_t = np.array([k_res_lab_direct(q / (2 * math.sqrt(1 - 6 * d3 * s_)), s_, d3) for q, s_ in zip(T["Q"], T["ks"])])
        lhs = T["Q"][-1] * T["ks"][-1] - T["Q"][0] * T["ks"][0]
        rhs_int = np.sum(0.5 * (kr_t[1:] + kr_t[:-1]) * np.diff(T["Q"]))
        self.assertLess(abs(lhs - rhs_int), 1e-6 * abs(lhs))

    def test_rate_decays_and_c_is_positive(self):
        d3 = 0.10
        rs = lambda dt: 50.0 * dt ** -3 * math.exp(-math.pi * model.q_res(dt))
        T = ode.integrate(d3, 1.98, -0.02, 10.0, 160.0, rs)
        self.assertTrue(np.all(np.diff(T["R"]) < 0))
        self.assertTrue(np.all(T["c_inst"][5:-5] > 0))

    def test_calibration(self):
        d3 = 0.08
        rs = lambda dt: math.exp(-math.pi * model.q_res(dt))
        s = ode.calibrate_scale(2e-5, 1.99, -0.01, d3, rs)
        self.assertAlmostEqual(ode.rate_now(1.99, -0.01, d3, rs, scale=s), 2e-5, places=15)


class TestCatalogue(unittest.TestCase):
    def test_unique_ids_and_groups(self):
        c = experiments.catalogue()
        self.assertEqual(len({r["id"] for r in c}), len(c))
        self.assertEqual({r["group"] for r in c}, {"table", "prepared", "recoil", "scale", "check"})
        for r in c:
            solver.cfg_full({k: v for k, v in r.items() if k not in ("id", "group")})


if __name__ == "__main__":
    res = unittest.TextTestRunner(verbosity=2).run(unittest.defaultTestLoader.loadTestsFromModule(sys.modules[__name__]))
    n, m = res.testsRun, len(res.failures) + len(res.errors)
    print("RESULT: %d tests, %d failed" % (n, m))
    sys.exit(0 if m == 0 else 1)
