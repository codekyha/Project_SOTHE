#!/usr/bin/env python3
"""zero_sector_tests.py -- the rank, perturbation, convergence, grid-scaling, scan and dispersion tests behind the section "BdG
discretization and convergence" (S4) of the supplement of the article "Linear spectrum, radiative loss and emission channels of solitons with third-order dispersion underlying optical analogue horizons" (H. Oguz, 2026), and a check of each numerical
statement of that section against them.

Operator. About the stationary third-order-dispersion (TOD) soliton phi the BdG operator is
    M = [[H + T, -diag(phi^2)], [diag(conj(phi)^2), -H + T]],   H = -D2/2 - 2 diag|phi|^2 + mu,   T = -i delta3 D3,   mu = N^2/2,
on the periodic grid tau_j = -tau_max + j 2 tau_max / n_tau (endpoint excluded), tau_max = 16, N = 1, with the fourth-order
central-difference matrices D2 and D3 of SOTHE-P2 (sothe_p2.suite.fd_circulant). It is the operator of
sothe_p2.kit.bdg_spectrum_bg. Its zero sector (phase and translation symmetry) is two 2x2 Jordan blocks, and the four computed
eigenvalues of the sector are the round-off floor of a defective eigenvalue, not a property of phi.

Part A (default): n_tau = 500, delta3 in {0.02, 0.05}. phi is the stationary soliton that SOTHE-P2 1.0.0 computed on UHeM
Altay (job 530047, stage S1, S1_g4a/bdg_backgrounds.npz of the results archive).
  Test 1  backward-error scaling: eigenvalues of M + delta ||M||_2 E, E complex Gaussian with ||E||_2 = 1 (three draws per
          delta, fixed seed). A semisimple zero eigenvalue moves by O(delta ||M||), a 2x2 Jordan block by O(sqrt(delta ||M||)).
          Log-log slope over delta in [1e-13, 1e-8].
  Test 2  rank: the smallest singular values of M and of M^2, relative to the largest.
  Test 3  null vectors: the phase mode (phi, -conj(phi)) and the translation mode (phi', conj(phi')), with the fourth-order first
          derivative; residual ||M v|| / (||M||_2 ||v||).
  Test 4  delta3 = 0, about the integrable profile psi0 = N sech(N tau), n_tau = 200 ... 600: the residual of psi0 in the discretized
          equation, the zero sector, the gap, the annulus and its Krein census.
  Test 5  delta3 = 0.005 ... 0.05, N = 1, again about psi0, which is not stationary for delta3 > 0: the four zero-sector modes split
          into a complex quartet. Its max|Im Omega| against sqrt(delta3), its convergence in n_tau = 200 ... 600 and in the box
          (tau_max = 8 ... 24 at the grid spacing of n_tau = 500), and the phase-matching wavenumber k_RR of the resonant branch.
Part B (--grid): the stationary soliton is recomputed (sothe_p2.kit.tod_soliton) at n_tau = 250 ... 1250, and the spectral radius
rho of M, the Newton residual, the zero sector and the annulus are recorded (a few minutes on two cores).
Part C (default): the numbers of the last three paragraphs of the section, the scan over N and delta3, the census at the operating
point and the dispersion of the stencils. They are read from block B5 of the archived run of SOTHE-P2-SWEEP (the table
B5_bdg/bdg_points.csv over N = 0.6 ... 2 and delta3 = 0, 0.02, 0.05, and the spectrum B5_bdg/spectra/bdg_N1_d0.05.csv at the
operating point) or computed in closed form (the symbols of the two stencils, which are compared with the matrices of sothe_p2, the
discrete branch of the positive-norm modes and its zero crossings). With --grid the operating point of the table is also compared with
the recomputation of Part B.

Each numerical statement of the section that rests on these numbers is then checked and printed; the first check is that the operator
built here is the one of sothe_p2. A printed value is accepted when the recomputed value rounds to it (it lies within half a unit of
its last digit); a printed bound is checked as a bound, and at the round-off level of the arithmetic, which depends on the BLAS
library, with a factor of 2 (2.5 for the null vector of the phase mode); an approximate law ("~", "of order") is accepted within
the tolerance that its line states. A check that holds only within the rounding of a printed value or bound, or only with the factor
allowed at the round-off level, says so (a NOTE line, and "literal": false in the result file). The script does not judge the argument built on the numbers, and it checks numbers, not the
qualitative statements of the section. Exit status 0 when every check holds.

Environment variables (or the options):
  SOTHE_SWEEP_PKG   folder of the SOTHE-P2-SWEEP package, which vendors sothe_p2 (default: ../SOTHE-P2-SWEEP)
  SOTHE_P2_RUN      run folder 20260927T175224Z, unpacked from SOTHE-P2_results_20260927T175224Z.tar.gz
  SOTHE_SWEEP_RUN   run folder 20260928T131100Z, unpacked from SOTHE-P2-SWEEP_results_20260928T131100Z.tar.gz
  SOTHE_ARTICLE_OUT output folder of zero_sector_tests.json (default: ./article_output)
Python 3.9 or newer with NumPy."""
import argparse
import csv
import json
import math
import os
import sys
import time

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
EPS = float(np.finfo(float).eps)
N_TAU, TAU_MAX, N_SOL = 500, 16.0, 1.0
TAGS = (("d3_0p02", 0.02), ("d3_0p05", 0.05))
DELTAS = (0.0, 1e-15, 1e-14, 1e-13, 1e-12, 1e-11, 1e-10, 1e-9, 1e-8)
SEED = 20260927
CONV_N = (200, 300, 400, 500, 600)
GRID_N = (250, 500, 750, 1000, 1250)
Q_DELTAS = (0.005, 0.01, 0.02, 0.05)
BOX_TAU = (8.0, 12.0, 16.0, 20.0, 24.0)
KERM2_BOUND = 4.7e-8      # the printed lower bound of the next singular value of M^2, relative to sigma_1
KERM2_ACCEPT = 4.7e-8    # accepted from here: the rounding limit of the printed bound (the bound taken literally when equal to KERM2_BOUND)
OP = (1.0, 0.05)          # the operating point (N, delta3) of the scan

# Where the printed numbers stand in the supplement (substrings of paper2_supplement.tex as it was audited; a reading aid: each
# check below is the numerical content of the passages listed under its key), in the order of the section.
TEX_QUOTES = {
    "d0_doublet": [r"$\mathrm{Re}\,\Omega\approx\pm5.2\times10^{-4}$", r"$|\mathrm{Im}\,\Omega|\approx1.3\times10^{-3}$"],
    "d0_res_h4": [r"since $\psi_0$ solves the discretized equation only to $\mathcal{O}(h^4)$"],
    "d0_krein": [r"($13/13$ at $\dtres=0$)"],
    "d0_collapse": [r"and they collapse under refinement"],
    "d0_real": [r"$\max|\mathrm{Im}\,\Omega|<10^{-15}$ on the annulus"],
    "d0_gap": [r"$2.3\times10^{-7}$ at $n_\tau=500$"],
    "d0_pair": [r"decays from $7.8\times10^{-3}$ to $8.7\times10^{-4}$ across"],
    "d0_h2": [r"as $h^2$, identifying it as a lattice artifact"],
    "newton": [r"whose Newton residual $\|F(\phi)\|_\infty$ lies below $10^{-12}$"],
    "kerM": [r"below $1.1\times10^{-16}\sigma_1$ and the next exceeds $2\times10^{-4}\sigma_1$", r"so $\dim\ker M=2$"],
    "phase": [r"$8\times10^{-17}$ in relative terms"],
    "translation": [r"($6\times10^{-8}$)"],
    "kerM2": [r"at most $1.4\times10^{-16}\sigma_1$ and the next is at least $4.7\times10^{-8}\sigma_1$", r"so $\dim\ker M^2=4$"],
    "slopes": [r"(log-log slopes $0.48$ and $0.51$ at $\dtres=0.02$ and $0.05$)"],
    "sector": [r"at $1$--$5\times10^{-7}$ with the eigensolvers used", r"they lie below $10^{-6}$"],
    "floor": [r"of order $\sqrt{\epsilon_{\rm mach}\rho}\approx4$--$6\times10^{-7}$ here"],
    "annulus": [r"$|\mathrm{Im}\,\Omega|<10^{-12}$"],
    "scaling": [r"$\rho$ scales as $n_\tau^{2.4}$ to $n_\tau^{2.6}$"],
    "floors": [r"Both floors grow with the grid", r"the bounds quoted here refer to $n_\tau=500$"],
    "q_quartet": [r"the same four modes split into a complex quartet"],
    "q_law": [r"$\max|\mathrm{Im}\,\Omega|\simeq0.68\sqrt{\dtres}$ at $N=1$"],
    "q_value": [r"($0.0957$ at $\dtres=0.02$)"],
    "q_conv": [r"converged to three digits across $n_\tau\in\{200,\dots,600\}$"],
    "q_box": [r"and box-independent"],
    "q_krr": [r"persists on grids that cannot represent the Cherenkov wavenumber $k_{\rm RR}\simeq25$", r"gives $k_{\rm RR}\simeq25$ at $\dtres=0.02$"],
    "q_phi": [r"and vanishes about $\phi$"],
    "scan_ok": [r"the Newton residual is below $10^{-12}$", r"$1.1\times10^{-13}$, and the rank test gives $\dim\ker M=2$ and $\dim\ker M^2=4$"],
    "scan_op": [r"standing tail is $3.6\times10^{-5}$ and the zero sector lies at $3.8\times10^{-7}$"],
    "scan_real": [r"spectrum stays real up to $N\dtres=0.07$", r"(tail $4.3\times10^{-3}$ at $(N,\dtres)=(1.4,0.05)$)"],
    "scan_complex": [r"($\max|\mathrm{Im}\,\Omega|=9.6\times10^{-3}$, $1.3\times10^{-2}$ and", r"$7.6\times10^{-4}$ at $N=1.6$, $1.8$ and $2$, $\dtres=0.05$)",
                     r"tails of $0.5\%$, $1.3\%$ and $8.7\%$"],
    "symbol_max": [r"and it never exceeds $4.61$"],
    "symbol_613": [r"therefore stays above $\mu$ at $\dtres=0.02$ and $n_\tau=500$", r"a zero crossing would require $n_\tau\ge613$"],
    "crossings": [r"the discrete branch crosses zero at $k=10.18$,", r"against $10.10$ in the continuum, and again at $k=42.1$"],
    "census": [r"the annulus holds $34$ eigenvalues, $17$ of each symplectic norm"],
    "neg_energy": [r"Four of them carry negative energy", r"at $\mathrm{Re}\,\Omega=-0.304$ and $-1.348$", r"$+0.304$ and $+1.348$; all $34$ are real to $9.8\times10^{-14}$"],
    "coef": [r"the coefficient is $0.677$ and $0.670$ at $N=1$ for $\dtres=0.02$ and $0.05$, and $0.189$ at",
             r"$N=0.6$ and $3.80$ at $N=2$ for $\dtres=0.02$, against $0.190$ and $3.85$ from the law", r"($3.13$ at $N=2$, $\dtres=0.05$)"],
    "law_2p5": [r"The scan reproduces this law to within $2.5\%$ for"],
}


# ---------------------------------------------------------------------------------------------------------------- setup
def load_package(pkg):
    """Import sothe_p2 from the folder that holds it (the SOTHE-P2-SWEEP package vendors it)."""
    pkg = os.path.abspath(pkg)
    if not os.path.isdir(os.path.join(pkg, "sothe_p2")):
        sys.exit("sothe_p2 not found in %s (set SOTHE_SWEEP_PKG, see the header)" % pkg)
    sys.path.insert(0, pkg)
    from sothe_p2 import kit, suite
    return kit, suite


def build_M(phi, d3, n, mu, D2, D3):
    H = -0.5 * D2 - 2.0 * np.diag(np.abs(phi) ** 2) + mu * np.eye(n)
    T = (-1j * d3) * D3
    return np.block([[H + T, -np.diag(phi ** 2)], [np.diag(np.conj(phi) ** 2), -H + T]])


def annulus(w, mu):
    return (np.abs(w.real) <= 3.0 * mu) & (np.abs(w.real) > 0.1 * mu)


def smallest(w, k=4):
    return w[np.argsort(np.abs(w), kind="stable")[:k]]


def rounds_to(x, printed):
    """True when x rounds to the printed decimal string, that is, lies within half a unit of its last digit (plain or e-notation)."""
    s = printed.strip().lower()
    mant, _, ex = s.partition("e")
    dec = len(mant.split(".")[1]) if "." in mant else 0
    return abs(float(x) - float(s)) <= 0.5 * 10.0 ** (int(ex or 0) - dec) * (1.0 + 1e-9)


# ---------------------------------------------------------------------------------------------------------------- part A
def part_a(Z, suite, kit):
    n, mu = N_TAU, 0.5 * N_SOL ** 2
    h = 2.0 * TAU_MAX / n
    D1, D2, D3 = (suite.fd_circulant(n, h, 4, m) for m in (1, 2, 3))
    rng = np.random.default_rng(SEED)
    out = {}
    for tag, d3 in TAGS:
        phi = Z["phi_" + tag]
        M = build_M(phi, d3, n, mu, D2, D3)
        nM = float(np.linalg.norm(M, 2))
        w = np.linalg.eigvals(M)
        ws = smallest(w)
        rho = float(np.max(np.abs(w)))
        ann = annulus(w, mu)
        F = 0.5 * (D2 @ phi) + np.abs(phi) ** 2 * phi + 1j * d3 * (D3 @ phi) - mu * phi
        wk = kit.bdg_spectrum_bg(phi, N_SOL, d3, n, TAU_MAX)["w"]
        r = {"operator_mismatch": float(np.max(np.abs(np.sort(np.abs(w))[4:] - np.sort(np.abs(wk))[4:])) / rho), "normM": nM, "rho": rho, "sqrt_eps_rho": math.sqrt(EPS * rho), "w_small_here": [[float(x.real), float(x.imag)] for x in ws],
             "max_w_small_here": float(np.max(np.abs(ws))), "max_w_small_altay": float(np.max(np.abs(Z["w_small_phi_" + tag]))),
             "annulus_maxIm_here": float(np.max(np.abs(w[ann].imag))), "annulus_maxIm_altay": float(np.max(np.abs(Z["w_phi_" + tag][annulus(Z["w_phi_" + tag], mu)].imag))),
             "newton_residual_archived_phi": float(np.max(np.abs(F)))}
        rows = []
        for delta in DELTAS:
            mx = []
            for _ in range(3):
                E = rng.standard_normal(M.shape) + 1j * rng.standard_normal(M.shape)
                E /= np.linalg.norm(E, 2)
                wp = np.linalg.eigvals(M + delta * nM * E)
                mx.append(float(np.max(np.sort(np.abs(wp))[:4])))
            rows.append((delta, float(np.median(mx)), min(mx), max(mx)))
        d = np.array([x[0] for x in rows[3:]])
        m = np.array([x[1] for x in rows[3:]])
        r["scaling"] = rows
        r["slope"] = float(np.polyfit(np.log10(d), np.log10(m), 1)[0])
        r["prefactor"] = float(np.median(m / np.sqrt(d * nM)))
        s1 = np.linalg.svd(M, compute_uv=False)
        s2 = np.linalg.svd(M @ M, compute_uv=False)
        r["sv_M_rel"] = (s1[-5:][::-1] / s1[0]).tolist()          # five smallest, smallest first
        r["sv_M2_rel"] = (s2[-6:][::-1] / s2[0]).tolist()         # six smallest, smallest first
        v_ph = np.concatenate([phi, -np.conj(phi)])
        dphi = D1 @ phi
        v_tr = np.concatenate([dphi, np.conj(dphi)])
        res = lambda v: float(np.linalg.norm(M @ v) / (nM * np.linalg.norm(v)))
        r["phase_residual"] = res(v_ph)
        r["translation_residual"] = res(v_tr)
        r["phase_residual_wrong_sign"] = res(np.concatenate([phi, np.conj(phi)]))
        out[tag] = r
    return out


def part_a_delta0(suite):
    """Test 4: delta3 = 0 about psi0: residual of psi0 in the discretized equation, zero sector, gap, annulus, Krein census."""
    mu = 0.5 * N_SOL ** 2
    rows = []
    for n in CONV_N:
        h = 2.0 * TAU_MAX / n
        tau = -TAU_MAX + h * np.arange(n)
        psi0 = N_SOL / np.cosh(N_SOL * tau)
        F = 0.5 * (suite.fd_circulant(n, h, 4, 2) @ psi0) + psi0 ** 3 - mu * psi0
        B = suite.bdg_spectrum(N_SOL, 0.0, n, TAU_MAX)
        ws = smallest(B["w"])
        kr = B["krein"][B["ann"]]
        rows.append({"n_tau": n, "h": h, "im_pair": float(np.max(np.abs(ws.imag))), "re_doublet": float(np.max(np.abs(ws.real))),
                     "annulus_maxIm": B["maxIm"], "gap_rel_dev": abs(B["gap"] - mu) / mu, "n_annulus": int(B["ann"].sum()),
                     "krein_plus": int(np.sum(kr > 0)), "krein_minus": int(np.sum(kr < 0)),
                     "resid_bulk": float(np.max(np.abs(F[np.abs(tau) <= 8.0]))), "resid_outer": float(np.max(np.abs(F[np.abs(tau) > 8.0])))})
    lh = np.log([r["h"] for r in rows])
    ex = lambda key: float(np.polyfit(lh, np.log([r[key] for r in rows]), 1)[0])
    return {"rows": rows, "h_exponent_of_im_pair": ex("im_pair"), "h_exponent_of_re_doublet": ex("re_doublet"), "h_exponent_of_resid_bulk": ex("resid_bulk")}


def k_rr(d3, mu):
    """Phase-matched (Cherenkov) wavenumber of the resonant branch: the largest real root of d3 k^3 - k^2/2 = mu, mu = N^2/2."""
    r = np.roots([d3, -0.5, 0.0, -mu])
    return float(max(x.real for x in r if abs(x.imag) < 1e-9))


def quartet_row(suite, d3, n, tau_max):
    """Test 5: the four smallest eigenvalues (the zero sector) about psi0 = N sech(N tau) at delta3 > 0."""
    B = suite.bdg_spectrum(N_SOL, d3, n, tau_max)
    s = smallest(B["w"])
    h = 2.0 * tau_max / n
    quad = sorted((int(np.sign(z.real)), int(np.sign(z.imag))) for z in s)
    return {"d3": d3, "n_tau": n, "tau_max": tau_max, "h": h, "nyquist": math.pi / h, "im": float(np.max(np.abs(s.imag))),
            "re": float(np.max(np.abs(s.real))), "min_part": float(min(np.min(np.abs(s.real)), np.min(np.abs(s.imag)))),
            "one_per_quadrant": quad == [(-1, -1), (-1, 1), (1, -1), (1, 1)]}


def part_a_quartet(suite):
    """Test 5: the complex quartet about psi0 (delta3 > 0, N = 1): value, law in delta3, convergence in n_tau and in the box."""
    h0 = 2.0 * TAU_MAX / N_TAU
    conv = {tag: [quartet_row(suite, d3, n, TAU_MAX) for n in CONV_N] for tag, d3 in TAGS}
    law = []
    for d3 in Q_DELTAS:
        tags = [t for t, v in TAGS if v == d3]
        law.append([r for r in conv[tags[0]] if r["n_tau"] == N_TAU][0] if tags else quartet_row(suite, d3, N_TAU, TAU_MAX))
    box = [quartet_row(suite, 0.02, int(round(2.0 * tm / h0)), tm) for tm in BOX_TAU]
    return {"law": law, "conv": conv, "box": box, "k_rr": {tag: k_rr(d3, 0.5 * N_SOL ** 2) for tag, d3 in TAGS}}


# ---------------------------------------------------------------------------------------------------------------- part B
def part_b(kit):
    rows = []
    for n in GRID_N:
        for tag, d3 in TAGS:
            t0 = time.time()
            S = kit.tod_soliton(N_SOL, d3, n, TAU_MAX)
            B = kit.bdg_spectrum_bg(S["phi"], N_SOL, d3, n, TAU_MAX)
            small = float(np.max(np.abs(B["w_small"])))
            rho = float(np.max(np.abs(B["w"])))
            rows.append({"n_tau": n, "d3": d3, "small": small, "rho": rho, "sqrt_eps_rho": math.sqrt(EPS * rho), "ratio_small": small / math.sqrt(EPS * rho),
                         "newton_res": float(S["res"][-1]), "converged": bool(S["converged"]), "ratio_newton": float(S["res"][-1]) / (EPS * rho),
                         "annulus_maxIm": B["maxIm"]})
            print("  n_tau %5d d3 %.2f  rho %8.1f  zero sector %.2e (%.2f sqrt(eps rho))  Newton %.1e%s  %.0f s" % (
                n, d3, rho, small, rows[-1]["ratio_small"], rows[-1]["newton_res"], "" if rows[-1]["converged"] else " (not converged)", time.time() - t0), flush=True)
    exps, exps_conv = {}, {}
    for tag, d3 in TAGS:
        sel = [r for r in rows if r["d3"] == d3]
        ln = np.log([r["n_tau"] for r in sel])
        exps[tag] = {"rho": float(np.polyfit(ln, np.log([r["rho"] for r in sel]), 1)[0]),
                     "zero_sector": float(np.polyfit(ln, np.log([r["small"] for r in sel]), 1)[0])}
        cv = [r for r in sel if r["converged"]]
        exps_conv[tag] = float(np.polyfit(np.log([r["n_tau"] for r in cv]), np.log([r["rho"] for r in cv]), 1)[0]) if len(cv) >= 3 else None
    return {"rows": rows, "exponents": exps, "rho_exponent_converged_rows": exps_conv}


# ---------------------------------------------------------------------------------------------------------------- part C
A2 = lambda t: (15.0 - 16.0 * np.cos(t) + np.cos(2.0 * t)) / 6.0                              # k_2eff^2 h^2 of the five-point stencil for -d^2/dtau^2
B3 = lambda t: 0.25 * (13.0 * np.sin(t) - 8.0 * np.sin(2.0 * t) + np.sin(3.0 * t))            # k_3eff h^3 of the seven-point stencil for d^3/dtau^3


def branch(k, d3, mu, h):
    """The discrete positive-norm branch  k_2eff^2/2 + mu - delta3 k_3eff  at the wavenumbers k (theta = k h)."""
    return 0.5 * A2(k * h) / h ** 2 + mu - d3 * B3(k * h) / h ** 3


def branch_min(d3, n, mu, lo=0.05, m=40001):
    """Smallest value of the discrete branch over theta in [lo, pi], and its wavenumber."""
    h = 2.0 * TAU_MAX / n
    th = np.linspace(lo, math.pi, m)
    g = branch(th / h, d3, mu, h)
    i = int(np.argmin(g))
    return float(g[i]), float(th[i] / h)


def branch_ratio_max(d3, n, m=40001):
    """Largest ratio delta3 k_3eff / (k_2eff^2 / 2) of the discrete branch over theta in (0, pi], and its wavenumber. The branch minus mu is
    (k_2eff^2 / 2) (1 - ratio), so it stays above mu exactly where the ratio is below 1."""
    h = 2.0 * TAU_MAX / n
    th = np.linspace(1e-3, math.pi, m)
    r = 2.0 * (d3 / h) * B3(th) / A2(th)
    i = int(np.argmax(r))
    return float(r[i]), float(th[i] / h)


def zero_crossings(d3, n, mu):
    """Zeros of the discrete branch on 0 < k <= pi/h (sign changes on a fine grid, refined by bisection)."""
    h = 2.0 * TAU_MAX / n
    k = np.linspace(1e-3, math.pi / h, 400001)
    g = branch(k, d3, mu, h)
    roots = []
    for i in np.nonzero(np.sign(g[:-1]) != np.sign(g[1:]))[0]:
        a, b = float(k[i]), float(k[i + 1])
        fa = float(g[i])
        for _ in range(80):
            c = 0.5 * (a + b)
            fc = float(branch(c, d3, mu, h))
            if (fc > 0) == (fa > 0):
                a, fa = c, fc
            else:
                b = c
        roots.append(0.5 * (a + b))
    return roots


def load_b5(run):
    """Block B5 of the archived sweep run: the scan table and the spectrum at the operating point (N = 1, delta3 = 0.05)."""
    pts = os.path.join(run, "B5_bdg", "bdg_points.csv")
    spec = os.path.join(run, "B5_bdg", "spectra", "bdg_N1_d0.05.csv")
    for p in (pts, spec):
        if not os.path.isfile(p):
            sys.exit("not found: %s (SOTHE_SWEEP_RUN must be the run folder 20260928T131100Z of the sweep results, see the header)" % p)
    keep = ("newton_converged", "newton_res", "tail_rel", "annulus_maxIm", "n_annulus", "krein_pos", "krein_neg", "zero_sector_max", "rho",
            "dim_ker_M", "dim_ker_M2", "psi0_annulus_maxIm", "psi0_annulus_over_sqrt_d3")
    table = []
    with open(pts, newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            r = {"N": round(float(row["N"]), 6), "d3": round(float(row["d3"]), 6)}
            for k in keep:
                v = float(row[k])
                r[k] = None if math.isnan(v) else v
            table.append(r)
    with open(spec, newline="", encoding="utf-8") as f:
        sp = [[float(r["Re_omega"]), float(r["Im_omega"]), int(round(float(r["krein"])))] for r in csv.DictReader(f)]
    return table, sp


def part_c(suite, sweep_run):
    """Part C: the scan, the census and the dispersion argument of the last three paragraphs of the section."""
    table, sp = load_b5(sweep_run)
    n, h = N_TAU, 2.0 * TAU_MAX / N_TAU
    th = 2.0 * np.pi * np.arange(n) / n
    D2, D3 = (suite.fd_circulant(n, h, 4, m) for m in (2, 3))
    mis2 = float(np.max(np.abs(np.fft.fft(D2[:, 0]) + A2(th) / h ** 2))) * h ** 2       # the symbol of D2 is -k_2eff^2
    mis3 = float(np.max(np.abs(np.fft.fft(D3[:, 0]) + 1j * B3(th) / h ** 3))) * h ** 3  # the symbol of D3 is -i k_3eff
    tg = np.linspace(0.0, math.pi, 2000001)
    b3 = B3(tg)
    i = int(np.argmax(b3))
    mu = 0.5 * N_SOL ** 2
    ratio, k_ratio = branch_ratio_max(0.02, N_TAU)
    ratio_612 = branch_ratio_max(0.02, 612)[0]
    first_n, gmin_at = None, {}
    for nn in range(100, 1000):
        gmin_at[nn] = branch_min(0.02, nn, mu)[0]
        if gmin_at[nn] <= 0.0 and first_n is None:
            first_n = nn
    return {"table": table, "spectrum": sp,
            "symbol": {"b3_max": float(b3[i]), "theta_of_b3_max": float(tg[i]), "k_of_b3_max": float(tg[i] / h), "mismatch_D2": mis2, "mismatch_D3": mis3},
            "branch": {"d3_0p02_n500_largest_ratio": ratio, "k_of_that_ratio": k_ratio, "d3_0p02_n612_largest_ratio": ratio_612, "first_n_with_zero_crossing": first_n,
                       "min_g_at_first_n": gmin_at.get(first_n), "min_g_one_below": gmin_at.get((first_n or 0) - 1), "crossings_d3_0p05_n500": zero_crossings(0.05, N_TAU, mu),
                       "continuum_root_d3_0p05": k_rr(0.05, mu)}}


# ---------------------------------------------------------------------------------------------------------------- statements
def sci(x):
    return "%.2e" % x


def spread(rows):
    """Relative spread (max - min) / mean of max|Im| over a list of quartet rows."""
    v = [r["im"] for r in rows]
    return (max(v) - min(v)) / float(np.mean(v))


def evaluate(A, D0, Q, B, C):
    """The numerical statements of the supplement, in the order of the section, each with the computed values; the first entry is the
    check of the operator. Returns [(id, statement, holds, detail, literal)]: literal is False where a check holds only within the
    rounding of its printed value or bound, or only with the factor allowed at the round-off level."""
    both = [A[t] for t, _ in TAGS]
    out = []

    def add(key, text, ok, detail, literal=True):
        out.append((key, text, bool(ok), detail, bool(literal)))

    add("operator", "the operator M built here has the spectrum of sothe_p2.kit.bdg_spectrum_bg about the same phi (the moduli outside the zero sector agree to 1e-10 rho)",
        all(r["operator_mismatch"] < 1e-10 for r in both), "%s and %s (relative to rho)" % tuple(sci(r["operator_mismatch"]) for r in both))
    R = {(r["N"], r["d3"]): r for r in C["table"]}
    nd = lambda r: r["N"] * r["d3"]
    rows0 = D0["rows"]
    ref = [r for r in rows0 if r["n_tau"] == 500][0]
    z1 = R[(1.0, 0.0)]
    # -- the first paragraph: the lattice artifact of the zero sector about psi0
    add("d0_doublet", "delta3 = 0, n_tau = 500: real doublet at +-5.2e-4 and imaginary pair at 1.3e-3 (rounded to the printed digits)",
        rounds_to(ref["re_doublet"], "5.2e-4") and rounds_to(ref["im_pair"], "1.3e-3"), "%s and %s" % (sci(ref["re_doublet"]), sci(ref["im_pair"])))
    add("d0_res_h4", "delta3 = 0: psi0 = N sech(N tau) solves the discretized equation only to O(h^4): its residual in the bulk |tau| <= 8, n_tau = 200 ... 600, scales as h^4 (fitted exponent within 0.1 of 4)",
        abs(D0["h_exponent_of_resid_bulk"] - 4.0) <= 0.1,
        "exponent %.2f; residual/h^4 from %.3f to %.3f (outer region |tau| > 8: %s at n_tau = 200, %s at 600, the periodic wrap)" % (
            D0["h_exponent_of_resid_bulk"], rows0[0]["resid_bulk"] / rows0[0]["h"] ** 4, rows0[-1]["resid_bulk"] / rows0[-1]["h"] ** 4,
            sci(rows0[0]["resid_outer"]), sci(rows0[-1]["resid_outer"])))
    add("d0_krein", "delta3 = 0: the annulus census is balanced, 13/13 (n_tau = 500 here, and N = 1 of the archived sweep run)",
        (ref["krein_plus"], ref["krein_minus"]) == (13, 13) and (z1["krein_pos"], z1["krein_neg"]) == (13, 13),
        "%d/%d here; %d/%d in the archive" % (ref["krein_plus"], ref["krein_minus"], z1["krein_pos"], z1["krein_neg"]))
    dec = lambda key: all(rows0[i + 1][key] < rows0[i][key] for i in range(len(rows0) - 1))
    add("d0_collapse", "delta3 = 0: the four zero-sector eigenvalues collapse under refinement (the real doublet and the imaginary pair decrease with n_tau = 200 ... 600, both as h^2 within 0.1)",
        dec("re_doublet") and dec("im_pair") and abs(D0["h_exponent_of_re_doublet"] - 2.0) <= 0.1 and abs(D0["h_exponent_of_im_pair"] - 2.0) <= 0.1,
        "real doublet %s to %s (h^%.2f), imaginary pair %s to %s (h^%.2f)" % (sci(rows0[0]["re_doublet"]), sci(rows0[-1]["re_doublet"]), D0["h_exponent_of_re_doublet"],
                                                                        sci(rows0[0]["im_pair"]), sci(rows0[-1]["im_pair"]), D0["h_exponent_of_im_pair"]))
    # -- the convergence paragraph: delta3 = 0
    add("d0_real", "delta3 = 0: annulus spectrum real, max|Im| below 1e-15 (round-off level: factor 2 accepted; n_tau = 500 here, and N = 1 of the archived sweep run)",
        ref["annulus_maxIm"] < 2e-15 and z1["annulus_maxIm"] < 2e-15, "%s here; %s in the archive" % (sci(ref["annulus_maxIm"]), sci(z1["annulus_maxIm"])),
        literal=ref["annulus_maxIm"] < 1e-15 and z1["annulus_maxIm"] < 1e-15)
    add("d0_gap", "delta3 = 0: continuum edge reproduces the gap mu to a relative 2.3e-7 at n_tau = 500 (rounded to the printed digits)", rounds_to(ref["gap_rel_dev"], "2.3e-7"), sci(ref["gap_rel_dev"]))
    a, b = rows0[0], rows0[-1]
    add("d0_pair", "delta3 = 0: imaginary pair of the zero sector decays from 7.8e-3 to 8.7e-4 across n_tau in {200, 600} (rounded to the printed digits)",
        rounds_to(a["im_pair"], "7.8e-3") and rounds_to(b["im_pair"], "8.7e-4"), "%s at %d, %s at %d" % (sci(a["im_pair"]), a["n_tau"], sci(b["im_pair"]), b["n_tau"]))
    add("d0_h2", "delta3 = 0: the decay is as h^2", abs(D0["h_exponent_of_im_pair"] - 2.0) < 0.1, "fitted exponent %.2f over n_tau = 200 ... 600" % D0["h_exponent_of_im_pair"])
    # -- the convergence paragraph: about phi
    add("newton", "Newton residual of phi below 1e-12 (the phi of the archived run)", all(r["newton_residual_archived_phi"] < 1e-12 for r in both),
        "%s, %s" % tuple(sci(r["newton_residual_archived_phi"]) for r in both))
    s_small = max(max(r["sv_M_rel"][0], r["sv_M_rel"][1]) for r in both)
    s_next = min(r["sv_M_rel"][2] for r in both)
    add("kerM", "two smallest singular values of M below 1.1e-16 sigma_1 (round-off level: factor 2 accepted), the next above 2e-4 sigma_1",
        s_small < 2 * 1.1e-16 and s_next > 2e-4, "largest of the two smallest %s; next %s (min over delta3)" % (sci(s_small), sci(s_next)),
        literal=s_small < 1.1e-16 and s_next > 2e-4)
    ph = max(r["phase_residual"] for r in both)
    add("phase", "phase mode (phi, -conj(phi)) is a null vector to 8e-17 (round-off level: below 2e-16)", ph < 2e-16, "%s; the other sign convention gives %s" % (
        sci(ph), sci(min(r["phase_residual_wrong_sign"] for r in both))), literal=ph < 8.5e-17)
    tr = max(r["translation_residual"] for r in both)
    add("translation", "translation mode (phi', conj(phi')) is a null vector to the O(h^4) accuracy of the difference operator, 6e-8 (rounded to the printed digit; the larger of the two values of delta3)",
        rounds_to(tr, "6e-8"), "%s and %s" % tuple(sci(r["translation_residual"]) for r in both))
    s2_small = max(max(r["sv_M2_rel"][:4]) for r in both)
    s2_next = min(r["sv_M2_rel"][4] for r in both)
    kerm2_text = "four smallest singular values of M^2 at most 1.4e-16 sigma_1 (round-off level: factor 2), the next at least %.1e sigma_1" % KERM2_BOUND
    if KERM2_ACCEPT != KERM2_BOUND:
        kerm2_text += " (accepted from %.2e, the rounding limit of the printed value)" % KERM2_ACCEPT
    add("kerM2", kerm2_text, s2_small < 2 * 1.4e-16 and s2_next >= KERM2_ACCEPT,
        "largest of the four %s; next %s (min over delta3; %s at delta3 = 0.02, %s at 0.05)%s" % (
            sci(s2_small), sci(s2_next), sci(both[0]["sv_M2_rel"][4]), sci(both[1]["sv_M2_rel"][4]),
            "" if s2_next >= KERM2_BOUND else "; the printed bound %.1e is not met literally, its rounding limit %.2e is" % (KERM2_BOUND, KERM2_ACCEPT)),
        literal=s2_small <= 1.4e-16 and s2_next >= KERM2_BOUND)
    sl = [r["slope"] for r in both]
    add("slopes", "random perturbations split the zero sector as delta^(1/2): log-log slopes 0.48 and 0.51 (accepted: within 0.01 of the printed values)",
        abs(sl[0] - 0.48) <= 0.01 and abs(sl[1] - 0.51) <= 0.01, "%.3f and %.3f (a semisimple zero would give 1)" % tuple(sl))
    sec_here = [r["max_w_small_here"] for r in both]
    sec_altay = [r["max_w_small_altay"] for r in both]
    add("sector", "the four eigenvalues of the sector lie at 1-5e-7 (the archived Altay values) and below 1e-6 (the archived values and this run, whose values depend on the BLAS library)",
        all(1e-7 <= x <= 5e-7 for x in sec_altay) and max(sec_altay + sec_here) < 1e-6,
        "Altay %s, %s; this run %s, %s" % (sci(sec_altay[0]), sci(sec_altay[1]), sci(sec_here[0]), sci(sec_here[1])))
    fl = [r["sqrt_eps_rho"] for r in both]
    add("floor", "sqrt(eps_mach rho) is about 4-6e-7", all(4e-7 <= x <= 6e-7 for x in fl), "%s and %s (rho %.1f and %.1f)" % (sci(fl[0]), sci(fl[1]), both[0]["rho"], both[1]["rho"]))
    an = [r["annulus_maxIm_here"] for r in both] + [r["annulus_maxIm_altay"] for r in both]
    add("annulus", "annulus spectrum real to round-off, |Im| below 1e-12", max(an) < 1e-12, "largest %s (this run and Altay, both delta3)" % sci(max(an)))
    if B is not None:
        e = [B["exponents"][t]["rho"] for t, _ in TAGS]
        ec = [B["rho_exponent_converged_rows"][t] for t, _ in TAGS]
        add("scaling", "the spectral radius rho scales as n_tau^2.4 to n_tau^2.6 (delta3 = 0.02 and 0.05; fit over n_tau = 250 ... 1250, rounded to the printed digits; the fit over the rows whose Newton solve converged gives the same)",
            rounds_to(e[0], "2.4") and rounds_to(e[1], "2.6") and all(x is not None and rounds_to(x, p) for x, p in zip(ec, ("2.4", "2.6"))),
            "exponents %.2f and %.2f (converged rows only: %s)" % (e[0], e[1], " and ".join("%.2f" % x if x is not None else "n/a" for x in ec)))
        det, okf = [], True
        for tag, d3 in TAGS:
            rows = sorted([r for r in B["rows"] if r["d3"] == d3], key=lambda r: r["n_tau"])
            r500 = [r for r in rows if r["n_tau"] == N_TAU][0]
            nc = ["%d (residual %s)" % (r["n_tau"], sci(r["newton_res"])) for r in rows if not r["converged"]]
            okf = okf and r500["converged"] and B["exponents"][tag]["zero_sector"] > 0 and rows[-1]["annulus_maxIm"] > rows[0]["annulus_maxIm"] and r500["small"] < 1e-6 and r500["annulus_maxIm"] < 1e-12
            det.append("delta3 = %.2f: sector ~ n_tau^%.2f, annulus max|Im| %s to %s, at n_tau = 500 %s and %s; Newton not converged at: %s" % (
                d3, B["exponents"][tag]["zero_sector"], sci(rows[0]["annulus_maxIm"]), sci(rows[-1]["annulus_maxIm"]), sci(r500["small"]), sci(r500["annulus_maxIm"]), ", ".join(nc) or "none"))
        add("floors", "both floors grow with the grid (n_tau = 250 ... 1250: the sector maximum with a positive fitted exponent, the annulus max|Im| larger at the end than at the start), and the bounds 1e-6 and 1e-12 hold at n_tau = 500 (where the Newton solve converged)",
            okf, "; ".join(det))
    # -- the quartet about psi0
    r02 = [r for r in Q["law"] if r["d3"] == 0.02][0]
    add("q_quartet", "about psi0 (N = 1, n_tau = 500) the four zero-sector modes form a complex quartet: one eigenvalue in each quadrant, real and imaginary parts above 1e-3, for delta3 in {0.005, 0.01, 0.02, 0.05}",
        all(r["one_per_quadrant"] and r["min_part"] > 1e-3 for r in Q["law"]), "delta3 = 0.02: +-%.5f +- %.5f i" % (r02["re"], r02["im"]))
    coef = [r["im"] / math.sqrt(r["d3"]) for r in Q["law"]]
    add("q_law", "max|Im Omega| = 0.68 sqrt(delta3) at N = 1 for delta3 in {0.005, 0.01, 0.02, 0.05} (a law with ~: coefficient within 0.015 of 0.68)",
        all(abs(c - 0.68) <= 0.015 for c in coef), "coefficients " + ", ".join("%.4f" % c for c in coef))
    add("q_value", "the quartet has max|Im Omega| = 0.0957 at delta3 = 0.02 (rounded to the printed digits)", rounds_to(r02["im"], "0.0957"), "%.5f" % r02["im"])
    sp = [spread(Q["conv"][t]) for t, _ in TAGS]
    add("q_conv", "the quartet is converged to three digits across n_tau in {200, ..., 600} (relative spread of max|Im| below 1e-3; delta3 = 0.02 and 0.05)",
        all(x < 1e-3 for x in sp), "relative spreads %s and %s" % tuple(sci(x) for x in sp))
    add("q_box", "the quartet is box-independent: tau_max in {8, 12, 16, 20, 24} at the grid spacing of n_tau = 500 (delta3 = 0.02) changes max|Im| by less than a relative 1e-3",
        spread(Q["box"]) < 1e-3, "relative spread %s" % sci(spread(Q["box"])))
    kr = Q["k_rr"]["d3_0p02"]
    r200 = [r for r in Q["conv"]["d3_0p02"] if r["n_tau"] == 200][0]
    add("q_krr", "delta3 k^3 - k^2/2 = N^2/2 gives k_RR = 25 at delta3 = 0.02 (a value with ~: within 0.5), and the grid n_tau = 200, whose Nyquist wavenumber pi/h is below k_RR, still has the quartet (max|Im| within a relative 1e-3 of its value at n_tau = 500)",
        abs(kr - 25.0) < 0.5 and r200["nyquist"] < kr and abs(r200["im"] / r02["im"] - 1.0) < 1e-3,
        "k_RR = %.2f (%.2f at delta3 = 0.05); pi/h = %.1f at n_tau = 200; max|Im| %.5f against %.5f" % (kr, Q["k_rr"]["d3_0p05"], r200["nyquist"], r200["im"], r02["im"]))
    im_phi = max(abs(z[1]) for t, _ in TAGS for z in A[t]["w_small_here"])
    add("q_phi", "about the stationary soliton phi the four sector eigenvalues have |Im| below 1e-6 (this run): the quartet is absent", im_phi < 1e-6,
        "largest %s, against %.4f about psi0 at delta3 = 0.02" % (sci(im_phi), r02["im"]))
    # -- the scan (block B5 of the archived sweep run)
    inr = [r for r in C["table"] if nd(r) <= 0.05 + 1e-9]
    res_max = max(r["newton_res"] for r in inr)
    im_max = max(r["annulus_maxIm"] for r in inr)
    add("scan_ok", "scan over N in {0.6, ..., 2} and delta3 in {0, 0.02, 0.05} (archived sweep run): at every point with N delta3 <= 0.05, the operating point included, the Newton solve has converged with a residual below 1e-12, "
        "the annulus spectrum is real to 1.1e-13, and dim ker M = 2 and dim ker M^2 = 4",
        len(C["table"]) == 24 and OP in R and nd(R[OP]) <= 0.05 + 1e-9 and all(r["newton_converged"] == 1 for r in inr) and res_max < 1e-12 and im_max <= 1.1e-13
        and all(r["dim_ker_M"] == 2 and r["dim_ker_M2"] == 4 for r in inr),
        "%d of %d points; largest Newton residual %s, largest annulus max|Im| %s" % (len(inr), len(C["table"]), sci(res_max), sci(im_max)))
    op = R[OP]
    rec = None
    if B is not None:
        rec = [r for r in B["rows"] if r["n_tau"] == N_TAU and r["d3"] == OP[1]][0]
    add("scan_op", "at the operating point (N = 1, delta3 = 0.05) the standing tail is 3.6e-5 and the zero sector lies at 3.8e-7 (rounded to the printed digits)%s" % (
        "; the recomputation of Part B (n_tau = 500) has the same spectral radius to a relative 1e-6" if rec else ""),
        rounds_to(op["tail_rel"], "3.6e-5") and rounds_to(op["zero_sector_max"], "3.8e-7") and (rec is None or abs(rec["rho"] / op["rho"] - 1.0) < 1e-6),
        "tail %s, zero sector %s (archive)%s" % (sci(op["tail_rel"]), sci(op["zero_sector_max"]), "; rho %.4f here, %.4f in the archive" % (rec["rho"], op["rho"]) if rec else ""))
    lim = R[(1.4, 0.05)]
    mid = sorted(nd(r) for r in C["table"] if 0.05 + 1e-9 < nd(r) <= 0.07 + 1e-9)
    add("scan_real", "beyond N delta3 = 0.05 the spectrum stays real up to N delta3 = 0.07 (annulus max|Im| below 1e-12 at every point with N delta3 <= 0.07), with a tail of 4.3e-3 at (N, delta3) = (1.4, 0.05) (rounded to the printed digits)",
        all(r["annulus_maxIm"] < 1e-12 for r in C["table"] if nd(r) <= 0.07 + 1e-9) and rounds_to(nd(lim), "0.07") and rounds_to(lim["tail_rel"], "4.3e-3"),
        "N delta3 = %s: largest max|Im| %s; tail %s at (1.4, 0.05)" % (", ".join("%.2f" % x for x in mid), sci(max(r["annulus_maxIm"] for r in C["table"] if nd(r) <= 0.07 + 1e-9)), sci(lim["tail_rel"])))
    cx = [R[(nn, 0.05)] for nn in (1.6, 1.8, 2.0)]
    add("scan_complex", "the spectrum acquires complex eigenvalues at N delta3 = 0.08, 0.09 and 0.10 (max|Im| = 9.6e-3, 1.3e-2 and 7.6e-4 at N = 1.6, 1.8 and 2, delta3 = 0.05), where the tails are 0.5 %, 1.3 % and 8.7 % (rounded to the printed digits)",
        all(rounds_to(nd(r), p) for r, p in zip(cx, ("0.08", "0.09", "0.10"))) and all(rounds_to(r["annulus_maxIm"], p) for r, p in zip(cx, ("9.6e-3", "1.3e-2", "7.6e-4")))
        and all(r["annulus_maxIm"] > 1e-6 for r in cx) and all(rounds_to(100.0 * r["tail_rel"], p) for r, p in zip(cx, ("0.5", "1.3", "8.7"))),
        "max|Im| %s; tails %s %%" % (", ".join(sci(r["annulus_maxIm"]) for r in cx), ", ".join("%.2f" % (100.0 * r["tail_rel"]) for r in cx)))
    # -- the dispersion of the stencils and the census
    S, Bc = C["symbol"], C["branch"]
    add("symbol_max", "the symbols of the two stencils are those of the matrices of sothe_p2 (n_tau = 500, agreement to 1e-10), and k_3eff h^3 = (13 sin t - 8 sin 2t + sin 3t)/4 never exceeds 4.61 (its maximum, rounded to the printed digits)",
        S["mismatch_D2"] < 1e-10 and S["mismatch_D3"] < 1e-10 and S["b3_max"] <= 4.61 and rounds_to(S["b3_max"], "4.61"),
        "maximum %.6f at theta = %.4f (k = %.1f at n_tau = 500); the symbols of D2 and D3 differ from the printed formulas by %s and %s" % (
            S["b3_max"], S["theta_of_b3_max"], S["k_of_b3_max"], sci(S["mismatch_D2"]), sci(S["mismatch_D3"])))
    add("symbol_613", "at delta3 = 0.02 and n_tau = 500 the discrete positive-norm branch k_2eff^2/2 + mu - delta3 k_3eff stays above mu (the ratio delta3 k_3eff / (k_2eff^2/2) is below 1 for all k), and a zero crossing first occurs at n_tau = 613",
        Bc["d3_0p02_n500_largest_ratio"] < 1.0 and Bc["first_n_with_zero_crossing"] == 613,
        "largest ratio %.3f at k = %.1f (%.3f at n_tau = 612); first zero crossing at n_tau = %s (smallest branch value %.3f, against %.3f at n_tau = %s)" % (
            Bc["d3_0p02_n500_largest_ratio"], Bc["k_of_that_ratio"], Bc["d3_0p02_n612_largest_ratio"], Bc["first_n_with_zero_crossing"], Bc["min_g_at_first_n"], Bc["min_g_one_below"],
            (Bc["first_n_with_zero_crossing"] or 0) - 1))
    rts = Bc["crossings_d3_0p05_n500"]
    add("crossings", "at delta3 = 0.05 (n_tau = 500) the discrete branch crosses zero at k = 10.18 (10.10 in the continuum) and again at k = 42.1, beyond the maximum of the stencil symbol (rounded to the printed digits)",
        len(rts) == 2 and rounds_to(rts[0], "10.18") and rounds_to(Bc["continuum_root_d3_0p05"], "10.10") and rounds_to(rts[1], "42.1") and rts[1] > S["k_of_b3_max"],
        "crossings %s; continuum %.3f; the symbol turns over at k = %.1f" % (", ".join("%.3f" % x for x in rts), Bc["continuum_root_d3_0p05"], S["k_of_b3_max"]))
    w_op = np.array([x[0] + 1j * x[1] for x in C["spectrum"]])
    kr_op = np.array([x[2] for x in C["spectrum"]])
    ann = annulus(w_op, 0.5 * OP[0] ** 2)
    npos, nneg = int(np.sum(kr_op[ann] > 0)), int(np.sum(kr_op[ann] < 0))
    add("census", "at the operating point the annulus holds 34 eigenvalues, 17 of each symplectic norm (the spectrum file and the table of the scan agree)",
        int(ann.sum()) == 34 and (npos, nneg) == (17, 17) and (op["n_annulus"], op["krein_pos"], op["krein_neg"]) == (34, 17, 17),
        "%d eigenvalues, %d and %d; %d in the spectrum file altogether (four in the zero sector)" % (int(ann.sum()), npos, nneg, len(w_op)))
    neg = sorted((round(float(x.real), 3), int(k)) for x, k in zip(w_op[ann], kr_op[ann]) if x.real * k < 0)
    im_ann = float(np.max(np.abs(w_op[ann].imag)))
    add("neg_energy", "four annulus modes carry negative energy (Re Omega times the symplectic sign is below 0): the positive-norm modes at Re Omega = -0.304 and -1.348 and their mirror images at +0.304 and +1.348; "
        "all 34 annulus eigenvalues are real to 9.8e-14 (rounded to the printed digits)",
        neg == [(-1.348, 1), (-0.304, 1), (0.304, -1), (1.348, -1)] and rounds_to(im_ann, "9.8e-14"),
        "modes %s; largest |Im| %s" % (", ".join("%+.3f (%+d)" % x for x in neg), sci(im_ann)))
    law = lambda nn: 0.68 * nn ** 2.5
    cf = lambda nn, dd: R[(nn, dd)]["psi0_annulus_over_sqrt_d3"]
    add("coef", "coefficient max|Im Omega| / sqrt(delta3) of the psi0 quartet in the scan: 0.677 and 0.670 at N = 1 for delta3 = 0.02 and 0.05, 0.189 at N = 0.6 and 3.80 at N = 2 for delta3 = 0.02, against 0.190 and 3.85 from the law 0.68 N^(5/2), "
        "and 3.13 at N = 2, delta3 = 0.05 (rounded to the printed digits)",
        rounds_to(cf(1.0, 0.02), "0.677") and rounds_to(cf(1.0, 0.05), "0.670") and rounds_to(cf(0.6, 0.02), "0.189") and rounds_to(cf(2.0, 0.02), "3.80")
        and rounds_to(law(0.6), "0.190") and rounds_to(law(2.0), "3.85") and rounds_to(cf(2.0, 0.05), "3.13"),
        "%.4f, %.4f, %.4f, %.4f; law %.4f and %.4f; %.4f" % (cf(1.0, 0.02), cf(1.0, 0.05), cf(0.6, 0.02), cf(2.0, 0.02), law(0.6), law(2.0), cf(2.0, 0.05)))
    dev = [(nd(r), r["psi0_annulus_over_sqrt_d3"] / law(r["N"]) - 1.0, r["N"], r["d3"]) for r in C["table"] if r["d3"] > 0]
    inside = [x for x in dev if x[0] <= 0.07 + 1e-9]
    outside = sorted(x for x in dev if x[0] > 0.07 + 1e-9)
    worst = max(inside, key=lambda x: abs(x[1]))
    add("law_2p5", "the scan reproduces the law 0.68 N^(5/2) sqrt(delta3) to within 2.5 % for N delta3 <= 0.07 (delta3 > 0), and falls below it, increasingly, as N delta3 approaches 0.1",
        abs(worst[1]) <= 0.025 and all(x[1] < 0 for x in outside) and all(outside[i + 1][1] < outside[i][1] for i in range(len(outside) - 1)),
        "largest deviation %.2f %% at (N, delta3) = (%g, %g) among %d points; beyond: %s" % (
            100.0 * worst[1], worst[2], worst[3], len(inside), ", ".join("%.1f %% at N delta3 = %.2f" % (100.0 * x[1], x[0]) for x in outside)))
    return out


# ---------------------------------------------------------------------------------------------------------------- report
def report(A, D0, Q, B, C):
    for tag, d3 in TAGS:
        r = A[tag]
        print("\n== delta3 = %.2f, n_tau = %d   ||M||_2 = %.1f   rho = %.1f   sqrt(eps rho) = %.2e" % (d3, N_TAU, r["normM"], r["rho"], r["sqrt_eps_rho"]))
        print("   zero sector (this run):", ", ".join("%.2e%+.2ei" % (z[0], z[1]) for z in r["w_small_here"]), "  max|w| %.2e; archived Altay run: %.2e" % (r["max_w_small_here"], r["max_w_small_altay"]))
        for x in r["scaling"]:
            print("   delta %.0e  max|w_small| median %.2e (range %.2e .. %.2e)" % x)
        print("   log-log slope over delta in [1e-13, 1e-8]: %.3f (0.5 = 2x2 Jordan blocks, 1 = semisimple); |w| ~ %.2f sqrt(delta ||M||)" % (r["slope"], r["prefactor"]))
        print("   smallest singular values of M   / sigma_1:", " ".join("%.2e" % x for x in r["sv_M_rel"]))
        print("   smallest singular values of M^2 / sigma_1:", " ".join("%.2e" % x for x in r["sv_M2_rel"]))
        print("   phase mode residual %.1e (wrong sign %.1e), translation mode residual %.1e" % (r["phase_residual"], r["phase_residual_wrong_sign"], r["translation_residual"]))
    print("\n== delta3 = 0 about psi0")
    for x in D0["rows"]:
        print("   n_tau %4d  imaginary pair %.3e  real doublet %.3e  annulus max|Im| %.1e  gap deviation %.2e  Krein %d/%d  residual of psi0 (|tau| <= 8) %.3e" % (
            x["n_tau"], x["im_pair"], x["re_doublet"], x["annulus_maxIm"], x["gap_rel_dev"], x["krein_plus"], x["krein_minus"], x["resid_bulk"]))
    print("   h-exponents: imaginary pair %.2f, real doublet %.2f, residual of psi0 %.2f" % (D0["h_exponent_of_im_pair"], D0["h_exponent_of_re_doublet"], D0["h_exponent_of_resid_bulk"]))
    print("\n== about psi0 at delta3 > 0, N = 1: the complex quartet of the zero sector")
    for r in Q["law"]:
        print("   delta3 %.3f  n_tau %d: max|Im| %.5f = %.4f sqrt(delta3), |Re| %.5f" % (r["d3"], r["n_tau"], r["im"], r["im"] / math.sqrt(r["d3"]), r["re"]))
    for tag, d3 in TAGS:
        print("   delta3 %.2f, n_tau = %s: max|Im| %s (relative spread %.1e)" % (d3, ", ".join(str(r["n_tau"]) for r in Q["conv"][tag]), ", ".join("%.5f" % r["im"] for r in Q["conv"][tag]), spread(Q["conv"][tag])))
    print("   delta3 0.02, tau_max = %s: max|Im| %s (relative spread %.1e)" % (", ".join("%g" % r["tau_max"] for r in Q["box"]), ", ".join("%.5f" % r["im"] for r in Q["box"]), spread(Q["box"])))
    print("   phase matching: k_RR = %.2f at delta3 = 0.02 and %.2f at delta3 = 0.05; Nyquist wavenumber pi/h = %.1f at n_tau = 200" % (Q["k_rr"]["d3_0p02"], Q["k_rr"]["d3_0p05"], Q["conv"]["d3_0p02"][0]["nyquist"]))
    if B is not None:
        print("\n== grid scan (Part B)")
        for tag, d3 in TAGS:
            e = B["exponents"][tag]
            print("   delta3 = %.2f: rho ~ n_tau^%.2f, zero sector ~ n_tau^%.2f" % (d3, e["rho"], e["zero_sector"]))
    S, Bc = C["symbol"], C["branch"]
    print("\n== scan, census and dispersion (Part C): block B5 of the archived sweep run, and closed forms")
    op = [r for r in C["table"] if (r["N"], r["d3"]) == OP][0]
    print("   scan table: %d points; operating point (N, delta3) = (1, 0.05): tail %.3e, zero sector %.3e, rho %.4f" % (len(C["table"]), op["tail_rel"], op["zero_sector_max"], op["rho"]))
    print("   stencil symbols: max of k_3eff h^3 = %.6f at theta = %.4f; the symbols of the matrices differ from the formulas by %.1e (D2) and %.1e (D3)" % (
        S["b3_max"], S["theta_of_b3_max"], S["mismatch_D2"], S["mismatch_D3"]))
    print("   discrete branch, delta3 = 0.02: first zero crossing at n_tau = %s; delta3 = 0.05, n_tau = 500: zeros at k = %s (continuum %.3f)" % (
        Bc["first_n_with_zero_crossing"], ", ".join("%.3f" % x for x in Bc["crossings_d3_0p05_n500"]), Bc["continuum_root_d3_0p05"]))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--grid", action="store_true", help="also run Part B, the grid scan")
    ap.add_argument("--sweep-pkg", default=os.environ.get("SOTHE_SWEEP_PKG") or os.path.join(HERE, "..", "SOTHE-P2-SWEEP"))
    ap.add_argument("--p2-run", default=os.environ.get("SOTHE_P2_RUN"), help="unpacked run folder 20260927T175224Z")
    ap.add_argument("--sweep-run", default=os.environ.get("SOTHE_SWEEP_RUN"), help="unpacked run folder 20260928T131100Z of the sweep results")
    ap.add_argument("--out", default=os.environ.get("SOTHE_ARTICLE_OUT") or os.path.join(os.getcwd(), "article_output"))
    a = ap.parse_args()
    if not a.p2_run:
        sys.exit("set SOTHE_P2_RUN or --p2-run (see the header)")
    if not a.sweep_run:
        sys.exit("set SOTHE_SWEEP_RUN or --sweep-run (see the header)")
    npz = os.path.join(a.p2_run, "S1_g4a", "bdg_backgrounds.npz")
    if not os.path.isfile(npz):
        sys.exit("not found: %s" % npz)
    kit, suite = load_package(a.sweep_pkg)
    Z = np.load(npz, allow_pickle=False)
    t0 = time.time()
    print("Part A: n_tau = %d, delta3 in {0.02, 0.05}, phi of the archived run (%s)" % (N_TAU, npz))
    A = part_a(Z, suite, kit)
    D0 = part_a_delta0(suite)
    Q = part_a_quartet(suite)
    B = None
    if a.grid:
        print("\nPart B: grid scan")
        B = part_b(kit)
    C = part_c(suite, a.sweep_run)
    report(A, D0, Q, B, C)
    res = evaluate(A, D0, Q, B, C)
    print("\nchecks (the operator, then the numerical statements of the supplement, section BdG discretization and convergence)")
    for key, text, ok, detail, literal in res:
        print("  [%s] %s -- %s" % ("FAIL" if not ok else ("PASS" if literal else "PASS, within the rounding of the printed value or the round-off factor"), text, detail))
    nbad = len([1 for r in res if not r[2]])
    loose = [r[0] for r in res if r[2] and not r[4]]
    if loose:
        print("\nNOTE: %d check(s) hold only within the rounding of the printed value or within the round-off factor, not as printed: %s" % (len(loose), ", ".join(loose)))
    print("\nRESULT: %s (%d of %d checks hold; %.0f s)" % ("PASS" if not nbad else "FAIL", len(res) - nbad, len(res), time.time() - t0))
    os.makedirs(a.out, exist_ok=True)
    with open(os.path.join(a.out, "zero_sector_tests.json"), "w", encoding="utf-8", newline="\n") as f:
        json.dump({"part_A": A, "part_A_delta3_0": D0, "part_A_quartet": Q, "part_B": B, "part_C": C,
                   "statements": [{"id": k, "statement": t, "holds": ok, "literal": lit, "detail": d} for k, t, ok, d, lit in res],
                   "environment": {"python": sys.version.split()[0], "numpy": np.__version__}}, f, indent=1, default=float)
        f.write("\n")
    return 0 if not nbad else 1


if __name__ == "__main__":
    sys.exit(main())
