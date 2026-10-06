"""Port of sothe_phase2_matlab v0.1.0 (the Phase-2 physics suite, matlab/suite/*.m) and of the v0.2.0
parameter overlay (matlab/suite_v020/p2_phase2_params_v020.m).

One function per MATLAB file, same grids, stencils, formulas and reductions; the MATLAB file name is given
in each docstring.  Grids use MATLAB's linspace formula (compat.mlinspace), so every grid column of the CSV
products is identical to the MATLAB R2025b one.  RECORDED only (HR-3)."""
import math
import os

import numpy as np

from . import compat as C
from .par import pmap

STAB_ENV = ("PHASE2_STAB_NTAU", "PHASE2_STAB_GRIDN", "PHASE2_STAB_GRIDD")
# the stability map feeds no oracle; stages S4 and S5 run it on a 2 x 2 grid at n_tau = 48 (as the MATLAB stages do)
SMALL_STAB = {"PHASE2_STAB_NTAU": "48", "PHASE2_STAB_GRIDN": "2", "PHASE2_STAB_GRIDD": "2"}


# ============================================================ parameters
def phase2_params(env=None):
    """p2_phase2_params.m v0.1.0: operating points and FROZEN targets (29/29 oracle set).

    ENV (default os.environ) supplies the optional stability-map overrides PHASE2_STAB_NTAU, PHASE2_STAB_GRIDN
    and PHASE2_STAB_GRIDD, read as the MATLAB file reads them with getenv.  Stages pass an explicit mapping, so
    nothing is written into the process environment."""
    env = os.environ if env is None else env
    P = {}
    P["N"] = 2.0
    P["d3"] = 0.05
    P["kg"] = 0.3
    P["k_op"] = 1.0
    P["drive0"] = 1.0
    P["omega"] = C.mlinspace(0.05, 1.60, 60)
    P["drive"] = C.mlinspace(0.70, 1.30, 9)
    P["kg_flow"] = C.mlinspace(0.0, 1.0, 21)
    P["tau_kin"] = {"lo": -20, "hi": 20, "n": 4096}
    P["bdg_N"] = 1.0
    P["bdg_d3"] = [0.0, 0.02]
    P["bdg_n_tau"] = 500
    P["bdg_tau_max"] = 16.0
    P["rate_d3"] = [0.0, 0.01, 0.02, 0.03, 0.05, 0.08]
    P["stab_N"] = C.mlinspace(0.6, 2.4, 13)
    P["stab_d3"] = C.mlinspace(0.0, 0.08, 13)
    P["stab_n_tau"] = 360
    P["stab_tau_max"] = 16.0
    ov = env.get("PHASE2_STAB_NTAU", "")
    if ov:
        P["stab_n_tau"] = int(C.mround(float(ov)))
    gN = env.get("PHASE2_STAB_GRIDN", "")
    if gN:
        P["stab_N"] = C.mlinspace(0.6, 2.4, int(C.mround(float(gN))))
    gD = env.get("PHASE2_STAB_GRIDD", "")
    if gD:
        P["stab_d3"] = C.mlinspace(0.0, 0.08, int(C.mround(float(gD))))
    P["fp_kg"] = C.mlinspace(0.0, 1.0, 61)
    P["fp_drive"] = C.mlinspace(0.0, 1.5, 61)
    P["cont_N"] = [1, 2, 3]
    T = {}
    T["kappa_kinematic"] = 1.6635880450803207
    T["T_H"] = 0.26476826064311587
    T["thermality_R2"] = 1.0
    T["thermality_slope"] = -3.7768877492000867
    T["flux_residual"] = 2.220446049250313e-16
    T["dS_tot"] = 1.814284948647733
    T["eta_GSL"] = 0.6093410449799196
    T["E_N_max_nats"] = 3.053868887353342
    T["E_N_at_3kappa"] = 1.6139903549101825e-04
    T["bdg_gap_over_mu_d3_0"] = 0.9999997734961216
    T["bdg_maxIm_d3_0"] = 1.1575664769089186e-16
    T["bdg_maxIm_d3_002"] = 0.09572703807379634
    T["kappa_drive_spread"] = 9.10e-05
    T["PT_nu_minus_max"] = 0.9071
    T["greybody_lo"] = 0.0071
    T["greybody_hi"] = 0.9774
    T["kappa_cont"] = [0.9999, 1.9994, 2.9979]
    T["rate_table"] = [1.16e-16, 0.06798, 0.09573, 0.11681, 0.14982, 0.18450]
    P["T"] = T
    P["tol_eig"] = 1e-6
    P["tol_det"] = 1e-9
    P["tol_exact"] = 1e-12
    P["seed"] = 20260628
    return P


def phase2_params_v020(env=None):
    """p2_phase2_params_v020.m 0.2.0: v0.1.0 parameters with the T-14 oracle changes (retired pins moved)."""
    P = phase2_params(env)
    P["suite_version"] = "0.2.0"
    T = dict(P["T"])
    X = {}
    X["bdg_maxIm_d3_002"] = {"value": T["bdg_maxIm_d3_002"],
                             "reason": "pins the psi0-linearization artifact (S1-1); retired with P1.4"}
    X["rate_table_d3_gt_0"] = {"value": T["rate_table"][1:], "delta3": [0.01, 0.02, 0.03, 0.05, 0.08],
                               "reason": "same artifact (S1-1); retired with P1.5"}
    X["greybody_lo"] = {"value": T["greybody_lo"], "reason": "Poeschl-Teller grey-body pin (S1-2); retired with P3.3"}
    X["greybody_hi"] = {"value": T["greybody_hi"], "reason": "Poeschl-Teller grey-body pin (S1-2)"}
    for k in ("bdg_maxIm_d3_002", "greybody_lo", "greybody_hi"):
        del T[k]
    T["rate_table"] = T["rate_table"][0]
    T["eta15"] = 0.5452
    T["eta_xi_ge_2"] = 0.500073523187179
    T["DeltaS_xi_ge_2"] = 0.002189246165536396
    T["nsafe_factor"] = 7 / 12
    T["loss_band"] = [2e-4, 6e-4]
    T["kaup_min_order"] = 3.7
    T["kaup_recorded"] = [4.121e-3, 2.771e-4, 1.765e-5, 1.108e-6]
    T["tod_newton_tol"] = 1e-12
    T["zero_sector_tol"] = 1e-6
    T["annulus_tol"] = 1e-10
    T["bg_regression_tol"] = 1e-10
    P["T"] = T
    P["T_retired"] = X
    P["tau_kin_dtau"] = (P["tau_kin"]["hi"] - P["tau_kin"]["lo"]) / (P["tau_kin"]["n"] - 1)
    return P


# ============================================================ BdG
def fd_circulant(n, h, order, m):
    """p2_fd_circulant.m: periodic central FD matrix of the m-th derivative, formal accuracy `order`."""
    half = (order + m - 1) // 2
    offs = np.arange(-half, half + 1)
    npt = offs.size
    pw = np.arange(npt).reshape(-1, 1)
    A = offs.reshape(1, -1).astype(float) ** pw          # A[p, j] = offs[j]^p
    b = np.zeros(npt)
    b[m] = math.factorial(m)
    w = np.linalg.solve(A, b)
    D = np.zeros((n, n))
    idx = np.arange(n)
    for t in range(npt):
        cols = np.mod(idx + offs[t], n)
        np.add.at(D, (idx, cols), w[t])
    return D / h ** m


def soliton_mu(N, tau=None):
    """p2_soliton_mu.m: mu = N^2/2 and psi0 = N sech(N tau)."""
    mu = 0.5 * N ** 2
    return mu, (None if tau is None else N * C.sech(N * np.asarray(tau)))


def _krein_annulus(M, n_tau, mu):
    w, V = np.linalg.eig(M)
    u = V[:n_tau, :]
    v = V[n_tau:, :]
    krein = np.sign(np.sum(np.abs(u) ** 2, axis=0) - np.sum(np.abs(v) ** 2, axis=0))
    win = np.abs(w.real) <= 3.0 * mu
    ann = win & (np.abs(w.real) > 0.1 * mu)
    return w, krein, win, ann


def bdg_spectrum(N, d3, n_tau, tau_max):
    """p2_bdg_spectrum.m: BdG spectrum about psi0 = N sech(N tau); periodic grid, endpoint excluded."""
    tau = -tau_max + np.arange(n_tau) * (2.0 * tau_max / n_tau)
    h = 2.0 * tau_max / n_tau
    D2 = fd_circulant(n_tau, h, 4, 2)
    D3 = fd_circulant(n_tau, h, 4, 3)
    psi2 = (N * C.sech(N * tau)) ** 2
    mu = 0.5 * N ** 2
    P2 = np.diag(psi2)
    H = -0.5 * D2 - 2.0 * P2 + mu * np.eye(n_tau)
    T = (-1j * d3) * D3
    M = np.block([[H + T, -P2], [P2, -H + T]])
    w, krein, win, ann = _krein_annulus(M, n_tau, mu)
    return {"w": w, "krein": krein, "mu": mu, "tau": tau,
            "maxIm": float(np.max(np.abs(w[ann].imag))), "gap": float(np.min(np.abs(w[ann].real))),
            "win": win, "ann": ann}


# ============================================================ kinematics
def beta0(k, d3):
    """p2_beta0.m: beta0(k) = k^2/2 - delta3 k^3."""
    return 0.5 * k ** 2 - d3 * k ** 3


def vg0(k, d3):
    """p2_vg0.m: vg0(k) = k - 3 delta3 k^2."""
    return k - 3.0 * d3 * k ** 2


def slow_light(k, d3, kg):
    """p2_slow_light.m: S = |beta0| / sqrt(beta0^2 + kappa_g^2)."""
    b = beta0(k, d3)
    return np.abs(b) / np.sqrt(b ** 2 + kg ** 2)


def vg_dressed(k, d3, kg):
    """p2_vg_dressed.m."""
    return vg0(k, d3) * slow_light(k, d3, kg)


def flow_profile(tau, c_match, kappa0, tau_h):
    """p2_flow_profile.m: v = c_match (1 + tanh(kappa0 (tau - tau_h)/c_match))."""
    return c_match * (1.0 + np.tanh(kappa0 * (tau - tau_h) / c_match))


def grad_uniform(v, h):
    """p2_grad_uniform.m: numpy.gradient on a uniform grid (2nd-order interior and ends)."""
    v = np.asarray(v, dtype=float).ravel()
    n = v.size
    g = np.zeros(n)
    g[1:n - 1] = (v[2:n] - v[0:n - 2]) / (2 * h)
    g[0] = (-3 * v[0] + 4 * v[1] - v[2]) / (2 * h)
    g[n - 1] = (3 * v[n - 1] - 4 * v[n - 2] + v[n - 3]) / (2 * h)
    return g


def surface_gravity(x, v, c):
    """p2_surface_gravity.m: |dv/dx| at the first linear-interpolated root of v - c."""
    x = np.asarray(x, dtype=float).ravel()
    v = np.asarray(v, dtype=float).ravel()
    g = v - c
    s = np.sign(g)
    ic = C.first_true(np.diff(s) != 0)
    if ic is None:
        return 0.0, float("nan")
    x0, x1, g0, g1 = x[ic], x[ic + 1], g[ic], g[ic + 1]
    x_h = x0 - g0 * (x1 - x0) / (g1 - g0)
    h = x[1] - x[0]
    dv = grad_uniform(v, h)
    return abs(C.interp1(x, dv, x_h)), float(x_h)


def horizon_strength(kg, drive):
    """p2_horizon_strength.m (diagnostic only; INV1 forbids it for kappa)."""
    return kg * drive * 1.0


def kinematic_kappa(N, d3, kg, drive, k_op=1.0):
    """p2_kinematic_kappa.m: drive-independent kinematic surface gravity (INV1: kinematic only)."""
    Sv = float(slow_light(k_op, d3, kg))
    c_match = max(abs(float(vg0(k_op, d3))) * Sv, 1e-3)
    kappa0 = N * Sv
    tau_h = -(drive - 1.0) / (2.0 * N)
    tau = C.mlinspace(-20, 20, 4096)
    v = flow_profile(tau, c_match, kappa0, tau_h)
    kappa, x_h = surface_gravity(tau, v, c_match)
    return kappa, {"Sv": Sv, "c_match": c_match, "kappa0": kappa0, "tau_h": tau_h, "x_h": x_h}


# ============================================================ thermal layer
def bogoliubov_thermal(omega, kappa):
    """p2_bogoliubov_thermal.m: ratio = exp(-2 pi omega/kappa), |alpha|^2 = 1/(1 - ratio), |beta|^2 = ratio |alpha|^2."""
    ratio = np.exp(-2.0 * np.pi * np.asarray(omega, dtype=float) / kappa)
    alpha2 = 1.0 / (1.0 - ratio)
    beta2 = alpha2 * ratio
    return alpha2, beta2, ratio


def flux_residual(alpha2, beta2):
    """p2_flux_residual.m."""
    return float(np.max(np.abs(alpha2 - beta2 - 1.0)))


def thermality_fit(omega, ratio):
    """p2_thermality_fit.m: unweighted LSQ of ln(ratio) on omega by the normal equations."""
    omega = np.asarray(omega, dtype=float).ravel()
    y = np.log(np.asarray(ratio, dtype=float).ravel())
    A = np.column_stack([omega, np.ones(omega.size)])
    coef = np.linalg.solve(A.T @ A, A.T @ y)
    slope, intercept = float(coef[0]), float(coef[1])
    yhat = A @ coef
    ss_res = np.sum((y - yhat) ** 2)
    ss_tot = np.sum((y - np.mean(y)) ** 2)
    R2 = float(1.0 - ss_res / ss_tot)
    return slope, intercept, R2, -2.0 * np.pi / slope


def greybody(omega, N):
    """p2_greybody.m: Poeschl-Teller transmission, clipped to [0, 1]."""
    V0 = 0.5 * N ** 2
    s = math.sqrt(1.0 + 8.0 * V0 / N ** 2)
    k = np.maximum(np.asarray(omega, dtype=float), 1e-9)
    sh = np.sinh(np.pi * k / N) ** 2
    G = sh / (sh + math.cos(math.pi * s / 2.0) ** 2)
    return np.minimum(np.maximum(G, 0.0), 1.0)


def partner_logneg(omega, kappa):
    """p2_partner_logneg.m: E_N = 2 asinh|beta|, E_N_cross = 2 ln(|alpha| + |beta|), nu_- = exp(-2 r)."""
    alpha2, beta2, _ = bogoliubov_thermal(omega, kappa)
    b = np.sqrt(beta2)
    a = np.sqrt(alpha2)
    r = np.arcsinh(b)
    return 2.0 * r, 2.0 * np.log(a + b), np.exp(-2.0 * r)


def no_horizon_null(kg, drive):
    """p2_no_horizon_null.m: |beta| = |tanh(kappa_g drive)|."""
    return np.abs(np.tanh(np.asarray(kg) * np.asarray(drive)))


# ============================================================ trajectory
def read_csv_numeric(path):
    """p2_read_csv_numeric.m: numeric matrix of a comma-separated file with one text header row."""
    rows = []
    first = True
    with open(path, encoding="utf-8") as f:
        for ln in f:
            ln = ln.rstrip("\r\n")
            if not ln.strip():
                continue
            parts = ln.split(",")
            if first:
                first = False
                if math.isnan(_num(parts[0])):
                    continue
            rows.append([_num(p) for p in parts])
    return np.array(rows, dtype=float)


def _num(s):
    s = s.strip()
    try:
        return float(s)
    except ValueError:
        return float("nan")


def entropy_reductions(csv_path):
    """p2_entropy_reductions.m: frozen GSL-trajectory reductions (P5)."""
    M = read_csv_numeric(csv_path)
    S_hor, S_rad, S_tot, eta, P = M[:, 1], M[:, 2], M[:, 3], M[:, 4], M[:, 5]
    m = ~np.isnan(eta)
    return {"dS_tot": float(S_tot[-1] - S_tot[0]),
            "mean_eta_GSL": float(np.sum(eta[m]) / np.sum(m)),
            "partition_res": float(np.max(np.abs(S_tot - (S_hor + S_rad)))),
            "norm_drift": float(np.max(np.abs(P / P[0] - 1.0))),
            "n": int(S_tot.size)}


# ============================================================ one physics pass
def _bdg_task(N, d3, n_tau, tau_max):
    return bdg_spectrum(N, d3, n_tau, tau_max)


_COMPUTE_CACHE = {}


def compute_all(P, csv_path):
    """p2_compute_all.m v0.1.0: one deterministic physics pass reused by the oracles and the figures.

    The BdG spectra (eigenvalue problems of size 2 n_tau) run in parallel; identical (N, delta3, n_tau)
    requests (the delta3 = 0 and 0.02 spectra of R.bdg and of the rate table) are computed once.
    The result is cached per process for the same parameters and trajectory file."""
    key = (repr({k: (v.tolist() if isinstance(v, np.ndarray) else v) for k, v in P.items() if k not in ("T", "T_retired")}),
           os.path.abspath(csv_path))
    if key in _COMPUTE_CACHE:
        return _COMPUTE_CACHE[key]
    R = {}
    kappa_kin, kin_aux = kinematic_kappa(P["N"], P["d3"], P["kg"], P["drive0"], P["k_op"])
    R["kappa_kin"] = kappa_kin
    R["kin_aux"] = kin_aux
    R["T_H"] = kappa_kin / (2 * np.pi)
    om = np.asarray(P["omega"], dtype=float).ravel()
    a2, b2, ratio = bogoliubov_thermal(om, kappa_kin)
    R.update(omega=om, alpha2=a2, beta2=b2, ratio=ratio)
    R["flux_residual"] = flux_residual(a2, b2)
    sl, ic, R2, kfit = thermality_fit(om, ratio)
    R.update(slope=sl, intercept=ic, R2=R2, kappa_fit=kfit)
    m0 = -2 * np.pi / 1.3
    sl_s, _, R2_s, _ = thermality_fit(om, np.exp(m0 * om))
    R["synfit_slope_err"] = abs(sl_s - m0)
    R["synfit_R2"] = R2_s
    R["Gamma"] = greybody(om, P["N"])
    EN, ENx, num = partner_logneg(om, kappa_kin)
    R.update(E_N=EN, E_N_cross=ENx, nu_minus=num)
    R["E_N_max"] = float(np.max(EN))
    EN3, _, _ = partner_logneg(3 * kappa_kin, kappa_kin)
    R["E_N_at_3kappa"] = float(EN3)
    R["EN_asymptote_lo"] = np.log(2 * kappa_kin / (np.pi * om))
    kd = np.array([kinematic_kappa(P["N"], P["d3"], P["kg"], d, P["k_op"])[0] for d in P["drive"]])
    R["drive"] = np.asarray(P["drive"], dtype=float).ravel()
    R["kappa_drive"] = kd
    R["kappa_drive_spread"] = float(np.max(kd) - np.min(kd))
    R["kg_flow"] = np.asarray(P["kg_flow"], dtype=float).ravel()
    R["kappa_flow"] = np.array([kinematic_kappa(P["N"], P["d3"], g, P["drive0"], P["k_op"])[0] for g in P["kg_flow"]])
    R["kappa_cont"] = np.array([kinematic_kappa(n, P["d3"], 0.0, P["drive0"], P["k_op"])[0] for n in P["cont_N"]])
    KG, DR = np.meshgrid(P["fp_kg"], P["fp_drive"])
    R["fp_kg"] = np.asarray(P["fp_kg"], dtype=float)
    R["fp_drive"] = np.asarray(P["fp_drive"], dtype=float)
    R["fp_beta"] = no_horizon_null(KG, DR)
    R["null_a"] = float(no_horizon_null(0.0, 1.0))
    R["null_b"] = float(no_horizon_null(0.5, 0.0))
    R["offnull"] = float(no_horizon_null(0.3, 0.4))
    # BdG: the canonical-soliton spectra, the radiation-rate table and the stability map, in parallel
    jobs = []
    for d3 in list(P["bdg_d3"]) + list(P["rate_d3"]):
        j = (P["bdg_N"], d3, P["bdg_n_tau"], P["bdg_tau_max"])
        if j not in jobs:
            jobs.append(j)
    for N in P["stab_N"]:
        for d3 in P["stab_d3"]:
            j = (float(N), float(d3), P["stab_n_tau"], P["stab_tau_max"])
            if j not in jobs:
                jobs.append(j)
    out = dict(zip(jobs, pmap(_bdg_task, jobs, progress="BdG spectra" if len(jobs) > 20 else None)))
    R["bdg"] = [{"d3": d3, "out": out[(P["bdg_N"], d3, P["bdg_n_tau"], P["bdg_tau_max"])]} for d3 in P["bdg_d3"]]
    R["rate_d3"] = np.asarray(P["rate_d3"], dtype=float)
    R["rate_table"] = np.array([out[(P["bdg_N"], d3, P["bdg_n_tau"], P["bdg_tau_max"])]["maxIm"] for d3 in P["rate_d3"]])
    Smap = np.zeros((len(P["stab_d3"]), len(P["stab_N"])))
    for iN, N in enumerate(P["stab_N"]):
        for iD, d3 in enumerate(P["stab_d3"]):
            Smap[iD, iN] = out[(float(N), float(d3), P["stab_n_tau"], P["stab_tau_max"])]["maxIm"]
    R["stab_N"] = np.asarray(P["stab_N"], dtype=float)
    R["stab_d3"] = np.asarray(P["stab_d3"], dtype=float)
    R["stab_map"] = Smap
    R["red"] = entropy_reductions(csv_path)
    _COMPUTE_CACHE[key] = R
    return R
