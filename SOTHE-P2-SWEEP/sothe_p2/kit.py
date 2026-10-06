"""Port of the counssug g4a_kit functions (matlab/kit/*.m) and of the suite v0.2.0 Kaup residual.

  tod_soliton        p2_tod_soliton.m 0.1.0       stationary TOD soliton by gauge-fixed Newton
  bdg_spectrum_bg    p2_bdg_spectrum_bg.m 0.1.0   BdG spectrum about an arbitrary background
  gsl_eta15          p2_gsl_eta15.m 0.1.0         Eq. (15) on a trajectory CSV, window values
  splitstep_loss     p2_splitstep_loss.m 0.1.1    direct integration of Eq. (1) with absorbing edges;
                     + sub-grid peak position (SOTHE-P2 default; the MATLAB source records the grid maximum)
  channel_census     p2_channel_census.m 0.1.0    T-12 channel structure in the comoving frame
  gauss_layer        p2_gauss_layer.m 0.1.0       Gaussian layer at a given kappa (T-06)
  fmt_sci            p2_fmt_sci.m 0.1.0           LaTeX number format of the G4C tokens
  kaup_residual      p2_kaup_residual.m 0.2.0     oracle KP-1
RECORDED only (HR-3)."""
import math

import numpy as np

from . import compat as C
from .suite import fd_circulant, _krein_annulus


def tod_soliton(N, d3, n_tau, tau_max, tol=1e-12, maxit=60):
    """p2_tod_soliton.m: solve 1/2 D2 phi + |phi|^2 phi + i d3 D3 phi - mu phi = 0 (periodic, endpoint
    excluded) by Newton on x = [Re phi; Im phi] with the phase and translation gauge rows appended;
    each step is the least-squares solution of the rectangular system (MATLAB backslash)."""
    h = 2.0 * tau_max / n_tau
    tau = -tau_max + np.arange(n_tau) * h
    D2 = fd_circulant(n_tau, h, 4, 2)
    D3 = fd_circulant(n_tau, h, 4, 3)
    mu = 0.5 * N ** 2
    I = np.eye(n_tau)
    psi0 = N * C.sech(N * tau)
    dpsi0 = (-N ** 2) * np.tanh(N * tau) * C.sech(N * tau)
    F0 = 0.5 * (D2 @ psi0) + psi0 ** 3 + (1j * d3) * (D3 @ psi0) - mu * psi0
    res_psi0 = float(np.max(np.abs(F0)))
    p = psi0.copy()
    q = np.zeros(n_tau)
    res = []
    converged = False
    g_phase = np.concatenate([np.zeros(n_tau), psi0])
    g_trans = np.concatenate([dpsi0, np.zeros(n_tau)])
    for _ in range(maxit):
        rho = p ** 2 + q ** 2
        Fr = 0.5 * (D2 @ p) - mu * p - d3 * (D3 @ q) + rho * p
        Fi = 0.5 * (D2 @ q) - mu * q + d3 * (D3 @ p) + rho * q
        F = np.concatenate([Fr, Fi])
        res.append(float(np.max(np.abs(F))))
        if res[-1] < tol:
            converged = True
            break
        J = np.block([[0.5 * D2 - mu * I + np.diag(3 * p ** 2 + q ** 2), -d3 * D3 + np.diag(2 * p * q)],
                      [d3 * D3 + np.diag(2 * p * q), 0.5 * D2 - mu * I + np.diag(p ** 2 + 3 * q ** 2)]])
        A = np.vstack([J, g_phase, g_trans])
        dx = np.linalg.lstsq(A, np.concatenate([-F, [0.0, 0.0]]), rcond=None)[0]
        p = p + dx[:n_tau]
        q = q + dx[n_tau:]
    phi = p + 1j * q
    k = (2 * np.pi / (n_tau * h)) * np.concatenate([np.arange(0, math.ceil(n_tau / 2)),
                                                    np.arange(-math.floor(n_tau / 2), 0)])
    Pw = np.abs(np.fft.fft(phi)) ** 2
    return {"phi": phi, "tau": tau, "mu": mu, "res": res, "res_psi0": res_psi0,
            "kmean": float(np.sum(k * Pw) / np.sum(Pw)), "converged": converged}


def bdg_spectrum_bg(phi, N, d3, n_tau, tau_max):
    """p2_bdg_spectrum_bg.m: BdG spectrum about a complex stationary background phi; w_small = the
    four eigenvalues of smallest modulus (zero-mode sector)."""
    h = 2.0 * tau_max / n_tau
    tau = -tau_max + np.arange(n_tau) * h
    D2 = fd_circulant(n_tau, h, 4, 2)
    D3 = fd_circulant(n_tau, h, 4, 3)
    mu = 0.5 * N ** 2
    phi = np.asarray(phi).ravel()
    H = -0.5 * D2 - 2.0 * np.diag(np.abs(phi) ** 2) + mu * np.eye(n_tau)
    T = (-1j * d3) * D3
    M = np.block([[H + T, -np.diag(phi ** 2)], [np.diag(np.conj(phi) ** 2), -H + T]])
    w, krein, win, ann = _krein_annulus(M, n_tau, mu)
    order = np.argsort(np.abs(w), kind="stable")
    return {"w": w, "krein": krein, "mu": mu, "tau": tau,
            "maxIm": float(np.max(np.abs(w[ann].imag))), "gap": float(np.min(np.abs(w[ann].real))),
            "w_small": w[order[:4]], "ann": ann}


def gsl_eta15(csv_path):
    """p2_gsl_eta15.m: eta15 = sum(max(dS, 0)) / sum|dS| on S_tot, the identity residual, the running
    mean of the eta_GSL column, and the window values for xi >= 1 and xi >= 2."""
    from .suite import read_csv_numeric
    M = read_csv_numeric(csv_path)
    xi, St, eta = M[:, 0], M[:, 3], M[:, 4]
    dS = np.diff(St)
    red = {}
    red["DeltaS"] = float(St[-1] - St[0])
    red["sum_abs"] = float(np.sum(np.abs(dS)))
    red["eta15"] = float(np.sum(np.maximum(dS, 0)) / np.sum(np.abs(dS)))
    red["identity_residual"] = float(abs(red["eta15"] - 0.5 * (1 + red["DeltaS"] / red["sum_abs"])))
    m = ~np.isnan(eta)
    red["eta_mean_running"] = float(np.sum(eta[m]) / np.sum(m))
    red["rise_first4"] = float(St[4] - St[0])
    for x0 in (1, 2):
        i0 = C.first_true(xi >= x0)
        if i0 is None:                                 # MATLAB: empty find -> [] and 0/0 = NaN
            red["DeltaS_xi_ge_%d" % x0] = float("nan")
            red["eta_xi_ge_%d" % x0] = float("nan")
            continue
        d = np.diff(St[i0:])
        red["DeltaS_xi_ge_%d" % x0] = float(St[-1] - St[i0])
        red["eta_xi_ge_%d" % x0] = float(np.sum(np.maximum(d, 0)) / np.sum(np.abs(d)))
    return red


def subgrid_peak(I2, ip, tau, h):
    """Vertex of the parabola through ln I at the grid maximum ip and its two (periodic) neighbours:
    tau_c = tau[ip] + (h/2) (ln I[ip-1] - ln I[ip+1]) / (ln I[ip-1] - 2 ln I[ip] + ln I[ip+1]).
    Exact for a Gaussian peak; |offset| <= h/2 whenever ip is a strict maximum.  Falls back to tau[ip]
    if the three logarithms are not strictly concave (flat or degenerate peak)."""
    n = I2.size
    im, i0, iq = I2[(ip - 1) % n], I2[ip], I2[(ip + 1) % n]
    if im <= 0 or i0 <= 0 or iq <= 0:
        return float(tau[ip])
    lm, l0, lq = math.log(im), math.log(i0), math.log(iq)
    den = lm - 2.0 * l0 + lq
    if not den < 0:
        return float(tau[ip])
    return float(tau[ip] + 0.5 * (lm - lq) / den * h)


def splitstep_loss(N, d3, xi_end=40.0, dxi=0.002, n=8192, L=400.0, noise=0.0, seed=1, keep_field=False, peak="subgrid"):
    """p2_splitstep_loss.m 0.1.1: symmetric split step for d_xi psi = i(1/2 psi_tt + |psi|^2 psi + i d3 psi_ttt),
    absorber exp(-sigma dxi) beyond |tau| > 0.35 L; every 0.5 in xi record the core norm Q, the peak
    position tau_c and the Hann-windowed spectral centroid k_s of the core; fit on xi >= 10.
    noise > 0 is not ported (MATLAB's legacy rand('seed') stream has no NumPy equivalent) and raises.

    Peak position (sub-grid since P2_pypack 2.2.0).  The MATLAB source records tau_c = tau[argmax |psi|^2], which moves in
    steps of h = L/n; window slopes and the recoil check inherit that staircase (Altay run 20260926T221552Z:
    the [30,40] drift at (2, 0.05) jitters by 1.6e-3 between n = 8192, 16384, 32768).  peak="subgrid" (default)
    records tau_c = subgrid_peak(...) instead, which is resolution-independent to ~1e-5; peak="grid" gives the
    0.1.1 behaviour.  Both are always returned: tau_c (the selected one) and tau_c_grid (the grid maximum).
    The core window, Q and k_s stay anchored at the grid maximum, so they are identical in both modes."""
    if peak not in ("subgrid", "grid"):
        raise ValueError("peak must be 'subgrid' or 'grid'")
    if noise > 0:
        raise NotImplementedError("noise > 0: the MATLAB legacy random stream cannot be reproduced in NumPy")
    h = L / n
    tau = -L / 2 + h * np.arange(n)
    k = (2 * np.pi / L) * np.concatenate([np.arange(0, n // 2), np.arange(-n // 2, 0)])
    lin = np.exp(1j * (-0.5 * k ** 2 + d3 * k ** 3) * dxi / 2)
    ta = 0.35 * L
    sig = np.zeros(n)
    m = np.abs(tau) > ta
    sig[m] = 8.0 * ((np.abs(tau[m]) - ta) / (L / 2 - ta)) ** 2
    absorb = np.exp(-sig * dxi)
    psi = (N * C.sech(N * tau)).astype(complex)
    ns = int(C.mround(xi_end / dxi))
    every = int(C.mround(0.5 / dxi))
    rec = []
    fft, ifft = np.fft.fft, np.fft.ifft
    for s in range(ns + 1):
        if s % every == 0:
            I2 = np.abs(psi) ** 2
            ip = int(np.argmax(I2))
            tc = tau[ip]
            core = np.abs(tau - tc) < 6.0 / N
            idx = np.flatnonzero(core)
            nc = idx.size
            w = np.zeros(n)
            w[idx] = 0.5 * (1 - np.cos(2 * np.pi * np.arange(nc) / (nc - 1)))
            Pk = np.abs(fft(psi * w)) ** 2
            ks = np.sum(k * Pk) / np.sum(Pk)
            ts = subgrid_peak(I2, ip, tau, h)
            rec.append((s * dxi, np.sum(I2[core]) * h, ts if peak == "subgrid" else tc, ks, tc))
        if s == ns:
            break
        psi = ifft(lin * fft(psi))
        psi = psi * np.exp(1j * np.abs(psi) ** 2 * dxi)
        psi = ifft(lin * fft(psi))
        psi = psi * absorb
    rec = np.array(rec)
    sel = rec[:, 0] >= 10.0
    pl = np.polyfit(rec[sel, 0], np.log(rec[sel, 1]), 1)
    pd = np.polyfit(rec[sel, 0], rec[sel, 2], 1)
    pg = np.polyfit(rec[sel, 0], rec[sel, 4], 1)
    out = {"xi": rec[:, 0], "Q": rec[:, 1], "tau_c": rec[:, 2], "k_s": rec[:, 3], "tau_c_grid": rec[:, 4],
           "peak": peak, "loss_rate": float(-pl[0]), "drift": float(pd[0]), "drift_grid": float(pg[0]),
           "Q_ratio": float(rec[-1, 1] / rec[0, 1])}
    if keep_field:
        out["tau"] = tau
        out["psi_end"] = psi
    return out


def _realroots(c):
    r = np.roots(c)
    return np.real(r[np.abs(np.imag(r)) < 1e-9])


def channel_census(N, k_op, d3, Vd, kg=0.3, band=(0.05, 1.6), nw=311):
    """p2_channel_census.m 0.1.0: real roots of omega = V k + Omega_pm(k) in the comoving frame, classified
    as low-k (|k| < 1/(3 d3)) or far-detuned; Phi_pair, Phi_up_low, omega_max_pair (closed form)."""
    mu = 0.5 * N ** 2
    f = np.float64                                     # IEEE division as in MATLAB: 1/0 = Inf, 0/0 = NaN
    with np.errstate(divide="ignore", invalid="ignore"):
        b0 = 0.5 * k_op ** 2 - d3 * k_op ** 3
        S = float(f(abs(b0)) / f(math.sqrt(b0 ** 2 + kg ** 2)))
        vg = k_op - 3 * d3 * k_op ** 2
        cm = abs(vg - Vd) * S
        V = 2 * cm
        kinfl = float(f(1) / f(3 * d3))
    w = C.mlinspace(band[0], band[1], nw)
    open_pair = np.zeros(nw, dtype=bool)
    open_up = np.zeros(nw, dtype=bool)
    for i in range(nw):
        rn = _realroots([-d3, -0.5, V, -mu - w[i]])
        ru = _realroots([-d3, 0.5, 0, mu - w[i]])
        open_pair[i] = bool(np.any(np.abs(rn) < kinfl))
        open_up[i] = bool(np.any(np.abs(ru) < kinfl))
    with np.errstate(divide="ignore", invalid="ignore"):
        ks = float(f(math.sqrt(1 + 12 * d3 * V) - 1) / f(6 * d3))
    wmax = V * ks - 0.5 * ks ** 2 - mu - d3 * ks ** 3
    ob = float(np.max(w[open_pair])) if np.any(open_pair) else 0
    return {"N": N, "k_op": k_op, "d3": d3, "Vd": Vd, "kg": kg, "S": S, "c_match": cm, "two_c": V,
            "omega_max_pair": wmax, "omega_max_pair_no_tod": V ** 2 / 2 - mu,
            "Phi_pair": float(np.mean(open_pair)), "Phi_up_low": float(np.mean(open_up)),
            "pair_open_below": ob, "k_star": ks}


def gauss_layer(kappa, wlo=0.05, whi=1.6, n=60):
    """p2_gauss_layer.m 0.1.0: detailed balance with unit flux at kappa; E_N, nu_-, nbar_max; tokens."""
    w = C.mlinspace(wlo, whi, n)
    x = 2 * np.pi * w / kappa
    beta2 = 1 / np.expm1(x)
    alpha2 = 1 + beta2
    EN = 2 * np.arcsinh(np.sqrt(beta2))
    EN_alt = 2 * np.log(np.sqrt(alpha2) + np.sqrt(beta2))
    nu = np.exp(-EN)
    nbar = 0.5 * (np.exp(EN) - 1)
    out = {"kappa": kappa, "w": w, "beta2": beta2, "EN": EN, "nu": nu, "nbar_max_w": nbar,
           "slope": -2 * np.pi / kappa, "EN_max": float(EN[0]), "nu_max": float(nu[-1]),
           "nbar_peak": float(nbar[0]), "nbar_band": float(nbar[-1]),
           "EN_at_3kappa": float(2 * np.arcsinh(np.sqrt(1 / np.expm1(6 * np.pi))))}
    out["IR_asymptote"] = math.log(2 * kappa / (math.pi * wlo))
    out["IR_dev"] = abs(out["EN_max"] - out["IR_asymptote"])
    out["closed_form_agreement"] = float(np.max(np.abs(EN - EN_alt)))
    out["wronskian_residual"] = float(np.max(np.abs(alpha2 - beta2 - 1)))
    out["tokens"] = {"G4A_EN_MAX": C.mfmt("%.4f", out["EN_max"]), "G4A_NU_MAX": C.mfmt("%.4f", out["nu_max"]),
                     "G4A_NBAR_PEAK": C.mfmt("%.2f", out["nbar_peak"]), "G4A_NBAR_BAND": C.mfmt("%.4f", out["nbar_band"])}
    return out


def fmt_sci(x, d):
    """p2_fmt_sci.m: d-decimal mantissa times a power of ten; plain decimal when 0.1 <= |x| < 1000."""
    x = float(x)
    if math.isnan(x):
        return "NaN\\times10^{NaN}"                   # what p2_fmt_sci.m returns for NaN
    if math.isinf(x):
        return "NaN\\times10^{Inf}"                   # ... and for +-Inf (m = Inf/Inf = NaN, e = Inf)
    if x == 0:
        return "0"
    e = math.floor(math.log10(abs(x)))
    m = x / 10 ** e
    if abs(C.mround(m * 10 ** d) / 10 ** d) >= 10:
        e = e + 1
        m = x / 10 ** e
    if -1 <= e <= 2:
        return "%.*f" % (max(d - e, 0), x)
    return "%.*f\\times10^{%d}" % (d, m, e)


def kaup_residual(k=0.7, tau_max=8.0, hs=(0.2, 0.1, 0.05, 0.025)):
    """p2_kaup_residual.m 0.2.0 (oracle KP-1): residual of the exact Kaup continuum eigenfunctions under the
    non-periodic 5-point fourth-order D2 at the interior points; the observed order must approach 4."""
    hs = list(hs)
    res = np.zeros(len(hs))
    Om = 0.5 * (1 + k ** 2)
    for i, h in enumerate(hs):
        n = int(C.mround(2 * tau_max / h)) + 1
        x = C.mlinspace(-tau_max, tau_max, n)
        U = np.exp(1j * k * x) * (k + 1j * np.tanh(x)) ** 2
        V = np.exp(1j * k * x) * C.sech(x) ** 2
        P2 = C.sech(x) ** 2
        j = np.arange(2, n - 2)

        def d2(f):
            return (-f[j + 2] + 16 * f[j + 1] - 30 * f[j] + 16 * f[j - 1] - f[j - 2]) / (12 * h ** 2)

        Hu = -0.5 * d2(U) - 2 * P2[j] * U[j] + 0.5 * U[j]
        Hv = -0.5 * d2(V) - 2 * P2[j] * V[j] + 0.5 * V[j]
        R1 = Hu - P2[j] * V[j] - Om * U[j]
        R2 = P2[j] * U[j] - Hv - Om * V[j]
        res[i] = np.max(np.concatenate([np.abs(R1), np.abs(R2)]))
    return {"k": k, "tau_max": tau_max, "h": hs, "res": res, "order": np.log2(res[:-1] / res[1:]),
            "recorded": [4.121e-3, 2.771e-4, 1.765e-5, 1.108e-6]}
