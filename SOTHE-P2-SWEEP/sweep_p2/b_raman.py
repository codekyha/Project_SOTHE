"""B6: Eq. (1) extended by intrapulse Raman scattering and self-steepening, integrated by RK4IP (Hult, J. Lightwave
Technol. 25, 3770 (2007)), in the convention of the manuscript.

The manuscript's tau is the retarded time reversed (tau = -T/T0 relative to Agrawal's convention; this is what makes
delta3 = beta3 / (6 |beta2| T0) with beta3 > 0 put the resonant radiation at k > 0, the blue side).  In that convention

  d_xi psi = i(1/2 psi_tt + i d3 psi_ttt) + i (1 - i s d_tau) [ psi ((1 - f_R)|psi|^2 + f_R (h * |psi|^2)) ],
  (h * I)(tau) = int_0^inf h(sig) I(tau + sig) dsig,        s = 1 / (omega_0 T0),

h the Blow-Wood response (tau1 = 12.2 fs, tau2 = 32 fs, f_R = 0.18) in units of T0; its transfer function is used in
closed form, G(k) = H_R(k / T0), H_R(w) = [(t1^2 + t2^2)/(t1^2 t2^2)] / [(1/t2 - i w)^2 + 1/t1^2], so the kernel is
exact at every T0.  The linear model is i tau_R psi d_tau|psi|^2 with tau_R = T_R / T0.  Signs are fixed by the tests:
Raman lowers the spectral centroid (red shift, away from the resonant radiation) at Gordon's rate (8/15) tau_R N^4, and
f_R = s = 0 reproduces kit.splitstep_loss.  Absorber, grid, recording and the sub-grid peak are those of kit.splitstep_loss."""
import math

import numpy as np

from sothe_p2 import compat as C
from sothe_p2 import kit

TAU1, TAU2, F_R = 12.2, 32.0, 0.18         # Blow & Wood, IEEE J. Quantum Electron. 25, 2665 (1989)
C_LIGHT = 299792458.0
M1_FS = 2 * TAU1 ** 2 * TAU2 / (TAU1 ** 2 + TAU2 ** 2)   # first moment of h_R: 8.122 fs, T_R = F_R * M1 = 1.462 fs
T_R_LIN_FS = 3.0                                        # the manuscript's estimate T_R / T0 = 0.06 at T0 = 50 fs


def model_params(model, T0_fs, lambda0_nm=1550.0):
    """(f_R, s, tauR_lin) of a named model at T0."""
    s_ss = lambda0_nm * 1e-9 / (2 * math.pi * C_LIGHT * T0_fs * 1e-15)
    if model == "none":
        return 0.0, 0.0, None
    if model == "bw_ss":
        return F_R, s_ss, None
    if model == "bw":
        return F_R, 0.0, None
    if model == "lin3":
        return 0.0, 0.0, T_R_LIN_FS / T0_fs
    if model == "ss":
        return 0.0, s_ss, None
    raise ValueError(model)


def tauR_eff(model, T0_fs):
    if model in ("bw_ss", "bw"):
        return F_R * M1_FS / T0_fs
    if model == "lin3":
        return T_R_LIN_FS / T0_fs
    return 0.0


def raman_G(k, T0_fs):
    w = np.asarray(k, dtype=float) / T0_fs
    return ((TAU1 ** 2 + TAU2 ** 2) / (TAU1 ** 2 * TAU2 ** 2)) / ((1.0 / TAU2 - 1j * w) ** 2 + 1.0 / TAU1 ** 2)


def rk4ip(N, d3, xi_end=40.0, dxi=0.002, n=8192, L=400.0, fR=0.0, s=0.0, T0_fs=None, tauR_lin=None, every=0.5, track=True):
    """Returns the records every 0.5 in xi.  track=True keeps the soliton inside the window: whenever the peak is more
    than 0.1 L from the centre at a record step, the field is rolled by the integer number of grid points that brings the
    peak back (exact on the periodic grid) and tau_c is reported in the fixed frame (offset added).  A Raman-shifted
    soliton at T0 = 50 fs would otherwise reach the absorber (|tau| > 0.35 L) before xi = 40."""
    h = L / n
    tau = -L / 2 + h * np.arange(n)
    k = (2 * np.pi / L) * np.concatenate([np.arange(0, n // 2), np.arange(-n // 2, 0)])
    half = np.exp(1j * (-0.5 * k ** 2 + d3 * k ** 3) * dxi / 2)
    ta = 0.35 * L
    sig = np.zeros(n)
    m = np.abs(tau) > ta
    sig[m] = 8.0 * ((np.abs(tau[m]) - ta) / (L / 2 - ta)) ** 2
    absf = np.exp(-sig * dxi)
    fft, ifft = np.fft.fft, np.fft.ifft
    G = raman_G(k, T0_fs) if fR > 0 else None
    ik = 1j * k

    def NL(p):
        I = np.abs(p) ** 2
        if tauR_lin is not None:
            Ie = I + tauR_lin * np.real(ifft(ik * fft(I)))
        elif fR > 0:
            Ie = (1 - fR) * I + fR * np.real(ifft(fft(I) * G))
        else:
            Ie = I
        P = p * Ie
        if s != 0:
            P = P - 1j * s * ifft(ik * fft(P))
        return 1j * P

    psi = (N * C.sech(N * tau)).astype(complex)
    ns = int(C.mround(xi_end / dxi))
    ev = int(C.mround(every / dxi))
    rec = []
    offset, n_roll = 0.0, 0
    for st in range(ns + 1):
        if st % ev == 0:
            I2 = np.abs(psi) ** 2
            ip = int(np.argmax(I2))
            tc = tau[ip]
            core = np.abs(tau - tc) < 6.0 / N
            idx = np.flatnonzero(core)
            nc = idx.size
            w = np.zeros(n)
            w[idx] = 0.5 * (1 - np.cos(2 * np.pi * np.arange(nc) / (nc - 1)))
            Pk = np.abs(fft(psi * w)) ** 2
            Pf = np.abs(fft(psi)) ** 2
            rec.append((st * dxi, float(np.sum(I2[core]) * h), kit.subgrid_peak(I2, ip, tau, h) + offset, float(np.sum(k * Pk) / np.sum(Pk)),
                        float(np.sum(k * Pf) / np.sum(Pf)), float(np.sum(I2) * h)))
            if track and abs(tc) > 0.1 * L:
                sh = int(C.mround(tc / h))
                psi = np.roll(psi, -sh)
                offset += sh * h
                n_roll += 1
        if st == ns:
            break
        Af = fft(psi)
        AI = ifft(half * Af)
        k1 = ifft(half * fft(dxi * NL(psi)))
        k2 = dxi * NL(AI + 0.5 * k1)
        k3 = dxi * NL(AI + 0.5 * k2)
        k4 = dxi * NL(ifft(half * fft(AI + k3)))
        psi = ifft(half * fft(AI + k1 / 6 + k2 / 3 + k3 / 3)) + k4 / 6
        psi = psi * absf
    rec = np.array(rec)
    return {"xi": rec[:, 0], "Q": rec[:, 1], "tau_c": rec[:, 2], "k_s": rec[:, 3], "k_full": rec[:, 4], "P": rec[:, 5],
            "n_roll": n_roll, "offset": offset}


def _rate(xi, y, lo, hi):
    sw = (xi >= lo) & (xi <= hi)
    if np.sum(sw) < 2:
        return float("nan")
    return float(np.polyfit(xi[sw], y[sw], 1)[0])


def windows_for(xi_end):
    base = [(0, 5), (5, 10), (10, 20), (20, 30), (30, 40)]
    w = [(lo, min(hi, xi_end)) for lo, hi in base if lo < xi_end]
    return [x for x in w if x[1] - x[0] >= 1.0]


def run_model(N, d3, T0_fs, model, P):
    """One propagation; records every P['every'] (default 0.5) in xi.  window_drift_05 uses the records at multiples
    of 0.5 only (the sampling of B3), for the comparison with the split step."""
    fR, s, tl = model_params(model, T0_fs, P.get("lambda0_nm", 1550.0))
    every = P.get("every", 0.5)
    X = rk4ip(N, d3, P["xi_end"], P["dxi"], P["n"], P["L"], fR=fR, s=s, T0_fs=T0_fs, tauR_lin=tl, every=every)
    W = windows_for(P["xi_end"])
    xi = X["xi"]
    sub = slice(None, None, max(1, int(C.mround(0.5 / every))))
    tre = tauR_eff(model, T0_fs)
    out = {"N": N, "d3": d3, "T0_fs": T0_fs, "model": model, "f_R": fR, "s": s, "tauR_lin": tl, "tauR_eff": tre,
           "gordon_rate": -(8.0 / 15.0) * tre * N ** 4, "windows": W,
           "window_drift": [_rate(xi, X["tau_c"], lo, hi) for lo, hi in W],
           "window_drift_05": [_rate(xi[sub], X["tau_c"][sub], lo, hi) for lo, hi in W], "every": every,
           "window_ks_rate": [_rate(xi, X["k_s"], lo, hi) for lo, hi in W],
           "window_kfull_rate": [_rate(xi, X["k_full"], lo, hi) for lo, hi in W],
           "loss_end": float(1 - X["Q"][-1] / X["Q"][0]), "k_s_end": float(X["k_s"][-1]), "k_full_end": float(X["k_full"][-1]),
           "tau_c_end": float(X["tau_c"][-1]), "xi_end": P["xi_end"], "dxi": P["dxi"], "n": P["n"], "L": P["L"],
           "n_roll": X["n_roll"], "offset": X["offset"],
           "series": {k: X[k].tolist() for k in ("xi", "Q", "tau_c", "k_s", "k_full", "P")}}
    fitw = (0.0, min(5.0, P["xi_end"]))
    out["early_fit_window"] = list(fitw)
    out["early_kfull_rate"] = _rate(xi, X["k_full"], *fitw)
    out["early_ks_rate"] = _rate(xi, X["k_s"], *fitw)
    out["initial_fit_window"] = [0.0, 0.5]
    out["initial_ks_rate"] = _rate(xi, X["k_s"], 0.0, 0.5 + 1e-9)
    out["initial_kfull_rate"] = _rate(xi, X["k_full"], 0.0, 0.5 + 1e-9)
    out["shift_total_ks"] = float(X["k_s"][-1] - X["k_s"][0])
    return out


def gordon_check(tauR_lin=None, T0_bw=None, xi_end=10.0, dxi=0.002, n=8192, L=400.0):
    """N = 1, delta3 = 0: full-field centroid rate over xi in [1, 5] against -(8/15) tau_R (Gordon 1986)."""
    if tauR_lin is not None:
        X = rk4ip(1.0, 0.0, xi_end, dxi, n, L, tauR_lin=tauR_lin)
        tre, label = tauR_lin, "linear tau_R = %g" % tauR_lin
    else:
        X = rk4ip(1.0, 0.0, xi_end, dxi, n, L, fR=F_R, T0_fs=T0_bw)
        tre, label = F_R * M1_FS / T0_bw, "Blow-Wood at T0 = %g fs" % T0_bw
    lo, hi = 1.0, min(5.0, xi_end)
    r = _rate(X["xi"], X["k_full"], lo, hi)
    g = -(8.0 / 15.0) * tre
    return {"label": label, "tauR_eff": tre, "rate_full": r, "rate_windowed": _rate(X["xi"], X["k_s"], lo, hi), "gordon": g,
            "ratio": r / g, "fit_window": [lo, hi], "photon_drift": float(abs(X["P"][min(len(X["P"]) - 1, 10)] / X["P"][0] - 1))}
