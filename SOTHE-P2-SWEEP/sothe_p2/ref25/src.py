"""Python ports of the MATLAB functions of Ref. [25] (SOTHE_pkg v1.0.0, src/*.m; matlab/ref25/ of P2_altay_pack).

  soliton_mask          soliton_mask.m          centroid-tracking super-Gaussian system/bath mask
  spectral_entropy      spectral_entropy.m      differential Shannon entropy of a masked, normalized spectrum
  eta_gsl               eta_gsl.m               running entropy-production efficiency, Eq. (8) of Ref. [25]
  find_rr_lobe          find_rr_lobe.m          resonant-radiation lobe centre on the smoothed log spectrum
  rr_phase_match_fit    rr_phase_match_fit.m    through-origin fit k_RR = c_RR / delta3
  gnlse_dimensionless   gnlse_dimensionless.m   symmetric SSFM, dimensionless GNLSE, Eq. (9) of Ref. [25]
  gnlse_dimensional     gnlse_dimensional.m     symmetric SSFM in physical units (corrected beta3 sign)
  run_robustness_sweep  run_robustness_sweep.m  the 5 x 5 x 4 x 3 = 300-configuration ensemble (stage S7)

These follow the MATLAB files statement by statement: MATLAB's linspace, round, trapz and cumtrapz
(compat.mlinspace, mround_array, mtrapz, mcumtrapz), find_rr_lobe's own reflect padding and 3-sigma Gaussian
(the validation port uses scipy.ndimage.gaussian_filter1d, which pads and truncates differently), the input
checks and error messages.  The package's unchanged Python validation port (ref25/gnlse_dimensionless.py) stays
the solver of stage S3, which must reproduce the archived trajectory with the arithmetic that wrote it.
Author of the MATLAB originals: H. Oguz, 2026 (MIT, ref25/LICENSE).  RECORDED only (HR-3)."""
import math
import time

import numpy as np

from .. import compat as C


class MatlabError(ValueError):
    """An error the MATLAB function raises (identifier and message kept)."""

    def __init__(self, ident, msg):
        super().__init__("%s: %s" % (ident, msg))
        self.identifier = ident


# ------------------------------------------------------------------ masks, entropy, efficiency
def soliton_mask(k, P, width, order=10):
    """[mask_S, mask_R, k_c] = soliton_mask(k, P, width, order): M_S = exp(-((k - k_c)/width)^order),
    k_c = trapz(k, k P) / trapz(k, P) (0 if the spectrum carries no weight), M_R = 1 - M_S."""
    if order is None:
        order = 10
    k = np.asarray(k, dtype=float).ravel()
    P = np.asarray(P, dtype=float).ravel()
    Ptot = C.mtrapz(k, P)
    if Ptot < 1e-20:
        k_c = 0.0
    else:
        k_c = float(C.mtrapz(k, k * P) / Ptot)
    mask_S = np.exp(-((k - k_c) / width) ** order)
    return mask_S, 1 - mask_S, k_c


def spectral_entropy(k, P, mask):
    """S = spectral_entropy(k, P, mask) in nats: -trapz over p > 1e-20 of p ln p, p = max(P mask, 0) / Z."""
    k = np.asarray(k, dtype=float).ravel()
    P = np.asarray(P, dtype=float).ravel()
    mask = np.asarray(mask, dtype=float).ravel()
    if k.size != P.size or k.size != mask.size:
        raise MatlabError("spectral_entropy:dim", "k, P, mask must all be vectors of the same length.")
    if np.any(P < -np.finfo(float).eps * np.max(P)):
        raise MatlabError("spectral_entropy:neg", "Power spectrum P must be non-negative.")
    P_sub = np.maximum(P * mask, 0)
    Z = C.mtrapz(k, P_sub)
    if Z < 1e-20 * np.max(P):
        return 0.0
    p = P_sub / Z
    valid = p > 1e-20
    return float(-C.mtrapz(k[valid], p[valid] * np.log(p[valid])))


def eta_gsl(xi, S_tot):
    """eta = eta_gsl(xi, S_tot): cumulative trapezoidal ratio of max(0, dS/dxi) to |dS/dxi| (centred differences,
    one-sided at the ends); eta(1) = 0."""
    xi = np.asarray(xi, dtype=float).ravel()
    S_tot = np.asarray(S_tot, dtype=float).ravel()
    N = xi.size
    if N != S_tot.size:
        raise MatlabError("eta_gsl:dim", "xi and S_tot must have the same length.")
    if N < 3:
        raise MatlabError("eta_gsl:short", "Need at least 3 samples to compute eta_GSL.")
    if np.any(np.diff(xi) <= 0):
        raise MatlabError("eta_gsl:order", "xi must be strictly increasing.")
    dS = np.zeros(N)
    dS[1:N - 1] = (S_tot[2:N] - S_tot[0:N - 2]) / (xi[2:N] - xi[0:N - 2])
    dS[0] = (S_tot[1] - S_tot[0]) / (xi[1] - xi[0])
    dS[N - 1] = (S_tot[N - 1] - S_tot[N - 2]) / (xi[N - 1] - xi[N - 2])
    num = C.mcumtrapz(xi, np.maximum(dS, 0))
    denom = C.mcumtrapz(xi, np.abs(dS))
    eta = np.zeros(N)
    nz = denom > 0
    eta[nz] = num[nz] / denom[nz]
    return eta


def find_rr_lobe(k, P, k_c, mask_width, side="positive", smoothing_sigma=3.0, search_radius_frac=0.15):
    """k_RR = find_rr_lobe(k, P, k_c, mask_width, 'side', ..., 'smoothing_sigma', ..., 'search_radius_frac', ...):
    argmax of the Gaussian-smoothed log spectrum outside k_c +- 3 mask_width, then the centroid of the lobe."""
    k = np.asarray(k, dtype=float).ravel()
    P = np.asarray(P, dtype=float).ravel()
    s = str(side).lower()
    if s == "positive":
        sel = k > k_c + 3.0 * mask_width
    elif s == "negative":
        sel = k < k_c - 3.0 * mask_width
    else:
        raise MatlabError("find_rr_lobe:side", "side must be 'positive' or 'negative'.")
    sel = sel & (P > np.max(P) * 1e-6)
    if not np.any(sel):
        return float("nan")
    logP = np.log(P + 1e-30)
    win = max(1, int(math.ceil(3 * smoothing_sigma)))
    xx = np.arange(-win, win + 1, dtype=float)
    g = np.exp(-xx ** 2 / (2 * smoothing_sigma ** 2))
    g = g / np.sum(g)
    pad = np.concatenate([logP[1:win + 1][::-1], logP, logP[-win - 1:-1][::-1]])   # MATLAB: flipud(logP(2:win+1)) ...
    logP_smooth = np.convolve(pad, g, mode="valid")
    P_search = np.full(P.shape, -np.inf)
    P_search[sel] = logP_smooth[sel]
    j0 = int(np.argmax(P_search))
    k0 = k[j0]
    r = search_radius_frac * abs(k0)
    lobe = (np.abs(k - k0) < r) & sel
    if not np.any(lobe):
        return float(k0)
    return float(C.mtrapz(k[lobe], k[lobe] * P[lobe]) / C.mtrapz(k[lobe], P[lobe]))


def rr_phase_match_fit(delta3_vec, k_RR_vec):
    """out = rr_phase_match_fit(delta3_vec, k_RR_vec): c_RR = sum(k_RR/d3) / sum(1/d3^2), R2 of the through-origin
    model, the invariant k_RR d3, its mean and spread, and the residuals."""
    d3 = np.asarray(delta3_vec, dtype=float).ravel()
    kR = np.asarray(k_RR_vec, dtype=float).ravel()
    if d3.size != kR.size:
        raise MatlabError("rr_phase_match_fit:dim", "delta3_vec and k_RR_vec must have the same length.")
    inv_d3 = 1 / d3
    c_RR = float(np.sum(kR * inv_d3) / np.sum(inv_d3 ** 2))
    y_hat = c_RR * inv_d3
    SS_res = np.sum((kR - y_hat) ** 2)
    SS_tot = np.sum((kR - np.mean(kR)) ** 2)
    inv = kR * d3
    inv_m = float(np.mean(inv))
    return {"c_RR": c_RR, "R2": float(1 - SS_res / SS_tot), "invariant": inv, "invariant_mean": inv_m,
            "invariant_spread": float((np.max(inv) - np.min(inv)) / inv_m), "residuals": kR - y_hat}


# ------------------------------------------------------------------ solvers
def _check(name, cond):
    if not cond:
        raise MatlabError("gnlse_dimensionless:param", "The value of '%s' is invalid." % name)


def _save_ids(n_steps, n_save):
    return C.mround_array(C.mlinspace(1, n_steps, n_save)).astype(int)


def gnlse_dimensionless(N_sol=3.5, delta3=0.02, xi_max=12.0, n_steps=6000, Nt=2 ** 14, tau_window=20.0, n_save=200,
                        pulse_shape="sech", mask_width=3.0, mask_order=10, verbose=True, keep_fields=True):
    """out = gnlse_dimensionless('Name', Value, ...): symmetric split-step integration of
    d_xi psi = i ((1/2) psi_tautau + |psi|^2 psi + i delta3 psi_tautautau), with D(k) = -(i/2) k^2 + i delta3 k^3,
    entropies S_hor, S_rad on the adaptive mask at n_save snapshots, eta_gsl on S_tot.  Returns a dict with the
    fields of the MATLAB struct.  keep_fields=False drops psi_hist and spec_hist (memory; numbers unchanged)."""
    _check("N_sol", N_sol > 0)
    _check("delta3", delta3 >= 0)
    _check("xi_max", xi_max > 0)
    _check("n_steps", n_steps > 0)
    _check("Nt", Nt > 0)
    _check("tau_window", tau_window > 0)
    _check("n_save", n_save > 0)
    _check("mask_width", mask_width > 0)
    _check("mask_order", mask_order > 0)
    n_steps = int(n_steps)
    n_save = int(n_save)
    Nt = int(Nt)
    TauW = float(tau_window)
    dtau = TauW / Nt
    idx = np.arange(-Nt // 2, Nt // 2, dtype=float)
    tau = idx * dtau
    k = np.fft.fftshift(idx * (2 * np.pi / TauW))                # math-natural order for D
    k_disp = np.fft.fftshift(k)                                  # monotone order for diagnostics
    shape = str(pulse_shape).lower()
    if shape == "sech":
        psi = N_sol * C.sech(tau)
    elif shape == "gaussian":
        psi = N_sol * np.exp(-tau ** 2 / 2)
    elif shape == "super_gaussian":
        psi = N_sol * np.exp(-(tau ** 2 / 2) ** 2)
    elif shape == "chirped":
        chirp = 0.5
        psi = N_sol * C.sech(tau) * np.exp(-1j * chirp * tau ** 2 / 2)
    else:
        raise MatlabError("gnlse_dimensionless:shape", 'Unknown pulse_shape "%s".' % pulse_shape)
    psi_f = np.fft.fft(psi)
    D = -1j * 0.5 * k ** 2 + 1j * delta3 * k ** 3
    save_id = _save_ids(n_steps, n_save)
    xi_vec = C.mlinspace(0, xi_max, n_save)
    psi_hist = np.zeros((Nt, n_save), dtype=complex) if keep_fields else None
    spec_hist = np.zeros((Nt, n_save)) if keep_fields else None
    S_hor = np.zeros(n_save)
    S_rad = np.zeros(n_save)
    photon_N = np.zeros(n_save)
    k_c_hist = np.zeros(n_save)
    dxi = xi_max / n_steps
    half_step = np.exp(D * dxi / 2)
    if verbose:
        print("gnlse_dimensionless: N=%.2f, delta3=%.3f, xi_max=%.2f, Nt=%d" % (N_sol, delta3, xi_max, Nt))
    fft, ifft = np.fft.fft, np.fft.ifft
    c = 0
    next_progress = 0.1
    for n in range(1, n_steps + 1):
        psi_f = psi_f * half_step
        psi_t = ifft(psi_f)
        psi_t = psi_t * np.exp(1j * np.abs(psi_t) ** 2 * dxi)
        psi_f = fft(psi_t)
        psi_f = psi_f * half_step
        if c < n_save and n == save_id[c]:
            psi_t_save = ifft(psi_f)
            P = np.abs(np.fft.fftshift(psi_f)) ** 2
            if keep_fields:
                psi_hist[:, c] = psi_t_save
                spec_hist[:, c] = P
            mS, mR, kc = soliton_mask(k_disp, P, mask_width, mask_order)
            S_hor[c] = spectral_entropy(k_disp, P, mS)
            S_rad[c] = spectral_entropy(k_disp, P, mR)
            photon_N[c] = C.mtrapz(tau, np.abs(psi_t_save) ** 2)
            k_c_hist[c] = kc
            c += 1
        if verbose and n / n_steps >= next_progress:
            print("  ... %3d%% done" % int(C.mround(100 * n / n_steps)))
            next_progress += 0.1
    S_tot = S_hor + S_rad
    eta = eta_gsl(xi_vec, S_tot)
    out = {"xi": xi_vec, "tau": tau, "k": k_disp, "psi_hist": psi_hist, "spec_hist": spec_hist, "S_hor": S_hor,
           "S_rad": S_rad, "S_tot": S_tot, "eta_gsl": eta, "Delta_S_tot": float(S_tot[-1] - S_tot[0]),
           "eta_GSL_final": float(eta[-1]), "photon_number": photon_N, "k_c": k_c_hist,
           "params": {"N_sol": N_sol, "delta3": delta3, "xi_max": xi_max, "n_steps": n_steps, "Nt": Nt,
                      "tau_window": tau_window, "n_save": n_save, "pulse_shape": pulse_shape, "mask_width": mask_width,
                      "mask_order": mask_order, "verbose": verbose}}
    if verbose:
        print("Done. Delta S_tot = %+.3f nats, eta_GSL(xi_f) = %.3f" % (out["Delta_S_tot"], out["eta_GSL_final"]))
        print("Photon-number drift: %.2e (relative)" % (photon_N[-1] / photon_N[0] - 1))
    return out


def gnlse_dimensional(beta2=-15.0, beta3=0.10, gamma=0.1, lambda0=800, N_sol=3.5, T0=0.050, n_LD=4.0, n_steps=6000,
                      Nt=2 ** 14, T_window=20.0, n_save=200, mask_width=3.0, mask_order=8, verbose=True, keep_fields=True):
    """out = gnlse_dimensional('Name', Value, ...): dA/dz = i(|beta2|/2) A_TT + (beta3/6) A_TTT + i gamma |A|^2 A,
    D(omega) = i(beta2/2 omega^2 - beta3/6 omega^3); units ps, km, W; mask width mask_width/T0 [1/ps]."""
    LD = T0 ** 2 / abs(beta2)
    P0 = (N_sol ** 2 * abs(beta2)) / (gamma * T0 ** 2)
    LNL = 1 / (gamma * P0)
    Nt = int(Nt)
    n_steps = int(n_steps)
    n_save = int(n_save)
    dt = T_window / Nt
    idx = np.arange(-Nt // 2, Nt // 2, dtype=float)
    T = idx * dt
    omega = np.fft.fftshift(idx * (2 * np.pi / T_window))
    omega_disp = np.fft.fftshift(omega)
    A = math.sqrt(P0) * C.sech(T / T0)
    A_f = np.fft.fft(A)
    D = 1j * (beta2 / 2 * omega ** 2 - beta3 / 6 * omega ** 3)
    L_prop = n_LD * LD
    dz = L_prop / n_steps
    half_step = np.exp(D * dz / 2)
    save_id = _save_ids(n_steps, n_save)
    z_vec = C.mlinspace(0, L_prop, n_save)
    A_hist = np.zeros((Nt, n_save), dtype=complex) if keep_fields else None
    spec_hist = np.zeros((Nt, n_save)) if keep_fields else None
    S_hor = np.zeros(n_save)
    S_rad = np.zeros(n_save)
    photon_N = np.zeros(n_save)
    if verbose:
        print("gnlse_dimensional: N=%.2f, T0=%.0f fs, n_LD=%.1f" % (N_sol, T0 * 1000, n_LD))
        print("  LD = %.4f km, LNL = %.4f km, P0 = %.2f kW" % (LD, LNL, P0 / 1000))
    fft, ifft = np.fft.fft, np.fft.ifft
    c = 0
    next_progress = 0.1
    for n in range(1, n_steps + 1):
        A_f = A_f * half_step
        A_t = ifft(A_f)
        A_t = A_t * np.exp(1j * gamma * np.abs(A_t) ** 2 * dz)
        A_f = fft(A_t)
        A_f = A_f * half_step
        if c < n_save and n == save_id[c]:
            A_t_save = ifft(A_f)
            P = np.abs(np.fft.fftshift(A_f)) ** 2
            if keep_fields:
                A_hist[:, c] = A_t_save
                spec_hist[:, c] = P
            mS, mR, _ = soliton_mask(omega_disp, P, mask_width / T0, mask_order)
            S_hor[c] = spectral_entropy(omega_disp, P, mS)
            S_rad[c] = spectral_entropy(omega_disp, P, mR)
            photon_N[c] = C.mtrapz(T, np.abs(A_t_save) ** 2)
            c += 1
        if verbose and n / n_steps >= next_progress:
            print("  ... %3d%% done" % int(C.mround(100 * n / n_steps)))
            next_progress += 0.1
    S_tot = S_hor + S_rad
    eta = eta_gsl(z_vec, S_tot)
    out = {"z": z_vec, "T": T, "omega": omega_disp, "nu": omega_disp / (2 * np.pi), "A_hist": A_hist, "spec_hist": spec_hist,
           "S_hor": S_hor, "S_rad": S_rad, "S_tot": S_tot, "eta_gsl": eta, "Delta_S_tot": float(S_tot[-1] - S_tot[0]),
           "eta_GSL_final": float(eta[-1]), "photon_number": photon_N, "scales": {"LD": LD, "LNL": LNL, "P0": P0},
           "params": {"beta2": beta2, "beta3": beta3, "gamma": gamma, "lambda0": lambda0, "N_sol": N_sol, "T0": T0,
                      "n_LD": n_LD, "n_steps": n_steps, "Nt": Nt, "T_window": T_window, "n_save": n_save,
                      "mask_width": mask_width, "mask_order": mask_order, "verbose": verbose}}
    if verbose:
        print("Done. Delta S_tot = %+.3f nats, eta_GSL(z_f) = %.3f" % (out["Delta_S_tot"], out["eta_GSL_final"]))
        print("Photon-number drift: %.2e (relative)" % (photon_N[-1] / photon_N[0] - 1))
    return out


# ------------------------------------------------------------------ the robustness sweep
SWEEP_GRID_N = [2.5, 3.0, 3.5, 4.0, 4.5]                 # (2.5:0.5:4.5)'
SWEEP_SHAPES = ["sech", "gaussian", "super_gaussian", "chirped"]
SWEEP_MASK_WIDTHS = [2.5, 3.0, 4.0]                      # tight, nominal, wide
SWEEP_MASK_ORDERS = [12, 10, 8]                          # sharp, medium, soft


def _sweep_single(N_sol, d3, shape, mw, mo, xi_max, Nt, n_steps):
    """run_single of run_robustness_sweep.m."""
    t0 = time.time()
    r = gnlse_dimensionless(N_sol=N_sol, delta3=d3, pulse_shape=shape, mask_width=mw, mask_order=mo, xi_max=xi_max,
                            Nt=Nt, n_steps=n_steps, verbose=False, keep_fields=False)
    return {"Delta_S_tot": r["Delta_S_tot"], "eta_GSL_final": r["eta_GSL_final"],
            "photon_drift": float(np.max(np.abs(r["photon_number"] / r["photon_number"][0] - 1))), "wall_s": time.time() - t0}


def run_robustness_sweep(fast=False, Nt=2 ** 13, n_steps=4000, xi_max=10.0, workers_map=None, verbose=True):
    """out = run_robustness_sweep('fast', ..., 'Nt', ..., 'n_steps', ..., 'xi_max', ...): the 300 configurations
    (N, delta3, pulse shape, mask scheme) in MATLAB's ndgrid linear order; DeltaS and eta as 5 x 5 x 4 x 3 arrays,
    their ranges, pass_count / pass_rate of DeltaS > 0 and eta > 0.5.  workers_map(fn, arglist) runs the
    configurations (par.pmap for process parallelism; the MATLAB 'parallel' flag used parfor).  The .mat save of
    the MATLAB function is left to the caller (stage S7 writes .npz, .json and .csv)."""
    if fast:
        Nt, n_steps, xi_max = 2 ** 12, 1500, 6.0
        if verbose:
            print("[fast mode] Nt=%d, n_steps=%d, xi_max=%.1f" % (Nt, n_steps, xi_max))
    grid_N = np.array(SWEEP_GRID_N)
    grid_d3 = C.mlinspace(0.01, 0.10, 5)
    nN, nD, nS, nM = grid_N.size, grid_d3.size, len(SWEEP_SHAPES), len(SWEEP_MASK_WIDTHS)
    Nc = nN * nD * nS * nM
    if verbose:
        print("Robustness sweep: %d configurations" % Nc)
    iN, iD, iS, iM = [a.ravel(order="F") for a in np.meshgrid(np.arange(nN), np.arange(nD), np.arange(nS), np.arange(nM), indexing="ij")]
    args = [(float(grid_N[iN[c]]), float(grid_d3[iD[c]]), SWEEP_SHAPES[iS[c]], SWEEP_MASK_WIDTHS[iM[c]], SWEEP_MASK_ORDERS[iM[c]],
             float(xi_max), int(Nt), int(n_steps)) for c in range(Nc)]
    t0 = time.time()
    if workers_map is None:
        res = []
        for c, a in enumerate(args):
            res.append(_sweep_single(*a))
            if verbose and ((c + 1) % 30 == 0 or c + 1 == Nc):
                el = time.time() - t0
                print("  %3d/%d  (%.0f s elapsed, ~%.0f s remaining)" % (c + 1, Nc, el, el / (c + 1) * (Nc - c - 1)))
    else:
        res = workers_map(_sweep_single, args)
    DeltaS = np.zeros((nN, nD, nS, nM))
    eta = np.zeros((nN, nD, nS, nM))
    drift = np.zeros((nN, nD, nS, nM))
    for c in range(Nc):
        DeltaS[iN[c], iD[c], iS[c], iM[c]] = res[c]["Delta_S_tot"]
        eta[iN[c], iD[c], iS[c], iM[c]] = res[c]["eta_GSL_final"]
        drift[iN[c], iD[c], iS[c], iM[c]] = res[c]["photon_drift"]
    passed = (DeltaS.ravel(order="F") > 0) & (eta.ravel(order="F") > 0.5)
    out = {"grid_N": grid_N, "grid_d3": grid_d3, "shapes": list(SWEEP_SHAPES), "mask_widths": np.array(SWEEP_MASK_WIDTHS),
           "mask_orders": np.array(SWEEP_MASK_ORDERS), "DeltaS": DeltaS, "eta": eta,
           "DeltaS_range": [float(np.min(DeltaS)), float(np.max(DeltaS))], "eta_range": [float(np.min(eta)), float(np.max(eta))],
           "pass_count": int(np.sum(passed)), "pass_rate": float(np.mean(passed)), "runtime_s": time.time() - t0,
           "settings": {"fast": bool(fast), "Nt": int(Nt), "n_steps": int(n_steps), "xi_max": float(xi_max), "n_save": 200,
                        "tau_window": 20.0}, "photon_drift_max": drift,
           "config_order": "MATLAB ndgrid linear order: N fastest, then delta3, pulse shape, mask scheme",
           "rows": [{"N_sol": a[0], "delta3": a[1], "pulse_shape": a[2], "mask_width": a[3], "mask_order": a[4],
                     "DeltaS_tot": res[c]["Delta_S_tot"], "eta_GSL_final": res[c]["eta_GSL_final"],
                     "photon_drift_max": res[c]["photon_drift"], "wall_s": res[c]["wall_s"],
                     "pass": bool(res[c]["Delta_S_tot"] > 0 and res[c]["eta_GSL_final"] > 0.5)} for c, a in enumerate(args)]}
    if verbose:
        print("\n=== Robustness sweep summary ===")
        print("  Delta S_tot range: [%.3f, %.3f] nats   (paper: [0.82, 2.42])" % tuple(out["DeltaS_range"]))
        print("  eta_GSL range:     [%.3f, %.3f]         (paper: [0.52, 0.62])" % tuple(out["eta_range"]))
        print("  Pass rate (DeltaS>0 AND eta>0.5): %d/%d = %.1f%%" % (out["pass_count"], Nc, 100 * out["pass_rate"]))
    return out


# ------------------------------------------------------------------ checks of the other Ref. [25] functions (stage S7)
PM_DELTA3 = [0.06, 0.07, 0.08, 0.09, 0.10]


def _nominal_check(fast):
    """Test 1 of SOTHE_pkg validation/validate_paper_claims.py, through the ports of the MATLAB functions."""
    kw = dict(Nt=2 ** 12, n_steps=1500, n_save=80, xi_max=6.0) if fast else \
        dict(Nt=2 ** 13, n_steps=3000, n_save=100, xi_max=10.0, tau_window=30.0)
    r = gnlse_dimensionless(N_sol=3.5, delta3=0.02, pulse_shape="sech", mask_width=3.0, mask_order=10, verbose=False,
                            keep_fields=False, **kw)
    neg = float(np.mean(np.diff(r["S_tot"]) < 0))
    drift = float(abs(r["photon_number"][-1] / r["photon_number"][0] - 1))
    res = {"settings": kw, "Delta_S_tot": r["Delta_S_tot"], "eta_GSL": r["eta_GSL_final"], "neg_frac": neg, "drift": drift}
    res["passed"] = bool(r["Delta_S_tot"] > 0 and r["eta_GSL_final"] > 0.5 and 0.30 < neg < 0.65 and drift < 1e-3)
    return res


def _pm_single(d3, fast):
    kw = dict(Nt=2 ** 12, n_steps=1500, n_save=60, xi_max=6.0) if fast else \
        dict(Nt=2 ** 13, n_steps=3000, n_save=60, xi_max=10.0, tau_window=30.0)
    r = gnlse_dimensionless(N_sol=3.5, delta3=d3, verbose=False, **kw)
    P_end = r["spec_hist"][:, -1]
    k = r["k"]
    _, _, k_c = soliton_mask(k, P_end, r["params"]["mask_width"], r["params"]["mask_order"])
    return float(find_rr_lobe(k, P_end, k_c, r["params"]["mask_width"]))


def _dimensional_check(fast):
    kw = dict(Nt=2 ** 12, n_steps=1500) if fast else {}
    r = gnlse_dimensional(verbose=False, keep_fields=False, **kw)
    return {"settings": dict(r["params"]), "Delta_S_tot": r["Delta_S_tot"], "eta_GSL": r["eta_GSL_final"],
            "photon_drift": float(np.max(np.abs(r["photon_number"] / r["photon_number"][0] - 1))), "scales": r["scales"]}


def _check_task(name, arg):
    if name == "nominal":
        return _nominal_check(arg)
    if name == "dimensional":
        return _dimensional_check(arg)
    d3, fast = arg
    return _pm_single(d3, fast)


def ref25_checks(fast=False, workers_map=None):
    """The functions no sweep uses, exercised as SOTHE_pkg's own validation script uses them (validate_paper_claims.py
    Test 1 and Test 2, same settings and pass rules) plus one gnlse_dimensional run at the defaults of the MATLAB file
    (the physical parameters of Ref. [25] Sec. 3.1).  Returns a dict; the pass flags are information, not stage criteria."""
    jobs = [("nominal", fast), ("dimensional", fast)] + [("pm", (d3, fast)) for d3 in PM_DELTA3]
    res = workers_map(_check_task, jobs) if workers_map else [_check_task(*j) for j in jobs]
    nominal, dim, kRR = res[0], res[1], np.array(res[2:], dtype=float)
    fit = rr_phase_match_fit(PM_DELTA3, kRR)
    pm = {"delta3": list(PM_DELTA3), "k_RR": kRR, "k_RR_times_delta3": kRR * np.array(PM_DELTA3), "c_RR": fit["c_RR"],
          "R2": fit["R2"], "invariant_mean": fit["invariant_mean"], "invariant_spread": fit["invariant_spread"],
          "claims_of_ref25": {"c_RR": 1.01, "R2": 0.93, "invariant_mean": 1.02, "invariant_spread_below": 0.05}}
    pm["passed"] = bool(0.8 < fit["invariant_mean"] < 1.2 and fit["R2"] > 0.5)
    return {"protocol": "SOTHE_pkg v1.0.0 validation/validate_paper_claims.py (Tests 1 and 2) with the ports of the MATLAB "
                        "functions (find_rr_lobe with its own reflect padding and 3-sigma Gaussian); gnlse_dimensional at the "
                        "defaults of gnlse_dimensional.m", "fast": bool(fast),
            "nominal": nominal, "phase_matching": pm, "dimensional": dim}
