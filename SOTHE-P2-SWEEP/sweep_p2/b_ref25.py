"""B8: variants of the Ref. [25] robustness sweep (FINDING S7; REVIEW_G5 W2; OD-P8).

The shipped script (run_robustness_sweep.m of SOTHE_pkg v1.0.0, ported statement by statement in
sothe_p2/ref25/src.py) runs 5 N x 5 delta3 x 4 pulses x 3 masks at N_t = 2^13, 4000 steps, xi = 10.  The version of
record (Class. Quantum Grav. 43, 135014) states the production resolution N_t = 2^14 (6000 steps), a super-Gaussian
N exp(-tau^8/2) ("order m = 4"; the script has N exp(-tau^4/4)) and three masking schemes whose code is not in the
package.  The variants (grid.ref25_variant) separate these differences one at a time.

gnlse_ext() is src.gnlse_dimensionless() without the field history, with one more input pulse, 'super_gaussian_m4'
= N exp(-tau^8/2).  For the four pulses of the script it performs the same floating-point operations in the same order
as src.gnlse_dimensionless (oracle O-B8-ext-identity: bitwise equal DeltaS, eta and photon numbers)."""
import time

import numpy as np

from sothe_p2 import compat as C
from sothe_p2.ref25 import src

SRC_SHAPES = ("sech", "gaussian", "super_gaussian", "chirped")
RANGES_PAPER1 = {"DeltaS": (0.82, 2.42), "eta": (0.52, 0.62)}


def _initial(shape, tau, N_sol):
    s = str(shape).lower()
    if s == "sech":
        return N_sol * C.sech(tau)
    if s == "gaussian":
        return N_sol * np.exp(-tau ** 2 / 2)
    if s == "super_gaussian":
        return N_sol * np.exp(-(tau ** 2 / 2) ** 2)
    if s == "chirped":
        chirp = 0.5
        return N_sol * C.sech(tau) * np.exp(-1j * chirp * tau ** 2 / 2)
    if s == "super_gaussian_m4":
        return N_sol * np.exp(-tau ** 8 / 2)
    raise ValueError("unknown pulse shape %s" % shape)


def gnlse_ext(N_sol, delta3, xi_max, n_steps, Nt, pulse_shape, mask_width, mask_order, tau_window=20.0, n_save=200):
    n_steps, n_save, Nt = int(n_steps), int(n_save), int(Nt)
    TauW = float(tau_window)
    dtau = TauW / Nt
    idx = np.arange(-Nt // 2, Nt // 2, dtype=float)
    tau = idx * dtau
    k = np.fft.fftshift(idx * (2 * np.pi / TauW))
    k_disp = np.fft.fftshift(k)
    psi = _initial(pulse_shape, tau, N_sol)
    psi_f = np.fft.fft(psi)
    D = -1j * 0.5 * k ** 2 + 1j * delta3 * k ** 3
    save_id = src._save_ids(n_steps, n_save)
    xi_vec = C.mlinspace(0, xi_max, n_save)
    S_hor = np.zeros(n_save)
    S_rad = np.zeros(n_save)
    photon_N = np.zeros(n_save)
    dxi = xi_max / n_steps
    half_step = np.exp(D * dxi / 2)
    fft, ifft = np.fft.fft, np.fft.ifft
    c = 0
    for n in range(1, n_steps + 1):
        psi_f = psi_f * half_step
        psi_t = ifft(psi_f)
        psi_t = psi_t * np.exp(1j * np.abs(psi_t) ** 2 * dxi)
        psi_f = fft(psi_t)
        psi_f = psi_f * half_step
        if c < n_save and n == save_id[c]:
            psi_t_save = ifft(psi_f)
            P = np.abs(np.fft.fftshift(psi_f)) ** 2
            mS, mR, kc = src.soliton_mask(k_disp, P, mask_width, mask_order)
            S_hor[c] = src.spectral_entropy(k_disp, P, mS)
            S_rad[c] = src.spectral_entropy(k_disp, P, mR)
            photon_N[c] = C.mtrapz(tau, np.abs(psi_t_save) ** 2)
            c += 1
    S_tot = S_hor + S_rad
    eta = src.eta_gsl(xi_vec, S_tot)
    return {"Delta_S_tot": float(S_tot[-1] - S_tot[0]), "eta_GSL_final": float(eta[-1]), "photon_number": photon_N}


def sweep_single(N_sol, d3, shape, mw, mo, xi_max, Nt, n_steps):
    """One configuration: src._sweep_single for the script's pulses, gnlse_ext for the VoR super-Gaussian."""
    if shape in SRC_SHAPES:
        return src._sweep_single(N_sol, d3, shape, mw, mo, xi_max, Nt, n_steps)
    t0 = time.time()
    r = gnlse_ext(N_sol, d3, xi_max, n_steps, Nt, shape, mw, mo)
    return {"Delta_S_tot": r["Delta_S_tot"], "eta_GSL_final": r["eta_GSL_final"],
            "photon_drift": float(np.max(np.abs(r["photon_number"] / r["photon_number"][0] - 1))), "wall_s": time.time() - t0}


def configs(v):
    """The configurations of a variant in MATLAB ndgrid linear order (N fastest, then delta3, pulse, mask)."""
    out = []
    for im, (mw, mo) in enumerate(v["masks"]):
        for isx, shape in enumerate(v["shapes"]):
            for iD, d3 in enumerate(v["d3"]):
                for iN, N in enumerate(v["N"]):
                    out.append({"iN": iN, "iD": iD, "iS": isx, "iM": im, "N_sol": float(N), "delta3": float(d3), "pulse_shape": shape,
                                "mask_width": float(mw), "mask_order": int(mo)})
    return out


def chunk(v, iN, iS, iM):
    """The delta3 column of one (N, pulse, mask): the unit of work of a B8 task."""
    rows = []
    for c in configs(v):
        if c["iN"] == iN and c["iS"] == iS and c["iM"] == iM:
            r = sweep_single(c["N_sol"], c["delta3"], c["pulse_shape"], c["mask_width"], c["mask_order"], v["xi_max"], v["Nt"], v["n_steps"])
            c.update(DeltaS_tot=r["Delta_S_tot"], eta_GSL_final=r["eta_GSL_final"], photon_drift_max=r["photon_drift"], wall_s=r["wall_s"],
                     passed=bool(r["Delta_S_tot"] > 0 and r["eta_GSL_final"] > 0.5))
            rows.append(c)
    return rows


def summarize(v, rows):
    """Ranges, pass counts, counts outside the ranges of Paper 1, subset ranges, extremes."""
    rows = sorted(rows, key=lambda r: (r["iM"], r["iS"], r["iD"], r["iN"]))
    dS = np.array([r["DeltaS_tot"] for r in rows])
    et = np.array([r["eta_GSL_final"] for r in rows])
    lo, hi = RANGES_PAPER1["DeltaS"]
    elo, ehi = RANGES_PAPER1["eta"]
    out = {"variant": v["name"], "note": v["note"], "Nt": v["Nt"], "n_steps": v["n_steps"], "xi_max": v["xi_max"],
           "N": v["N"], "d3": v["d3"], "shapes": v["shapes"], "masks": v["masks"], "n_configs": len(rows),
           "DeltaS_range": [float(dS.min()), float(dS.max())], "eta_range": [float(et.min()), float(et.max())],
           "pass_count": int(sum(r["passed"] for r in rows)),
           "n_DeltaS_below": int(np.sum(dS < lo)), "n_DeltaS_above": int(np.sum(dS > hi)),
           "n_eta_below": int(np.sum(et < elo)), "n_eta_above": int(np.sum(et > ehi)),
           "photon_drift_max": float(max(r["photon_drift_max"] for r in rows)), "wall_s_sum": float(sum(r["wall_s"] for r in rows))}
    sub = {}
    for key, vals in (("pulse", v["shapes"]), ("mask", ["%g/%d" % (a, b) for a, b in v["masks"]])):
        for i, name in enumerate(vals):
            sel = [r for r in rows if (r["iS"] if key == "pulse" else r["iM"]) == i]
            sub["%s=%s" % (key, name)] = {"DeltaS": [min(r["DeltaS_tot"] for r in sel), max(r["DeltaS_tot"] for r in sel)],
                                          "eta": [min(r["eta_GSL_final"] for r in sel), max(r["eta_GSL_final"] for r in sel)]}
    out["subsets"] = sub
    imin, imax = int(np.argmin(dS)), int(np.argmax(dS))
    out["extremes"] = {"DeltaS_min": {k: rows[imin][k] for k in ("N_sol", "delta3", "pulse_shape", "mask_width", "mask_order", "DeltaS_tot", "eta_GSL_final")},
                       "DeltaS_max": {k: rows[imax][k] for k in ("N_sol", "delta3", "pulse_shape", "mask_width", "mask_order", "DeltaS_tot", "eta_GSL_final")}}
    return out, rows


def compare_reference(rows, ref_csv):
    """V1 against job 530047 (stage S7 of SOTHE-P2): row by row, same order."""
    with open(ref_csv, encoding="utf-8") as f:
        lines = [l.strip().split(",") for l in f if l.strip()]
    head, body = lines[0], lines[1:]
    iS, iE = head.index("DeltaS_tot"), head.index("eta_GSL_final")
    ref = {(float(b[0]), float(b[1]), b[2], float(b[3]), int(b[4])): (float(b[iS]), float(b[iE])) for b in body}
    dmax, emax, n, tight, worst = 0.0, 0.0, 0, 0, None
    for r in rows:
        key = (r["N_sol"], r["delta3"], r["pulse_shape"], r["mask_width"], r["mask_order"])
        if key in ref:
            n += 1
            dS, de = abs(r["DeltaS_tot"] - ref[key][0]), abs(r["eta_GSL_final"] - ref[key][1])
            tight += (dS <= 1e-6 and de <= 1e-6)
            if dS > dmax:
                worst = list(key)
            dmax, emax = max(dmax, dS), max(emax, de)
    return {"n_matched": n, "n_reference": len(ref), "DeltaS_max_abs_dev": dmax, "eta_max_abs_dev": emax, "n_within_1e-6": tight,
            "worst_config": worst}


def ext_identity(Nt=2 ** 10, n_steps=200, xi_max=2.0):
    """gnlse_ext against src.gnlse_dimensionless for the four pulses of the script (small grid): bitwise equality."""
    res = []
    for shape in SRC_SHAPES:
        a = src.gnlse_dimensionless(N_sol=3.0, delta3=0.05, xi_max=xi_max, n_steps=n_steps, Nt=Nt, pulse_shape=shape,
                                    mask_width=3.0, mask_order=10, verbose=False, keep_fields=False)
        b = gnlse_ext(3.0, 0.05, xi_max, n_steps, Nt, shape, 3.0, 10)
        res.append({"shape": shape, "DeltaS_equal": a["Delta_S_tot"] == b["Delta_S_tot"], "eta_equal": a["eta_GSL_final"] == b["eta_GSL_final"],
                    "photon_equal": bool(np.array_equal(a["photon_number"], b["photon_number"])),
                    "DeltaS": a["Delta_S_tot"], "eta": a["eta_GSL_final"]})
    return {"Nt": Nt, "n_steps": n_steps, "xi_max": xi_max, "rows": res,
            "all_equal": all(r["DeltaS_equal"] and r["eta_equal"] and r["photon_equal"] for r in res)}
