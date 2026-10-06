"""Symmetric split-step integration of the GNLSE with third-order dispersion, with the diagnostics of the
steady-rate and recoil analyses.

    d_xi psi = i ( 1/2 psi_tt + |psi|^2 psi + i delta3(xi) psi_ttt ) + V psi_tau

V is the velocity of the computational frame (tau' = tau - V xi); V = 0 is the laboratory frame.  One step is
half a linear step, the Kerr step, half a linear step, then the absorber exp(-sigma dxi) with
sigma = s0 ((|tau| - a)/(L/2 - a))^2 for |tau| > a, a = 0.35 L, s0 = 8 (the scheme, grid and absorber of the
radiative-loss runs of SOTHE-P2 1.0.0, kit.splitstep_loss: with V = 0, a constant delta3 and the default grid
the two programs perform the same arithmetic).  delta3(xi) is constant, or rises as
delta3 (1 - cos(pi xi / xi_ramp)) / 2 over [0, xi_ramp] and is constant afterwards ('ramp'); within a step it is
evaluated at the midpoint of the step.

Diagnostics every `sample` units of xi (the field itself is not stored):
  Q        norm of the core window |tau - tau[ip]| < 6/N about the grid maximum ip (as SOTHE-P2)
  Qbox     norm of the region without absorption, |tau| < a
  tc       peak position: vertex of the parabola through ln|psi|^2 at ip and its two neighbours (frame
           coordinate; add V xi for the laboratory frame)
  amp      peak amplitude from the same parabola
  ks_hann  spectral centroid of the Hann-windowed core (SOTHE-P2 definition)
  ks_mom   momentum centroid of the core: int_core Im(psi* psi_tau) / Q
  Mbox     momentum of the region without absorption
  phase    arg psi at the peak (unwrapped later); d(phase)/dxi is the phase rotation rate of the soliton
  E        Hamiltonian sum_k beta0(k)|psi_k|^2 - (1/2) int |psi|^4 over the whole grid (diagnostic)
  Iband    mean |b|^2 over the flat part of a tail window [tc - d1, tc - d0] (behind the soliton, where the
           resonant radiation goes), b the field filtered to |k - k_c| < bw, k_c the predicted resonant wavenumber
  Iband_near, Iband_far   the same over the nearer and the farther half of the flat part
  kband    power-weighted mean wavenumber inside the band
The emitted norm flux follows in the analysis as F = Iband |v_g0(kband) - V_soliton| (radiation of a steady
emission is a plane wave of uniform intensity behind the soliton).
"""
import math
import time

import numpy as np

try:                                   # SciPy's pocketfft is faster here; NumPy's is the fallback
    from scipy import fft as _fft
except ImportError:                    # pragma: no cover
    _fft = np.fft

from . import model


DEFAULTS = {
    "N": 1.0, "d3": 0.05, "schedule": "const", "xi_ramp": 0.0, "xi_end": 40.0, "dxi": 0.002, "n": 8192, "L": 400.0,
    "V": 0.0, "absorb_frac": 0.35, "absorb_strength": 8.0, "sample": 0.25, "kfilter": 0.0,
    "tail_d0": 20.0, "tail_d1": 120.0, "band_halfwidth": 2.0, "core_halfwidth": 6.0, "seed": None,
}


def cfg_full(cfg):
    c = dict(DEFAULTS)
    unknown = set(cfg) - set(DEFAULTS) - {"id", "group", "note"}
    if unknown:
        raise ValueError("unknown configuration keys: %s" % ", ".join(sorted(unknown)))
    c.update(cfg)
    if c["schedule"] not in ("const", "ramp"):
        raise ValueError("schedule must be 'const' or 'ramp'")
    if c["schedule"] == "ramp" and not c["xi_ramp"] > 0:
        raise ValueError("a ramp needs xi_ramp > 0")
    return c


def d3_at(c, xi):
    if c["schedule"] == "ramp" and xi < c["xi_ramp"]:
        return c["d3"] * 0.5 * (1.0 - math.cos(math.pi * xi / c["xi_ramp"]))
    return c["d3"]


def grid(c):
    n, L = int(c["n"]), float(c["L"])
    h = L / n
    tau = -L / 2 + h * np.arange(n)
    k = (2 * np.pi / L) * np.concatenate([np.arange(0, n // 2), np.arange(-n // 2, 0)])
    a = c["absorb_frac"] * L
    sig = np.zeros(n)
    m = np.abs(tau) > a
    sig[m] = c["absorb_strength"] * ((np.abs(tau[m]) - a) / (L / 2 - a)) ** 2
    return tau, k, h, sig, a


def subgrid(I2, ip, tau, h):
    """Vertex (position and value) of the parabola through ln I2 at ip and its two periodic neighbours."""
    n = I2.size
    im, i0, iq = I2[(ip - 1) % n], I2[ip], I2[(ip + 1) % n]
    if im <= 0 or i0 <= 0 or iq <= 0:
        return float(tau[ip]), float(i0)
    lm, l0, lq = math.log(im), math.log(i0), math.log(iq)
    den = lm - 2.0 * l0 + lq
    if not den < 0:
        return float(tau[ip]), float(i0)
    off = 0.5 * (lm - lq) / den
    lmax = l0 - (lq - lm) ** 2 / (8.0 * den)
    return float(tau[ip] + off * h), math.exp(lmax)


def tukey(m, frac=0.2):
    w = np.ones(m)
    t = int(frac * m / 2)
    if t > 0:
        r = 0.5 * (1 - np.cos(np.pi * np.arange(t) / t))
        w[:t] = r
        w[-t:] = r[::-1]
    return w, t


def measure(psi, c, tau, k, h, a, d3_now, fft=_fft.fft, ifft=_fft.ifft):
    N = c["N"]
    n = psi.size
    I2 = np.abs(psi) ** 2
    ip = int(np.argmax(I2))
    tgrid = tau[ip]
    core = np.abs(tau - tgrid) < c["core_halfwidth"] / N
    idx = np.flatnonzero(core)
    Q = float(np.sum(I2[core]) * h)
    box = np.abs(tau) < a
    Qbox = float(np.sum(I2[box]) * h)
    tc, Imax = subgrid(I2, ip, tau, h)
    nc = idx.size
    w = np.zeros(n)
    w[idx] = 0.5 * (1 - np.cos(2 * np.pi * np.arange(nc) / (nc - 1)))
    Pk = np.abs(fft(psi * w)) ** 2
    ks_hann = float(np.sum(k * Pk) / np.sum(Pk))
    dpsi = ifft(1j * k * fft(psi))
    P = np.imag(np.conj(psi) * dpsi)
    ks_mom = float(np.sum(P[core]) * h / Q)
    Mbox = float(np.sum(P[box]) * h)
    ps = fft(psi)
    E = float(np.sum(model.beta0(k, d3_now) * np.abs(ps) ** 2) * h / n - 0.5 * np.sum(I2 ** 2) * h)
    phase = float(np.angle(psi[ip])) + ks_mom * (tc - float(tgrid))      # phase at the sub-grid centre
    # tail window behind the soliton
    lo, hi = tc - c["tail_d1"], tc - c["tail_d0"]
    lo, hi = max(lo, -a + 5.0), min(hi, a - 5.0)
    out = {"Q": Q, "Qbox": Qbox, "tc": tc, "tgrid": float(tgrid), "amp": math.sqrt(Imax), "ks_hann": ks_hann,
           "ks_mom": ks_mom, "Mbox": Mbox, "phase": phase, "E": E, "Iband": float("nan"), "Iband_near": float("nan"),
           "Iband_far": float("nan"), "kband": float("nan"), "kc": float("nan"),
           "prof": [float("nan")] * NBINS, "dist": [float("nan")] * NBINS}
    if hi - lo > 20.0 and d3_now > 0:
        sel = np.flatnonzero((tau >= lo) & (tau <= hi))
        m = sel.size
        tw, t = tukey(m, 0.2)
        ww = np.zeros(n)
        ww[sel] = tw
        Nloc = model.amplitude_from_norm(Q, ks_mom, d3_now)
        try:
            kc = model.k_res_lab(Nloc, ks_mom, d3_now, "drift")
        except Exception:
            kc = float("nan")
        if math.isfinite(kc):
            spec = fft(psi * ww)
            band = np.abs(k - kc) < c["band_halfwidth"]
            bspec = np.where(band, spec, 0)
            b = ifft(bspec)
            flat = sel[t:m - t] if m - 2 * t > 4 else sel
            Ib = np.abs(b[flat]) ** 2
            half = flat.size // 2
            pw = np.abs(bspec) ** 2
            prof = [float(np.mean(x)) for x in np.array_split(Ib, NBINS)]          # far ... near
            dist = [float(tc - np.mean(tau[x])) for x in np.array_split(flat, NBINS)]
            out.update(prof=prof, dist=dist)
            out.update(Iband=float(np.mean(Ib)), Iband_far=float(np.mean(Ib[:half])), Iband_near=float(np.mean(Ib[half:])),
                       kband=float(np.sum(k[band] * pw[band]) / np.sum(pw[band])) if np.sum(pw[band]) > 0 else float("nan"),
                       kc=float(kc))
    return out


KEYS = ["xi", "d3", "Q", "Qbox", "tc", "tgrid", "amp", "ks_hann", "ks_mom", "Mbox", "phase", "E", "Iband", "Iband_near",
        "Iband_far", "kband", "kc"]
NBINS = 8


def integrate(cfg, progress=False, keep_field=False):
    """Run one configuration; returns (series: dict of arrays, info: dict)."""
    c = cfg_full(cfg)
    tau, k, h, sig, a = grid(c)
    N = c["N"]
    dxi = c["dxi"]
    psi = (N / np.cosh(N * tau)).astype(complex)
    absorb = np.exp(-sig * dxi)
    kf = None
    if c["kfilter"] and c["kfilter"] > 0:
        kmax = np.max(np.abs(k))
        kf = np.exp(-(np.abs(k) / (c["kfilter"] * kmax)) ** 36)
    base = -0.5 * k ** 2 + c["V"] * k
    k3 = k ** 3
    ns = int(round(c["xi_end"] / dxi))
    every = int(round(c["sample"] / dxi))
    if abs(every * dxi - c["sample"]) > 1e-12 or abs(ns * dxi - c["xi_end"]) > 1e-9:
        raise ValueError("sample and xi_end must be multiples of dxi")
    fft, ifft = _fft.fft, _fft.ifft
    rows, profs, dists = [], [], []
    rot = np.empty(psi.size, dtype=complex)
    lin, lin_d3 = None, None
    t0 = time.time()
    for s in range(ns + 1):
        xi = s * dxi
        if s % every == 0:
            d3n = d3_at(c, xi)
            mrow = measure(psi, c, tau, k, h, a, d3n)
            mrow["xi"] = xi
            mrow["d3"] = d3n
            rows.append([mrow[key] for key in KEYS])
            profs.append(mrow["prof"])
            dists.append(mrow["dist"])
            if progress and (s // every) % 200 == 0:
                print("  xi %.1f  Q %.10f  tc %.4f  ks %.5f  (%.0f s)" % (xi, mrow["Q"], mrow["tc"], mrow["ks_mom"], time.time() - t0), flush=True)
        if s == ns:
            break
        d3m = d3_at(c, xi + 0.5 * dxi)
        if lin is None or d3m != lin_d3:
            lin = np.exp(1j * (base + d3m * k3) * dxi / 2)
            lin_d3 = d3m
        psi = ifft(lin * fft(psi))
        ph = (psi.real * psi.real + psi.imag * psi.imag) * dxi
        rot.real = np.cos(ph)
        rot.imag = np.sin(ph)
        psi = psi * rot
        psi = ifft(lin * fft(psi))
        psi = psi * absorb
        if kf is not None:
            psi = ifft(kf * fft(psi))
    arr = np.array(rows)
    series = {key: arr[:, i] for i, key in enumerate(KEYS)}
    series["Iband_prof"] = np.array(profs)
    series["Iband_dist"] = np.array(dists)
    info = {"cfg": c, "wall_s": time.time() - t0, "steps": ns, "h": h, "kmax": float(np.max(np.abs(k))), "absorb_start": a}
    if keep_field:
        info["tau"] = tau
        info["psi_end"] = psi
    return series, info
