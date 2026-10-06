"""B5: BdG spectrum about the stationary TOD soliton phi (kit.tod_soliton + kit.bdg_spectrum_bg of SOTHE-P2 1.0.0),
with the zero-sector rank tests of REVIEW_G5 / A-2 (singular values of M and M^2) and the psi0-linearization artifact
(suite.bdg_spectrum about psi0, excluded by the manuscript): psi0_annulus_maxIm is the quantity of the v0.1.0 rate table
(max |Im| for 0.1 mu < |Re w| <= 3 mu); psi0_quartet_* are the four eigenvalues of smallest modulus.
Scope: the manuscript keeps the BdG blocks at N delta3 <= 0.05; larger N delta3 is reported, not claimed."""
import math

import numpy as np

from sothe_p2 import kit, suite


def _M(phi, N, d3, n_tau, tau_max):
    h = 2.0 * tau_max / n_tau
    D2 = suite.fd_circulant(n_tau, h, 4, 2)
    D3 = suite.fd_circulant(n_tau, h, 4, 3)
    mu = 0.5 * N ** 2
    H = -0.5 * D2 - 2.0 * np.diag(np.abs(phi) ** 2) + mu * np.eye(n_tau)
    T = (-1j * d3) * D3
    return np.block([[H + T, -np.diag(phi ** 2)], [np.diag(np.conj(phi) ** 2), -H + T]])


def _kernel_dim(sv_rel, floor=1e-9):
    """Number of singular values below the largest multiplicative gap among the smallest ones that crosses `floor`."""
    s = np.sort(np.asarray(sv_rel))[:12]
    small = int(np.sum(s < floor))
    return small


def bdg_point(N, d3, n_tau, tau_max):
    S = kit.tod_soliton(N, d3, n_tau, tau_max)
    phi, tau = S["phi"], S["tau"]
    B = kit.bdg_spectrum_bg(phi, N, d3, n_tau, tau_max)
    w, kr, ann = B["w"], B["krein"], B["ann"]
    mu = 0.5 * N ** 2
    rho = float(np.max(np.abs(w)))
    eps = float(np.finfo(float).eps)
    M = _M(phi, N, d3, n_tau, tau_max)
    s1 = np.linalg.svd(M, compute_uv=False)
    s2 = np.linalg.svd(M @ M, compute_uv=False)
    sv1 = (s1[::-1] / s1[0])[:8]
    sv2 = (s2[::-1] / s2[0])[:8]
    out = {"N": N, "d3": d3, "n_tau": n_tau, "tau_max": tau_max, "mu": mu,
           "newton_converged": bool(S["converged"]), "newton_res": float(S["res"][-1]), "newton_iters": len(S["res"]),
           "kmean": float(S["kmean"]), "tail_rel": float(np.max(np.abs(phi[np.abs(tau) > 0.8 * tau_max])) / N),
           "annulus_maxIm": B["maxIm"], "gap_over_mu": B["gap"] / mu, "n_annulus": int(np.sum(ann)),
           "krein_pos": int(np.sum(kr[ann] > 0)), "krein_neg": int(np.sum(kr[ann] < 0)),
           "zero_sector": [complex(z) for z in B["w_small"]], "zero_sector_max": float(np.max(np.abs(B["w_small"]))),
           "rho": rho, "sqrt_eps_rho": math.sqrt(eps * rho), "sv_M_rel_smallest": sv1.tolist(), "sv_M2_rel_smallest": sv2.tolist(),
           "dim_ker_M": _kernel_dim(sv1, 1e-12), "dim_ker_M2": _kernel_dim(sv2, 1e-12),
           "sv_M_gap": float(sv1[2] / max(sv1[1], 1e-300)), "sv_M2_gap": float(sv2[4] / max(sv2[3], 1e-300))}
    if d3 > 0:
        P = suite.bdg_spectrum(N, d3, n_tau, tau_max)
        w0 = P["w"]
        o = np.argsort(np.abs(w0), kind="stable")[:4]
        q = float(np.max(np.abs(w0[o].imag)))
        out.update(psi0_quartet=[complex(z) for z in w0[o]], psi0_quartet_maxIm=q, psi0_quartet_over_sqrt_d3=q / math.sqrt(d3),
                   psi0_annulus_maxIm=P["maxIm"], psi0_annulus_over_sqrt_d3=P["maxIm"] / math.sqrt(d3))
    spectra = {"w": w, "krein": kr, "phi": phi, "tau": tau}
    return out, spectra
