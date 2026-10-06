"""The sweep: grids, protocols and block list (full run and --fast smoke test).

Every protocol value that SOTHE-P2 1.0.0 fixes is kept (kinematic grid 4096 points on [-20, 20], split-step loss
dxi = 0.002, n = 8192, L = 400, xi = 40, BdG n_tau = 500 on [-16, 16], band omega in [0.05, 1.6] with 60 points,
census 311 points).  The sweep adds axes, never changes a protocol."""
from sothe_p2 import compat as C

BLOCKS = {
    "B1": "kinematic layer over (N, k_op, delta3, kappa_g): closed form, extraction, drive scan, kappa_g flow",
    "B2": "entanglement layer (closed form) and the Hawking window, per grid point (reduction)",
    "B3": "radiative loss and drift of the launched soliton (split step, absorbing edges, convergence runs)",
    "B4": "comoving channel census with the measured drift epochs; Hawking window and escape threshold",
    "B5": "BdG about the stationary TOD soliton phi: spectrum, Krein census, zero-sector rank tests, psi0 quartet",
    "B6": "Raman and self-steepening extension (RK4IP GNLSE): SSFS against Gordon, drift, census at the operating points",
    "B7": "mode-conversion scattering with the kinematic flow (delta3 = 0): S-matrix, n_H against thermal (D-06; not a Paper-2 claim)",
    "B8": "Ref. [25] robustness sweep variants: shipped script, production resolution, VoR super-Gaussian, draft delta3 grid, xi = 12",
    "B9": "extraction-bias parity tables (Supplement Tables 2 and 3)",
    "B10": "anchor products at (N, k_op) = (1, 1.5): paper2-style CSV files, summary, figures",
}
ALL_BLOCKS = ("B1", "B2", "B3", "B4", "B5", "B6", "B7", "B8", "B9", "B10")
DEPS = {"B2": ("B1", "B4"), "B4": ("B3",), "B10": ("B1", "B2", "B3", "B4", "B5", "B7")}


def expand_blocks(blocks):
    """The requested blocks with the blocks they read (DEPS), in the order of ALL_BLOCKS."""
    want = set()
    stack = list(blocks)
    while stack:
        b = stack.pop()
        if b not in ALL_BLOCKS:
            raise ValueError("unknown block %s (B1..B10)" % b)
        if b not in want:
            want.add(b)
            stack.extend(DEPS.get(b, ()))
    return [b for b in ALL_BLOCKS if b in want]

RECORDED_LOSS_PTS = [(1.0, 0.02), (1.0, 0.03), (1.0, 0.04), (1.0, 0.05), (1.0, 0.06), (1.0, 0.08), (1.0, 0.10),
                     (2.0, 0.025), (2.0, 0.05)]                      # job 530047, SOTHE-P2 S2 (g4c.PTS)
OP = {"N": 2.0, "d3": 0.05, "kg": 0.3, "k_op": 1.0}                 # Eq. (6) of the manuscript
ANCHOR = {"N": 1.0, "d3": 0.05, "kg": 0.3, "k_op": 1.5}              # REVIEW_G5 S1-1 option A
FINDINGS_OMEGA = [0.2, 0.6, 1.0, 1.6]                                 # FINDINGS.md table (N = 2, k_op = 3, delta3 = 0)
FINDINGS_RATIO = [0.47966, 0.48300, 0.49122, 0.50422]


def _u(xs):
    out = []
    for x in xs:
        if all(abs(x[0] - y[0]) > 1e-12 or abs(x[1] - y[1]) > 1e-12 for y in out):
            out.append(x)
    return out


def config(fast=False):
    """The full sweep (fast=False) or the smoke test (fast=True).  Plain data (JSON-able)."""
    if not fast:
        N = [0.6, 0.8, 1.0, 1.2, 1.4, 1.6, 1.8, 2.0]
        KOP = [1.0, 1.25, 1.5, 1.75, 2.0, 2.5, 3.0]
        cfg = {
            "fast": False,
            "N": N, "k_op": KOP, "d3_kin": [0.0, 0.02, 0.05], "kg": [0.0, 0.3, 0.6],
            "omega": C.mlinspace(0.05, 1.6, 60).tolist(), "drive": C.mlinspace(0.7, 1.3, 9).tolist(),
            "kg_flow": C.mlinspace(0.0, 1.0, 21).tolist(), "tau_kin_n": 4096, "census_nw": 311,
            "loss": {"points": _u([(n, d) for n in N for d in (0.02, 0.05)] + RECORDED_LOSS_PTS),
                     "xi_end": 40.0, "dxi": 0.002, "n": 8192, "L": 400.0, "conv": True},
            "bdg": {"points": [(n, d) for n in N for d in (0.0, 0.02, 0.05)], "n_tau": 500, "tau_max": 16.0},
            "raman": {"points": [(2.0, 0.05), (1.0, 0.05)], "T0_fs": [50.0, 100.0, 200.0, 500.0, 1000.0],
                      "models": ["none", "bw_ss", "bw", "lin3", "ss"], "xi_end": 40.0, "dxi": 0.002, "n": 8192, "L": 400.0,
                      "lambda0_nm": 1550.0, "every": 0.1, "gordon": {"tauR_lin": [0.005, 0.01, 0.02], "T0_bw": [1000.0], "xi_end": 10.0}},
            "scatter": {"points": [(n, k) for n in N for k in KOP], "kg": 0.3, "L": 25.0, "h": 0.01,
                        "conv": [(1.0, 1.5), (1.0, 3.0), (2.0, 1.5), (2.0, 3.0)], "conv_h": 0.005, "conv_L": 35.0,
                        "kaup_omega": [0.8, 1.2, 2.0], "kaup_h": [0.02, 0.01]},
            "ref25": {"variants": ["V1_shipped", "V2_prod", "V3_sg8_shipped", "V3_sg8_prod", "V4_d3draft_prod", "V5_xi12_prod"]},
            "parity": {"odd": [4095, 8191, 16383, 32767, 65535], "even": [4096, 8192, 16384, 32768, 65536],
                       "cont_N": [1.0, 2.0, 3.0]},
        }
    else:
        N = [1.0, 2.0]
        KOP = [1.0, 1.5]
        cfg = {
            "fast": True,
            "N": N, "k_op": KOP, "d3_kin": [0.05], "kg": [0.3],
            "omega": C.mlinspace(0.05, 1.6, 12).tolist(), "drive": C.mlinspace(0.7, 1.3, 5).tolist(),
            "kg_flow": C.mlinspace(0.0, 1.0, 6).tolist(), "tau_kin_n": 4096, "census_nw": 61,
            "loss": {"points": [(1.0, 0.05), (2.0, 0.05)], "xi_end": 40.0, "dxi": 0.002, "n": 2048, "L": 100.0, "conv": False},
            "bdg": {"points": [(1.0, 0.05)], "n_tau": 200, "tau_max": 16.0},
            "raman": {"points": [(2.0, 0.05)], "T0_fs": [50.0], "models": ["none", "bw_ss"], "xi_end": 6.0, "dxi": 0.004,
                      "n": 2048, "L": 100.0, "lambda0_nm": 1550.0, "every": 0.1,
                      "gordon": {"tauR_lin": [0.01], "T0_bw": [], "xi_end": 4.0}},
            "scatter": {"points": [(n, k) for n in N for k in KOP], "kg": 0.3, "L": 20.0, "h": 0.02,
                        "conv": [], "conv_h": 0.01, "conv_L": 25.0, "kaup_omega": [1.2], "kaup_h": [0.02]},
            "ref25": {"variants": ["V1_shipped_fast"]},
            "parity": {"odd": [4095, 8191], "even": [4096, 8192], "cont_N": [1.0, 2.0, 3.0]},
        }
    cfg["op"] = dict(OP)
    cfg["anchor"] = dict(ANCHOR)
    cfg["findings"] = {"N": 2.0, "k_op": 3.0, "omega": FINDINGS_OMEGA, "ratio": FINDINGS_RATIO}
    return cfg


def ref25_variant(name):
    """Settings of one Ref. [25] sweep variant (B8)."""
    shipped = {"Nt": 2 ** 13, "n_steps": 4000, "xi_max": 10.0}
    prod = {"Nt": 2 ** 14, "n_steps": 6000, "xi_max": 10.0}
    d3_shipped = C.mlinspace(0.01, 0.10, 5).tolist()
    std = {"N": [2.5, 3.0, 3.5, 4.0, 4.5], "d3": d3_shipped, "shapes": ["sech", "gaussian", "super_gaussian", "chirped"],
           "masks": [(2.5, 12), (3.0, 10), (4.0, 8)]}
    V = {
        "V1_shipped": dict(std, **shipped, note="run_robustness_sweep.m as shipped (SOTHE_pkg v1.0.0); control against job 530047 S7"),
        "V2_prod": dict(std, **prod, note="production resolution stated by Ref. [25] (N_t = 2^14, 6000 steps)"),
        "V3_sg8_shipped": dict(std, **shipped, shapes=["super_gaussian_m4"],
                               note="super-Gaussian of the version of record, N exp(-tau^8/2), shipped resolution"),
        "V3_sg8_prod": dict(std, **prod, shapes=["super_gaussian_m4"],
                            note="super-Gaussian of the version of record, N exp(-tau^8/2), production resolution"),
        "V4_d3draft_prod": dict(std, **prod, d3=[0.02, 0.04, 0.06, 0.08, 0.10],
                                note="delta3 grid of the Paper-1 draft text (main.md), production resolution"),
        "V5_xi12_prod": dict(std, **dict(prod, xi_max=12.0), note="production resolution with xi_max = 12 (the reference-run length)"),
        "V1_shipped_fast": dict(std, Nt=2 ** 12, n_steps=1500, xi_max=6.0, N=[3.0, 4.0], d3=[0.01, 0.1], shapes=["sech"],
                                masks=[(3.0, 10)], note="smoke test (the 'fast' settings of run_robustness_sweep.m), 4 configurations"),
    }
    if name not in V:
        raise KeyError(name)
    v = dict(V[name])
    v["name"] = name
    return v
