# SOTHE — Supplementary Material

Companion data and code for the revision of:

> H. Oguz, *Generalized Thermodynamics of Solitonic Event Horizons in Dispersive Field Theories*, **Classical and Quantum Gravity** (manuscript CQG-114816, revision submitted 2026).

This is a **lightweight** supplementary package: it contains exactly what is needed to reproduce the paper's three headline claims (coarse-grained generalized second law, Hamiltonian photon-number conservation, Cherenkov phase-matching $k_{\mathrm{RR}}\delta_3\approx 1$), the figures that appear in the paper, and a small set of pre-computed numerical results. A larger codebase will accompany the published version on GitHub.

A summary of what the package corrects relative to the original 2025 submission is in `NOTES.md`. A walkthrough of the physics with all formulae spelled out is in `EXPLANATION.md`.

---

## 1. Directory layout

```
SOTHE_Supplementary/
├── README.md                         this file (entry point)
│
├── codes/                            MATLAB sources (primary) + Python validator
│   ├── spectral_entropy.m            Shannon entropy on a masked spectrum (Eq. 2/5)
│   ├── soliton_mask.m                adaptive super-Gaussian soliton-following mask
│   ├── eta_gsl.m                     coarse-grained GSL efficiency (Eq. 8)
│   ├── gnlse_dimensionless.m         symmetric SSFM solver, corrected TOD sign (Eq. 9)
│   ├── find_rr_lobe.m                robust Cherenkov-peak locator
│   ├── rr_phase_match_fit.m          through-origin fit k_RR = c_RR / delta3 (Eq. 15)
│   ├── example_single_run.m          one-click nominal-point reproducer
│   ├── example_phase_matching_scan.m one-click 5-point delta_3 scan
│   ├── validate_paper_claims.py      pure-Python independent validator (~5 min)
│   └── make_validation_figure.py     reproduces figures/validation_figure.png
│
├── figures/                          publication figures
│   ├── fig_evolution_optical_event_horizon.png  Fig. 1 (bipartite partition + entropy density)
│   ├── fig_thermo_analysis_sol.png              Fig. 2 (spectral evolution + GSL + Cherenkov)
│   ├── fig_run1.png                             Run 1: reformulated-entropy GSL
│   ├── fig_run2.png                             Run 2: truncated vs full GNLSE robustness
│   ├── fig_run3.png                             Run 3: (N, delta3) parameter sweep
│   ├── fig_run4.png                             Run 4: thermal-fit robustness (justifies withdrawal)
│   ├── fig_run5.png                             Run 5: pulse-shape sensitivity
│   └── validation_figure.png                    Python validator output, fresh from this package
│
├── results/                          machine-readable numerical results
│   ├── validation_run.json           full structured output of make_validation_figure.py
│   ├── validation_run.txt            human-readable summary
│   ├── entropy_trajectory.csv        xi, S_hor, S_rad, S_tot, eta_GSL, P_norm
│   ├── phase_matching_scan.csv       delta3, k_RR, k_RR*delta3
│   ├── summary_run1.txt              entropy GSL across mask schemes (truncated GNLSE)
│   └── summary_run2.txt              truncated vs full GNLSE comparison (Raman/beta4)
│
└── docs/
    └── reproducibility_checklist.md  what to run to regenerate every figure and number
```

---

## 2. Quick start (60 seconds)

If you have **Python** with NumPy and Matplotlib:

```bash
cd codes
python3 validate_paper_claims.py --fast        # ~1 min, sanity check
python3 make_validation_figure.py              # ~20 s, regenerates figures/validation_figure.png
```

Expected output (your numbers will differ at the 1% level due to grid noise):

```
[1] Nominal-point GSL & soliton breathing
  Delta S_tot         =   +1.81 nats          (paper: [0.82, 2.42])
  eta_GSL(xi_f)       =    0.55                (paper: [0.52, 0.62])
  fraction dS/dxi<0   =     48 %               (paper qualitative: ~50%)
  photon-number drift = +1.4e-12               (target |drift| < 1e-3)
  -> PASS

[2] Cherenkov phase-matching scaling
  c_RR             = 0.99    (paper: 1.01)
  mean k_RR*delta3 = 1.02    (paper: 1.02)
  -> PASS
```

If you have **MATLAB** (R2021a or later, no toolboxes required):

```matlab
addpath('codes');
example_single_run;                            % ~30 s, nominal-point GSL
example_phase_matching_scan;                   % ~5 min, full delta_3 scan
```

---

## 3. What each script does in one sentence

| Script | Role |
|---|---|
| `spectral_entropy.m` | Shannon entropy of a normalized spectral distribution on a mask. Same routine for $S_{\mathrm{hor}}$ and $S_{\mathrm{rad}}$. |
| `soliton_mask.m` | Builds a super-Gaussian mask centred on the moving soliton centroid. |
| `eta_gsl.m` | Coarse-grained efficiency $\eta_{\mathrm{GSL}} = \sum_+|\Delta S|/\sum|\Delta S|$. $1/2$ is the random-walk baseline. |
| `gnlse_dimensionless.m` | Symmetric split-step solver with the integrability-consistent TOD sign $D(k)=i(-k^2/2+\delta_3 k^3)$. |
| `find_rr_lobe.m` | Locates the Cherenkov peak in the final spectrum, with sub-grid refinement. |
| `rr_phase_match_fit.m` | Through-origin LS fit of $k_{\mathrm{RR}}=c_{\mathrm{RR}}/\delta_3$. |
| `example_single_run.m` | Single nominal-point reproducer with a one-line headline-numbers report. |
| `example_phase_matching_scan.m` | 5-point $\delta_3$ scan reproducing the paper's phase-matching law. |
| `validate_paper_claims.py` | License-free independent Python implementation; CI-friendly. |
| `make_validation_figure.py` | Wraps the above and writes `results/` + `figures/validation_figure.png`. |

---

## 4. Headline numbers (from `results/validation_run.txt`)

| Quantity | Paper value | Python validator (this package) |
|---|---|---|
| $\Delta S_{\mathrm{tot}}$ (nats) | $[0.82, 2.42]$ across 300 configs | $+1.81$ at the nominal point |
| $\eta_{\mathrm{GSL}}(\xi_f)$ | $[0.52, 0.62]$ across 300 configs | $0.55$ at the nominal point |
| Photon-number drift | $< 10^{-3}$ | $1.4\times 10^{-12}$ |
| $c_{\mathrm{RR}}$ | $1.01$ | $0.99$ |
| $\langle k_{\mathrm{RR}}\delta_3\rangle$ | $1.02$ | $1.02$ |

The Python harness uses lighter grids ($N_t=2^{13}$, $n_{\text{steps}}=4000$) than the production MATLAB runs ($N_t=2^{14}$, $n_{\text{steps}}=6000$), which accounts for the residual sub-percent differences.

---

## 5. License

The MATLAB and Python sources are released under the **MIT License**. The figures, datasets and prose are © the author, all rights reserved under the journal's standard publication agreement.

---

## 6. Contact

H. Oguz · Department of Computer Technologies, Vocational School, Istanbul Okan University · hasan.oguz@okan.edu.tr

For supplementary-package bug reports only (not paper content): please include the output of `validate_paper_claims.py`, your NumPy/MATLAB versions, and a minimal reproducer.
