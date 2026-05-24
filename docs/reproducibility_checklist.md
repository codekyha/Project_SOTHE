# Reproducibility checklist

This document maps each figure and each numerical claim in the supplementary package to the exact script that produces it.

## 1. Software dependencies

| Component | Minimum tested version | Notes |
|---|---|---|
| MATLAB | R2021a | No toolboxes required. R2024a+ recommended for `exportgraphics` PDF output. |
| Python | 3.10 | 3.12 used in the shipped validation run. |
| NumPy | 1.22 (with `np.trapz`) or any 2.x | Compat alias handled in code. |
| SciPy | 1.7 | Optional, only used implicitly by Matplotlib internals. |
| Matplotlib | 3.5 | `Agg` backend used by `make_validation_figure.py`; no display required. |

## 2. Headline-number reproduction (~5 minutes total)

```bash
cd codes
python3 validate_paper_claims.py             # full validator, both tests
python3 make_validation_figure.py            # regenerates figures + results
```

The validator prints PASS for both tests. The exact numbers will vary at the third decimal place due to the soft prior in the Cherenkov locator and the inherent breathing noise of $S_{\mathrm{tot}}(\xi)$.

## 3. Per-figure reproduction

| Figure file | Reproducer | Approx. runtime |
|---|---|---|
| `figures/validation_figure.png` | `codes/make_validation_figure.py` | ~20 s |
| `figures/fig_thermo_analysis_sol.png` | MATLAB: production version of `gnlse_dimensionless.m` + a wrapper that draws the 3-panel figure with $S$ trajectories and the Cherenkov scan inset; same physics as `example_*.m` but with the publication aesthetic. | ~3 min |
| `figures/fig_evolution_optical_event_horizon.png` | MATLAB: standalone explanatory diagram with a synthetic two-Gaussian spectrum to illustrate the bipartite partitioning + entropy-density integrand. Static-data figure; not numerically tied to the GSL simulation. | ~1 s |
| `figures/fig_run1.png` ... `fig_run5.png` | These five are the originally-produced outputs of the five robustness runs documented in `results/summary_run1.txt` and `summary_run2.txt`. The production scripts are part of the GitHub repository release and not bundled here to keep the supplementary lightweight. | each: 5–40 min |

## 4. Per-claim reproduction

| Claim | Where | Script |
|---|---|---|
| $\Delta S_{\mathrm{tot}}>0$ at the nominal point | `results/validation_run.txt`, line "Delta S_tot" | `make_validation_figure.py` |
| $\eta_{\mathrm{GSL}}>1/2$ at the nominal point | same | same |
| Photon-number drift $<10^{-3}$ | `results/validation_run.txt` | same |
| $k_{\mathrm{RR}}\delta_3\approx 1$ across 5-point scan | `results/phase_matching_scan.csv`, `results/validation_run.txt` | same |
| Coarse-graining range $\Delta S_{\mathrm{tot}}\in[0.82,2.42]$ across 300 configs | `results/summary_run1.txt` | production sweep (not in package) |
| Boltzmann fit is window-dependent (justifies withdrawal) | `figures/fig_run4.png` | production Run 4 (not in package) |

## 5. What to expect on re-run

Running `validate_paper_claims.py` twice will not give bit-identical output even with fixed seeds, because the FFT-based SSFM accumulates floating-point rounding that depends on machine micro-architecture. The expected variance is below the third decimal place for $c_{\mathrm{RR}}$ and below 5% for $\Delta S_{\mathrm{tot}}$. All headline numbers should remain comfortably inside the paper's reported ranges.

If a re-run *does* miss those ranges, the most likely cause is a numpy/scipy version that has changed FFT internals. As a quick sanity check, the photon-number drift should be at the $10^{-12}$ level or smaller; if it isn't, the integrator is misconfigured (typically a dispersion step that has the wrong half-factor).

End of `reproducibility_checklist.md`.
