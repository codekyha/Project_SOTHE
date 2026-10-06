# SOTHE-P2-RAD: the files of a results folder

`python rad.py run --out results` writes `results/<UTC stamp>/`:

| File | Content |
|---|---|
| `catalogue.json` | the runs of this folder with their parameters (`sothe_rad/experiments.py`) |
| `RUN.json` | package and version, start and end (UTC), command, workers, platform and library versions, wall times |
| `runs/<id>/meta.json` | the full configuration of the run (defaults filled in), wall time, steps, grid spacing |
| `runs/<id>/series.npz` | the diagnostics every `sample` units of $\xi$ (below) |

Arrays of `series.npz` (one value per sample unless stated): `xi`, `d3` (third-order coefficient at that $\xi$), `Q` (core
norm), `Qbox` (norm in the region without absorption), `tc` (peak position in the computational frame; add $V\xi$ for the
laboratory frame), `tgrid` (grid maximum), `amp` (peak amplitude), `ks_hann` (Hann-windowed spectral centroid of the core,
the SOTHE-P2 definition), `ks_mom` (momentum centroid of the core), `Mbox` (momentum in the region without absorption),
`phase` (phase of the field at the peak), `E` (Hamiltonian of the whole grid), `Iband`, `Iband_near`, `Iband_far`
(band-filtered intensity in the tail window: whole flat part, nearer and farther half), `kband` (power-weighted wavenumber in
the band), `kc` (predicted resonant wavenumber that centres the band), `Iband_prof` (eight bins, far to near) and
`Iband_dist` (their distances behind the peak).

`python rad.py analyse` adds `analysis/` and `SUMMARY.md`:

| File | Content |
|---|---|
| `steady_rates.csv` | prepared solitons over the hold phase: $\tilde\delta$ (mean, minimum, maximum), the rate from the flux (mean and relative spread), from the core and from the box, the flag `steady` (spread below 5 percent), band and measured resonant wavenumbers, model wavenumbers (drift-corrected and without drift), velocity, phase rate, first-order rate and the ratio to it |
| `rate_function.json` | $C$, $p$, $g$ of $R_s=C\tilde\delta^{\,p}e^{g\tilde\delta-\pi q_r}$ (calibration on the prepared solitons), residuals overall and per run, the fits to the steady runs alone (two and three parameters, drift-corrected and no-drift root), the ratio-to-first-order law, the slope of $\ln R$ against $1/(N\delta_3)$ |
| `universal_points.csv` | every rate used: its run, `role` (`calibration`: prepared, hold phase; `test`: launched, emitted after $\xi_e=20$; `cross-check`: core slope of the prepared runs), estimator, state and ratio to the rate function |
| `recoil_model_vs_gnlse.csv` | launched runs at $\xi=20,40,80,120,160$ (as far as each run goes): norm lost since $\xi=10$ and total loss, $k_s$ and velocity, of the run and of the model; for the runs to $\xi=160$ the two calibrated variants |
| `momentum_balance.csv` | runs to $\xi=160$: median ratio of $dk_s/d\xi$ to $-(k_b-k_s)R$ with its 10 and 90 percent quantiles |
| `window_exponents.csv` | $c$ over $[10,20]$, $[20,40]$, $[40,80]$, $[80,160]$ for run and model, and the model's $1/K$ at the window midpoint |
| `loss_table_reproduction.csv` | the launched points of the loss table of the article: loss over $[0,40]$, drift over $[30,40]$, printed values |
| `numerical_checks.csv` | checks against their base runs |
| `model_tracks.npz` | model trajectories for the figure |
| `summary.json` | everything above in one file, with the comparison of the launched runs with the rate function per run and per group |

`python rad.py figures` adds `figures/fig_emission_recoil.{pdf,eps,png}` and `figures/figure_data.npz` (every plotted array):
(a) $R/N^2$ against $1/\tilde\delta$ (steady runs, calibration and test rates, rate function, first-order rate); (b) loss, (c) $k_s$
and (d) the window exponents $c$ of three launched runs ($\delta_3=0.08$, $0.10$, $0.12$) with the recoil model.
