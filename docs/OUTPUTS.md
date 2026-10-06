# OUTPUTS: every file of a SOTHE-P2 run folder

A run folder is `runs/<RUNID>/`. `RUNID` is the UTC start time, followed by `_FAST` for `--fast` and `_local` for runs without SLURM. The results tarball `SOTHE-P2_results_<RUNID>.tar.gz` holds the whole folder.

## Top level

| File | Content |
|---|---|
| `RUN.env` | how the run was started: stages, fast, workers or cpus, account, partition, module, pack version |
| `JOBID` | SLURM job id (SLURM runs only) |
| `job.sbatch` | the job script (SLURM runs; `submit --local` too) |
| `job.log`, `slurm-<id>.out`, `slurm-<id>.err` | the console of the whole run; still written after the tarball, so not in the MANIFEST |
| `SUMMARY.md` | per stage: ran, criteria, key numbers, token values, parity with MATLAB, the sweep against Ref. [25], the list of open decisions |
| `summary.json` | the same in JSON |
| `tokens_proposed.json` | for each manuscript placeholder a report supplies: value, source file and field, status. Nothing edits a manuscript |
| `MANIFEST.sha256` | sha256 of every file except the live logs; `sha256sum -c MANIFEST.sha256` passes in place and in an extracted tarball |

## Every stage folder

| File | Content |
|---|---|
| `stage_status.json` | `ok` (the stage ran), `pass` (its criteria held), `checks`, `outputs`, `notes`, `error`, `wall_s`, provenance |
| `stage_exit.json` | exit code (0 = ran, 1 = failed) and wall time |
| `provenance.json` | Python, NumPy, SciPy and Matplotlib versions, BLAS/LAPACK, threads, workers, host, SLURM job id, package version, status |
| `console.log` | the stage's console output |

## S1_g4a: gate G4a (C-05)

| File | Content |
|---|---|
| `g4a_report.json` | the fields of `run_g4a_canonical.m` 0.1.1:<br>T-01: the BdG rate about $\psi_0$;<br>T-02: Newton residuals of the TOD soliton;<br>T-03: annulus and zero sector about $\phi$;<br>T-04: operating-point loss rate and drift, with `T04_drift_grid` as the grid-peak variant;<br>T-05: Eq. (15);<br>T-06: the Gaussian layer and its tokens.<br>Also `G4a_all_pass`, `runtime` |
| `entropy_trajectory.csv` | copy of `data/entropy_trajectory.csv` (the input of T-05) |
| `T04_series.csv` | `xi, Q, tau_c, k_s, tau_c_grid` of the operating-point propagation, every 0.5 in $\xi$ |
| `bdg_backgrounds.npz` | $\phi$ at $\delta_3=0.02,\,0.05$ with the Newton histories; BdG spectra, Krein signs and zero sectors about $\phi$ and about $\psi_0$ |

## S2_g4c: G4c (C-05b, T-10, T-12)

| File | Content |
|---|---|
| `g4c_loss_report.json` | fields of `run_g4c_loss_scan.m` 0.1.2 plus the peak-estimator fields.<br>`points`: 9 points $(N,\delta_3)$ with window rates and drifts (`window_drift`, sub-grid; `window_drift_grid`), $c$, losses, recoil deviation (both peaks), Newton tail, convergence runs.<br>`scaling`, `fit_diagnostic`.<br>`channels`: lab frame, five epochs, $(N, k_{\rm op})$ grid.<br>`tokens` (sub-grid peak), `tokens_grid_peak`, `tokens_changed_by_peak`, `protocol`, `runtime` |
| `series/{base,half,n2}_N<N>_d<delta3>.csv` | all 27 propagations (base, $\Delta\xi/2$, $2n$; fast run: base only); columns as `T04_series.csv` |
| `newton_solutions.npz` | the 9 periodic Newton solutions with their histories |

## S3_g12: provenance of the archived trajectory (G-12, T-19)

| File | Content |
|---|---|
| `g12_report.json` | runs A to D: settings, $\Delta S_{\rm tot}$, Eq. (15) values, photon-number drift, maximum deviation from the archive, verdict.<br>`src_port_check`: run A through `ref25/src.py`.<br>`mask_rule`, `decision_needed` |
| `traj_<tag>_N<N>_d<delta3>_Nt<Nt>_ns<steps>.csv` | each trajectory in the archive's format: `xi, S_hor, S_rad, S_tot, eta_GSL, P_norm` |

## S4_oracles: oracle ledgers (T-14)

| File | Content |
|---|---|
| `oracle_ledger_v010.json` | 29 rows: id, name, pass, value, reference, note |
| `oracle_ledger_v020.json` | 35 rows with class (ACCURACY, STRUCTURAL, REGRESSION) and change (kept, modified, added); retired ids and their pins; diagnostics: Eq. (15), N-safe law, envelope, TOD solitons, loss, Kaup |
| `physics_pass.npz` | the arrays of the physics pass (small stability grid), flat keys as in `phase2_bundle.npz` |

## S5_figures: the manuscript figures (T-13, T-18)

| File | Content |
|---|---|
| `fig1_thermality.pdf`, `.eps`, `.png` | Fig. 1 (15 cm wide) |
| `fig3_gapless_flow.pdf`, `.eps`, `.png` | Fig. 2 (8.5 cm) |
| `fig2_bdg_krein.pdf`, `.eps`, `.png` | Fig. 3 (15 cm), with the zero-sector inset about $\phi$ and $\psi_0$ |
| `fig5_radiative_loss.pdf`, `.eps`, `.png` | Fig. 4 (8.5 cm; needs `S2_g4c/g4c_loss_report.json` or `P2_G4C_REPORT`) |
| `fig6_partner_entanglement.pdf`, `.eps`, `.png` | Fig. 5 (8.5 cm) |
| `captions.md`, `captions.tex` | the caption of each figure (explanation only; the analysis belongs to the paper text), in Markdown and as `\caption{...}` lines |
| `figs_manifest.json` | version, provenance, status, the captions, and for each figure its number and the numbers it shows |
| `figure_data.npz` | every plotted array |

The figures are written at their final size in the publication style of `sothe_p2/pubstyle.py` (IOP guidelines; see the README). PDF and EPS are vector and carry no creation date; PNG files are 600 dpi.

## S6_phase2: phase-2 data products (suite v0.1.0)

| File | Content |
|---|---|
| `paper2_data/thermality_ratio.csv` | `omega, ln_ratio, ratio, alpha2, beta2` (60 rows) |
| `paper2_data/kappa_drive_independence.csv` | `drive, kappa_kin` (9) |
| `paper2_data/kinematic_kappa_flow.csv` | `kappa_g, kappa_kin, T_H` (21) |
| `paper2_data/greybody.csv` | `omega, Gamma` (60) |
| `paper2_data/partner_log_negativity.csv` | `omega, E_N, E_N_cross, nu_minus` (60) |
| `paper2_data/radiation_rate_table.csv` | `delta3, maxImOmega` (6) |
| `paper2_data/bdg_spectrum_d3_000.csv`, `bdg_spectrum_d3_020.csv` | `Re_omega, Im_omega, krein` in the window $\lvert{\rm Re}\,\omega\rvert\le3\mu$ |
| `paper2_data/false_positive_map.csv` | `kappa_g, drive, beta_abs` (61 x 61, column-major) |
| `paper2_data/stability_map.csv` | `N, delta3, maxImOmega` (13 x 13, column-major) |
| `paper2_data/phase2_summary.json` | the fields of `run_phase2_all.m`, then `port` and `status` |
| `paper2_figures/fig1_thermality` ... `fig6_partner_entanglement` (`.pdf`, `.eps`, `.png`) | the six figures of suite v0.1.0, in the publication style of S5 (fig1 and fig2 15 cm wide, the others 8.5 cm) |
| `paper2_figures/captions.md`, `captions.tex` | their captions |
| `phase2_bundle.npz`, `phase2_bundle.mat` | the physics pass. The `.mat` holds the structs `R` and `P` as MATLAB saves them and is written only when SciPy is installed |
| `oracle_ledger_v010_fullgrid.json` | the 29-oracle ledger on the full-grid pass |
| `parity_vs_matlab.json` | per file: header, rows, byte identity, and per column the maximum absolute and relative deviation. For the spectra: eigenvalue-set distance and Krein counts. For the summary JSON: field-by-field differences |

The $\delta_3>0$ BdG products (`radiation_rate_table.csv`, `bdg_spectrum_d3_020.csv`, `stability_map.csv`, fig2 (b), fig5) linearize about $\psi_0$ (G3u audit, S1-1). They are regression data of suite v0.1.0, not results.

## S7_sweep: Ref. [25] robustness sweep and checks

| File | Content |
|---|---|
| `robustness_sweep.json` | grids, settings, `DeltaS_range`, `eta_range`, `pass_count`, `pass_rate`, the claims of Ref. [25], breakdown by pulse shape and by mask scheme, and one row per configuration |
| `robustness_sweep.npz` | `DeltaS`, `eta`, `photon_drift_max` as $5\times5\times4\times3$ arrays (N, $\delta_3$, shape, mask), and the grids |
| `robustness_table.csv` | one row per configuration: `N_sol, delta3, pulse_shape, mask_width, mask_order, DeltaS_tot, eta_GSL_final, pass, photon_drift_max` |
| `ref25_checks.json` | the nominal run ($N=3.5$, $\delta_3=0.02$); the phase-matching scan ($k_{\rm RR}$ at $\delta_3=0.06\ldots0.10$, $c_{\rm RR}$, $R^2$, the invariant $k_{\rm RR}\delta_3$); one `gnlse_dimensional` run |

## Probe folders (`runs/probe_<UTC>/`)

| File | Content |
|---|---|
| `probe.txt` | the probe output, ending with `PROBE OK` or `PROBE FAILED` |
| `S0_probe/` | `probe.json` (smoke-test numbers), `stage_status.json`, `provenance.json` |
| `probe.sbatch`, `JOBID`, `slurm-<id>.*` | for SLURM probes |
