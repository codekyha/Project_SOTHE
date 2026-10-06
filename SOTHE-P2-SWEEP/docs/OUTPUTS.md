# SOTHE-P2-SWEEP: output files

A run folder `runs/<RUNID>/` (and the results tarball `SOTHE-P2-SWEEP_results_<RUNID>.tar.gz`) contains the files below. CSV files have one header row, comma separators and floats at `%.17g` (`nan` where a quantity is undefined, `1`/`0` for booleans). Every number is RECORDED (HR-3).

## Top level

| File | Content |
|---|---|
| `SUMMARY.md` | run facts, the oracle ledger, block-by-block summary tables (start here) |
| `ORACLES.json` | every oracle: id, block, class (`hard`, `claim`, `info`), statement, value, reference, tolerance, pass |
| `summary.json` | all reduced numbers (`results.<block>`), task timing per block, failed tasks, phase-C errors, exit code |
| `config.json` | the grids and protocols of the run (`cfg`) and the expanded block list |
| `provenance.json` | package version, host, SLURM job id, workers, Python/NumPy/SciPy/Matplotlib versions, BLAS, thread settings (one entry per session: submit, resume) |
| `plan_A.json`, `plan_B.json` | the task lists of phases A and B (id, block, function, arguments, cost estimate) |
| `RUN.env`, `JOBID`, `EXIT_CODE` | run settings, SLURM job id(s), exit code (0 ok, 2 hard oracle failed, 3 task failed) |
| `job.sbatch`, `resume.sbatch` | the job scripts |
| `job.log`, `slurm-<id>.out`, `slurm-<id>.err` | logs (not in the MANIFEST) |
| `MANIFEST.sha256` | SHA-256 of every other file |
| `tasks/<block>/<id>.json` | the raw result of every task (`result` field), with its arguments and wall time; `<id>.error.json` holds the traceback of a failed task |

## B1_kinematics/

- `kin_grid.csv`: `N, k_op, d3, kg, S, c_match, lambda, kappa_closed, kappa_num, deficit_rel, eq7pp_rel, eq7pp_dev, drive_spread_abs, drive_spread_rel, band_rel, band_abs, eq7pp_drive_max_dev, mirror_max, T_H, HP, mu` (one row per grid point; `kappa_num` at drive 1; `band_*` = $\lambda^2\Delta\tau^2/4$).
- `kin_flows.csv`: `N, k_op, d3, kappa_g, kappa_num, kappa_closed, T_H` (the $\kappa_g$ flows).

## B2_entanglement/

- `entanglement_grid.csv`: `N, k_op, d3, kg, kappa, T_H, EN_max_band, nu_max_band, IR_asymptote, EN_at_mu, lab_lo, lab_hi, lab_frac, EN_lab_lo, EN_lab_hi, late_lo, late_hi, late_frac, EN_late_lo, EN_late_hi` (window edges in the laboratory frame and in the frame of the last drift epoch $[30,40]$).
- `EN_curves_op.csv`, `EN_curves_anchor.csv`: `omega, E_N, nu_minus, nbar_max, E_N_eta0.9, E_N_eta0.5, E_N_eta0.1`.

## B3_loss/

- `loss_points.csv`: `N, d3, Nd3, loss_rate, loss_0_40, loss_first_hp, c_over_xi, kappa_hp, hawking_period, k_s_early, k_s_late, recoil_rel_dev, drift_leading_order, rel_dev, rel_dev_0_40, drift_conv_max_abs, converged`, then `window_drift_1..5`, `window_rate_1..5`, `window_ks_rate_1..5` (windows $[0,5],[5,10],[10,20],[20,30],[30,40]$).
- `series/loss_N<N>_d<d3>.csv`: `xi, Q, tau_c, k_s, tau_c_grid` (base run, every 0.5).
- `loss_report.json`: per-point reductions and the comparison with job 530047.

## B4_census/

- `census.csv`: `N, d3, k_op, kg, frame, Vd, S, c_match, two_c, omega_max_pair_used, omega_max_pair, omega_max_pair_no_tod, Phi_pair, Phi_up_low, pair_open_below, k_star, mu, kappa, window_lo, window_hi, window_frac, piN_over_S, n_at_mu, EN_at_mu, epoch_lo, epoch_hi` (`frame` = `lab` or `xi<lo>-<hi>`; `window_*` = the Hawking window, empty when closed).
- `n_at_mu`, `EN_at_mu` in the census files use the closed-form $\kappa_0=NS(k_{\rm op})$ (the formula of REVIEW_G5 S1-1); `entanglement_grid.csv` evaluates $E_N$ at the extracted $\kappa$ (they differ by the Eq. (7″) deficit, a few $10^{-4}$ relative).
- `hawking_window_d3_005.csv`: the S1-1 table at $\delta_3=0.05$, $\kappa_g=0.3$ for every $(N,k_{\rm op})$: `N, k_op, mu, kappa, piN_over_S, n_at_mu, EN_at_mu, lab_edge, lab_lo, lab_hi, lab_frac, late_Vd, late_edge, late_lo, late_hi, late_frac, late_Phi_pair`.

## B5_bdg/

- `bdg_points.csv`: `N, d3, n_tau, tau_max, mu, newton_converged, newton_res, newton_iters, kmean, tail_rel, annulus_maxIm, gap_over_mu, n_annulus, krein_pos, krein_neg, zero_sector_max, rho, sqrt_eps_rho, dim_ker_M, dim_ker_M2, sv_M_gap, sv_M2_gap, psi0_annulus_maxIm, psi0_annulus_over_sqrt_d3, psi0_quartet_maxIm, psi0_quartet_over_sqrt_d3`.
- `spectra/bdg_N<N>_d<d3>.csv`: `Re_omega, Im_omega, krein` for $|{\rm Re}\,\omega|\le3\mu$ (the format of `bdg_spectrum_d3_*.csv`).

## B6_raman/

- `raman_runs.csv`: `N, d3, T0_fs, model, f_R, s, tauR_eff, gordon_rate, initial_ks_rate, initial_kfull_rate, early_ks_rate, early_kfull_rate, shift_total_ks, loss_end, k_s_end, k_full_end, tau_c_end, n_roll`, then `window_drift_i`, `window_ks_rate_i`, `window_kfull_rate_i`.
- `raman_S1-3_table.csv`: `N, d3, model, T0_fs, tauR_eff, gordon_rate, initial_rate_ks, early_rate_ks, dk_first_HP, HP, k_op_for_HP, shift_total_ks, linear_extrapolation, recoil_dVd_dxi, times_recoil, late_drift, loss_end`.
- `raman_gordon.csv`: `label, tauR_eff, rate_full, rate_windowed, gordon, ratio, photon_drift`.
- `raman_census.csv`: `N, d3, model, T0_fs, k_op, xi_lo, xi_hi, Vd, omega_max_pair, Phi_pair, Phi_up_low, window_lo, window_hi, window_frac` (the comoving census along each run).
- `series/raman_N<N>_d<d3>_<model>_T<T0>.csv`: `xi, Q, tau_c, k_s, k_full, P` (every 0.1).

Models: `none` (Kerr + TOD), `ss` (+ self-steepening), `lin3` (+ linear Raman, $T_R=3$ fs), `bw` (+ Blow–Wood Raman), `bw_ss` (Blow–Wood + self-steepening).

## B7_scatter/

- `scatter_omega.csv`: `N, k_op, omega, ins, outs, two_sided, r_conf, n_H, planck, boltzmann, n_H_over_planck, r_conf_over_boltzmann, fluxbal, resid` (channel labels: side L/R + component u/v).
- `scatter_points.csv`: per $(N,k_{\rm op})$ summary, `... frac_two_sided, fluxbal_max, r_conf_min, r_conf_max, r_conf_spread_rel, n_H_over_planck_min, n_H_over_planck_max, kappa_fit_n_H`.
- `scatter_convergence.csv`: $|S|^2$ changes under $h\to h/2$ and $L\to35$.
- `findings_check.csv`: the FINDINGS table at $(2,3)$ against the solve.
- `kaup_null.json`: the Kaup and no-soliton checks.

The full $|S|^2$ matrices by label are in `tasks/B7/*.json` (`omega_records[*].S2`).

## B8_ref25/

- `ref25_<variant>.csv`: `N_sol, delta3, pulse_shape, mask_width, mask_order, DeltaS_tot, eta_GSL_final, pass, photon_drift_max, wall_s` in MATLAB ndgrid order (the columns of job 530047's `robustness_table.csv`).
- `ref25_variants.json`: per variant the ranges, pass count, counts outside the Paper-1 ranges, subset ranges (per pulse, per mask), extremes; for V1 the comparison with job 530047.
- `ext_identity.json`: the bitwise identity of the extended solver.

## B9_parity/

- `table2_parity.csv`: `parity, n_tau, deficit_over_lamdt2, eq7pp_over_lamdt2, s_over_dt, published, dev_published`.
- `table3_continuum.csv`: `N, extracted, eq7pp, relative, published_extracted, published_eq7pp`.

## B10_anchor/

- `paper2_data_anchor/thermality_ratio.csv` (`omega, ln_ratio, ratio, alpha2, beta2`), `kappa_drive_independence.csv` (`drive, kappa_kin`), `kinematic_kappa_flow.csv` (`kappa_g, kappa_kin, T_H`), `partner_log_negativity.csv` (`omega, E_N, E_N_cross, nu_minus`): the paper2 formats of `run_phase2_all`, at $(N,k_{\rm op},\delta_3,\kappa_g)=(1,1.5,0.05,0.3)$.
- `anchor_summary.json`: $\kappa$, $T_H$, $\mu$, Hawking period, thermality fit, flux residual, drive spread, $E_N$ at $0.05$, $\mu$ and $1.6$, the census per drift epoch, the loss at $(1,0.05)$ (also over one anchor Hawking period after $\xi=10$), the BdG block at $(1,0.05)$, the $\delta_3=0$ scattering summary at $(1,1.5)$.

## figures/

`sweep_fig1_hawking_window` … `sweep_fig6_loss_drift` as PDF and EPS (vector) and PNG (600 dpi), with `captions.md` and `captions.tex`. Diagnostic figures of recorded numbers in the publication style of SOTHE-P2 (8.5 or 15 cm wide, Computer Modern, no titles).
