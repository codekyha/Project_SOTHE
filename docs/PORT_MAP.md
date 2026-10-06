# PORT_MAP: P2_altay_pack 1.0.0 (MATLAB) to SOTHE-P2 1.0.0 (Python)

The MATLAB pack has 68 `.m` files in five folders. Every one of them has a Python counterpart, listed below with the file's version and the stage that runs it. The LaTeX tools of the MATLAB pack (`tex/`) are not part of SOTHE-P2.

Status of every number SOTHE-P2 writes: RECORDED. Under HR-3, only MATLAB R2025b output can be canonical.

## `matlab/kit`: the counssug g4a_kit (9 files)

| MATLAB file | Version | Python | Stage |
|---|---|---|---|
| `p2_bdg_spectrum_bg.m` | 0.1.0 | `kit.bdg_spectrum_bg` | S1, S4, S5 |
| `p2_channel_census.m` | 0.1.0 | `kit.channel_census` | S2 |
| `p2_fmt_sci.m` | 0.1.0 | `kit.fmt_sci` | S2 |
| `p2_gauss_layer.m` | 0.1.0 | `kit.gauss_layer` | S1, S5 |
| `p2_gsl_eta15.m` | 0.1.0 | `kit.gsl_eta15` | S1, S3, S4 |
| `p2_splitstep_loss.m` | 0.1.1 | `kit.splitstep_loss`, `kit.subgrid_peak` | S1, S2, S4, S5 |
| `p2_tod_soliton.m` | 0.1.0 | `kit.tod_soliton` | S1, S2, S4, S5 |
| `run_g4a_canonical.m` | 0.1.1 | `g4a.run` | S1 |
| `run_g4c_loss_scan.m` | 0.1.2 | `g4c.run` | S2 |

## `matlab/pack`: stage drivers, figures, provenance (13 files)

| MATLAB file | Version | Python | Stage |
|---|---|---|---|
| `make_paper2_figs.m` | 0.2.1 | `figures.make_paper2_figs` (publication style, no watermark) | S5 |
| `p2pack_json_write.m` | 1.0.0 | `jsonio.write` | all |
| `p2pack_provenance.m` | 1.0.0 | `runtime.provenance` | all |
| `p2pack_run_stage.m` | 1.0.0 | `stages.run_stage`, `cli` command `_job` | all |
| `p2pack_setup.m` | 1.0.0 | `p2.py` (path) and `pubstyle.setup` (headless Agg, publication settings) | all |
| `p2pack_status.m` | 1.0.0 | `runtime.Status` | all |
| `p2pack_status_write.m` | 1.0.0 | `runtime.Status.finish` | all |
| `stage_S0_probe.m` | 1.0.0 | `stages.stage_S0` | S0 |
| `stage_S1_g4a.m` | 1.0.0 | `stages.stage_S1` | S1 |
| `stage_S2_g4c.m` | 1.0.0 | `stages.stage_S2` | S2 |
| `stage_S3_g12.m` | 1.0.0 | `stages.stage_S3`, `g12.run` | S3 |
| `stage_S4_oracles.m` | 1.0.0 | `stages.stage_S4` | S4 |
| `stage_S5_figures.m` | 1.0.0 | `stages.stage_S5` | S5 |

## `matlab/ref25`: SOTHE_pkg v1.0.0 `src/` (Ref. [25]) (8 files)

| MATLAB file | Python | Stage |
|---|---|---|
| `eta_gsl.m` | `ref25.src.eta_gsl` | S7 |
| `find_rr_lobe.m` | `ref25.src.find_rr_lobe` | S7 (phase-matching check) |
| `gnlse_dimensional.m` | `ref25.src.gnlse_dimensional` | S7 (physical-units check) |
| `gnlse_dimensionless.m` | `ref25.src.gnlse_dimensionless` | S3 (cross-check), S7 |
| `rr_phase_match_fit.m` | `ref25.src.rr_phase_match_fit` | S7 (phase-matching check) |
| `run_robustness_sweep.m` | `ref25.src.run_robustness_sweep` | S7 |
| `soliton_mask.m` | `ref25.src.soliton_mask` | S3 (cross-check), S7 |
| `spectral_entropy.m` | `ref25.src.spectral_entropy` | S3 (cross-check), S7 |

The Python validation port of the same package (`validation/gnlse_dimensionless.py`, unchanged in `sothe_p2/ref25/`) stays the solver of stage S3. S3 has to reproduce the archived trajectory, and that port uses the arithmetic that wrote it (FINDING G-12).

## `matlab/suite`: sothe_phase2_matlab v0.1.0 (35 files)

| MATLAB file | Python | Stage |
|---|---|---|
| `fig1_thermality.m` | `suite_figs.fig1_thermality` | S6 |
| `fig2_bdg_krein.m` | `suite_figs.fig2_bdg_krein` | S6 |
| `fig3_gapless_flow.m` | `suite_figs.fig3_gapless_flow` | S6 |
| `fig4_false_positive.m` | `suite_figs.fig4_false_positive` | S6 |
| `fig5_stability_map.m` | `suite_figs.fig5_stability_map` | S6 |
| `fig6_partner_entanglement.m` | `suite_figs.fig6_partner_entanglement` | S6 |
| `p2_bdg_spectrum.m` | `suite.bdg_spectrum` | S0, S1, S4, S6 |
| `p2_beta0.m` | `suite.beta0` | S4, S6 |
| `p2_bogoliubov_thermal.m` | `suite.bogoliubov_thermal` | S4, S6 |
| `p2_compute_all.m` | `suite.compute_all` | S4, S5, S6 |
| `p2_entropy_reductions.m` | `suite.entropy_reductions` | S4, S6 |
| `p2_fd_circulant.m` | `suite.fd_circulant` | S1, S4, S6 |
| `p2_flow_profile.m` | `suite.flow_profile` | S4, S6 |
| `p2_flux_residual.m` | `suite.flux_residual` | S4, S6 |
| `p2_grad_uniform.m` | `suite.grad_uniform` | S4, S6 |
| `p2_greybody.m` | `suite.greybody` | S4, S6 |
| `p2_horizon_strength.m` | `suite.horizon_strength` (diagnostic, unused) | none |
| `p2_kinematic_kappa.m` | `suite.kinematic_kappa` | S4, S6 |
| `p2_label_contours.m` | `suite_figs.label_contours` | S6 |
| `p2_magma.m` | `suite_figs.magma` | S6 |
| `p2_no_horizon_null.m` | `suite.no_horizon_null` | S4, S6 |
| `p2_partner_logneg.m` | `suite.partner_logneg` | S4, S6 |
| `p2_phase2_params.m` | `suite.phase2_params` | S4, S5, S6 |
| `p2_read_csv_numeric.m` | `suite.read_csv_numeric` | S1, S3, S4, S6 |
| `p2_slow_light.m` | `suite.slow_light` | S4, S5, S6 |
| `p2_soliton_mu.m` | `suite.soliton_mu` | (library) |
| `p2_surface_gravity.m` | `suite.surface_gravity` | S4, S6 |
| `p2_thermality_fit.m` | `suite.thermality_fit` | S4, S6 |
| `p2_vg0.m` | `suite.vg0` | S4, S6 |
| `p2_vg_dressed.m` | `suite.vg_dressed` | (library) |
| `p2_viridis.m` | `suite_figs.viridis` | S6 |
| `p2_write_csv.m` | `phase2.write_csv` | S6 |
| `p2_write_json.m` | `phase2.write_json`, `phase2.json_text` | S6 |
| `run_oracles.m` | `oracles.run_oracles` | S4, S6 |
| `run_phase2_all.m` | `phase2.run_phase2_all` | S6 |

## `matlab/suite_v020`: the T-14 overlay (3 files)

| MATLAB file | Version | Python | Stage |
|---|---|---|---|
| `p2_kaup_residual.m` | 0.2.0 | `kit.kaup_residual` | S0, S4 |
| `p2_phase2_params_v020.m` | 0.2.0 | `suite.phase2_params_v020` | S4 |
| `run_oracles_v020.m` | 0.2.0 | `oracles.run_oracles_v020` | S4 |

Count: 9 + 13 + 8 + 35 + 3 = 68.

## The rest of the MATLAB pack

| MATLAB pack | SOTHE-P2 |
|---|---|
| `bin/altay.sh`, `bin/lib.sh` (install, probe, submit, status, show, repack, verify) | `p2.py` = `sothe_p2/cli.py`, same commands plus `run`, `list`, `stages`, `version` |
| `bin/job_main.sbatch`, `bin/job_probe.sbatch` | written for each run by `cli.write_job_script` (`runs/<RUN>/job.sbatch`, `probe.sbatch`) |
| `python/postprocess.py` | `sothe_p2/postprocess.py` (with S6 and S7 sections) |
| `python/g12_numpy_check.py` | `g12.run` (S3 runs the Python solver itself) |
| `python/ref25_port/gnlse_dimensionless.py` | `sothe_p2/ref25/gnlse_dimensionless.py` (unchanged) |
| `recorded/*.json` (Octave 8.4.0 reports, phase-0 validation run) | `reference/recorded_octave840/` |
| `config/site.env.example`, `data/entropy_trajectory.csv` | same paths (the trajectory is byte-identical) |
| `tex/` (LaTeX export and checks) | not included |

## MATLAB semantics reproduced (`sothe_p2/compat.py`)

| MATLAB | Python | Why it matters |
|---|---|---|
| `linspace(a, b, n)` = `a + (k*(b-a))/(n-1)`, ends exact | `mlinspace` | `numpy.linspace` computes `a + k*((b-a)/(n-1))` and differs in the last bit at 22 of the 60 frequencies; with `mlinspace` every grid column of the ten CSV products equals the MATLAB R2025b file |
| `round` (halves away from zero) | `mround`, `mround_array` | Python and NumPy round halves to even |
| `trapz`, `cumtrapz` operation order | `mtrapz`, `mcumtrapz` | Ref. [25] functions |
| `gradient`, `interp1` (NaN outside), `interp1(..., 'extrap')`, `conv(..., 'same')`, `find(..., 1)` | `mgradient`, `interp1`, `interp1_extrap`, `conv_same`, `first_true` | kit and G4c |
| `sprintf` spelling of NaN and Inf | `mfmt` | token strings |
| `p2_write_csv` (`%.17g`) and `p2_write_json` (2-space indent, `%d` for integers, `%.17g` otherwise) | `phase2.write_csv`, `phase2.write_json` | byte-identical output: re-writing the MATLAB R2025b files through them reproduces all ten CSV files and both JSON summaries exactly (`tests/test_writers.py`) |

## Differences that remain, and why

1. **Round-off.** The math library (`tanh`, `exp`, `pow`), LAPACK (`eig`, least squares), `polyfit` (NumPy least squares against MATLAB QR) and the FFT library differ between any two runtimes. Measured against the MATLAB R2025b products (`reference/matlab_R2025b_U1`):
   - grids are identical;
   - closed-form columns agree to about $10^{-14}$ relative;
   - eigenvalue-derived numbers agree to a few $10^{-13}$;
   - the zero sector of the $\delta_3=0$ spectrum agrees to between $7\times10^{-12}$ and $5\times10^{-11}$, depending on the runtime (Jordan blocks, $O(\sqrt{\varepsilon\lVert M\rVert})$ sensitivity).

   The fissioning $N=3.5$ GNLSE of S3 amplifies FFT rounding to about $10^{-6}$ (FINDING G-12).
2. **No watermark.** The MATLAB pack watermarked every figure of a runtime that is not canonical-eligible. SOTHE-P2 draws no watermark (PI decision, 2026-09-27). The status of the figures is recorded in `figs_manifest.json` and in `SUMMARY.md`.
3. **Sub-grid peak.** `kit.splitstep_loss` records the peak position as the vertex of the log-parabola through the grid maximum and its two neighbours. `peak="grid"` gives the MATLAB 0.1.1 definition. Both are in every G4c report (`tokens_grid_peak`, `tokens_changed_by_peak`).
4. **`noise > 0`** in `p2_splitstep_loss.m` is not ported: MATLAB's legacy `rand('seed')` stream has no NumPy equivalent. No stage uses it.
5. **Figures: publication style instead of the MATLAB layouts.** S5 and S6 plot the same data with the same symbols, marker shapes and colours as the MATLAB scripts, but in the publication style of `pubstyle.py`, which follows IOP Publishing's figure guidelines (PI instruction, 2026-09-27):
   - final widths 8.5 cm or 15 cm, Computer Modern at 8 and 9 pt;
   - PDF and EPS (vector), and PNG at 600 dpi (the MATLAB suite printed 200-dpi PNG only);
   - `captions.md` and `captions.tex` next to the figures, with explanation only (the analysis belongs to the paper text).

   The following parts of the MATLAB layouts are left out: every title; every text box and in-plot comment (fit values, spreads, region names, endpoint notes); the cyan null-locus lines of suite fig4, which lie on the frame at $\kappa_g=0$ and $d=0$; and the $\phi$ and $\psi_0$ letters and tick numbers of the Fig. 3 inset, whose range the caption gives. Four layout changes: the part labels (a), (b), (c) sit above the top-left corner of their panel, not inside it; in Fig. 4(c) of S5, the open circles ($N=1$) are drawn larger than and over the filled squares ($N=2$), so that the two points that coincide at $N\delta_3=0.05$ both show; an axis with values below $10^{-2}$ is labelled "$10^{k}$ quantity" instead of carrying an offset factor such as $\times10^{-3}$; and where a suite script leaves the axis limits automatic (fig1, fig3, fig6), the figure takes those of the corresponding `make_paper2_figs.m` 0.2.1 figure (in fig3, this puts the star at $\kappa=N$ inside the frame). Limits that a MATLAB script sets are kept, with one exception adopted by the PI on 2026-09-28 (P-4): Fig. 3 of S5 takes the vertical range $1.15\max|\mathrm{Im}\,\Omega|$ over both panels, so that the two zero-sector eigenvalues of panel (a), at $|\mathrm{Im}\,\Omega|=1.25\times10^{-3}$, lie inside the frame (MATLAB's $|\mathrm{Im}\,\Omega|\le1.15\times10^{-3}$ leaves them outside), and it draws the zero-mode sector and every non-real eigenvalue as Krein-neutral crosses with their own legend entry, since their Krein signs are round-off. MATLAB's TeX interpreter prints `\ln` literally in the legend of suite fig6; Matplotlib renders it.
6. **Bundle.** `run_phase2_all.m` saves `phase2_bundle.mat`. S6 always writes `phase2_bundle.npz`, and writes `phase2_bundle.mat` with the structs `R` and `P` when SciPy is installed.
7. **JSON of the stage reports** (`p2pack_json_write` used `jsonencode`). The content is the same; the whitespace differs.
8. **Parallelism.** Independent tasks run in a process pool, one BLAS thread per process. A parallel run gives the same numbers as a serial one.
9. **New stages.** No stage of the MATLAB pack called `run_phase2_all.m` or `run_robustness_sweep.m`. S6 and S7 run them.
10. **S3 solver.** The MATLAB stage used the MATLAB solver. SOTHE-P2 uses the Python validation port and adds run A through `ref25.src` as a cross-check.
