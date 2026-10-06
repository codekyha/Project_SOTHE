# Vendored code

`sothe_p2/` is the package of **SOTHE-P2 1.0.0** (release tarball sha256 `5ea00200…e6f491ee`, delivered 2026-09-27), copied unchanged: all 24 files are byte-identical to the entries `./sothe_p2/...` of `reference/SOTHE-P2_1.0.0_MANIFEST.sha256`, the MANIFEST of that release. `python tools/build_pack.py --check` verifies it; the pack build refuses to run if a vendored file differs.

The sweep imports from it:

| Module | Used for |
|---|---|
| `sothe_p2/kit.py` | `splitstep_loss` (B3), `tod_soliton`, `bdg_spectrum_bg` (B5), `channel_census` (B4, B6), `subgrid_peak` (B6) |
| `sothe_p2/suite.py` | the kinematic layer (B1, B7, B10), `bdg_spectrum` (B5, the $\psi_0$ artefact), the thermal layer (B10), `fd_circulant` |
| `sothe_p2/compat.py` | MATLAB-compatible `linspace`, `round`, `trapz`, `interp1`, `gradient`, `conv` |
| `sothe_p2/phase2.py` | `write_csv` (the paper2 CSV format, B10) |
| `sothe_p2/pubstyle.py` | the figure style |
| `sothe_p2/ref25/src.py` | the Ref. [25] functions and the robustness-sweep configuration (B8) |

The other modules of SOTHE-P2 (stages, figures, oracles, g4a, g4c, g12, cli) are present because the package is copied whole; the sweep does not run them. `sothe_p2/cli.py` is SOTHE-P2's own command line and is not used here (`sweep.py` is this pack's).

`reference/job530047/` holds three result files of the SOTHE-P2 1.0.0 run on Altay (job 530047): `g4c_loss_report.json` (stage S2), `g4a_report.json` (stage S1) and `robustness_table.csv` (stage S7). The hard oracles of blocks B3, B4, B5 (T-01 value) and B8 compare with them.
