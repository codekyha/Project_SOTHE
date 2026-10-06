# SOTHE-P2

SOTHE-P2 computes the numerical results of the second paper of the SOTHE project on solitonic event horizons in dispersive media:

- the Bogoliubov-de Gennes (BdG) spectra about the stationary third-order-dispersion (TOD) soliton;
- the radiative loss of the soliton core, from direct integration of the generalized nonlinear Schroedinger equation;
- the kinematic surface gravity, with its thermal and entanglement layers;
- the data products and the figures.

It is a Python (NumPy + Matplotlib) rewrite of the project's MATLAB pack P2_altay_pack 1.0.0. Every one of its 68 `.m` files has a Python counterpart ([docs/PORT_MAP.md](docs/PORT_MAP.md)). It runs on a laptop (Windows, Linux, macOS) or on a SLURM cluster. For UHeM Altay, see [README_UHeM_Altay.md](README_UHeM_Altay.md).

Version 1.0.0. Branch `p2` of <https://github.com/codekyha/Project_SOTHE>, tag `p2-v1.0.0`.
Zenodo DOI of this version: [10.5281/zenodo.23156910](https://doi.org/10.5281/zenodo.23156910)

> **Status of the numbers.** The reference implementation of the paper's numbers is MATLAB R2025b. This package reproduces it to round-off: grids are bit-identical, closed forms agree to about $10^{-14}$ and eigenvalues to about $10^{-13}$ (see [Agreement with MATLAB](#agreement-with-matlab)). Every report it writes marks its numbers as RECORDED Python output. The figures carry no watermark; their provenance is in the manifest next to them.

## What a run computes

A run is a sequence of stages. `python p2.py stages` lists them. Each stage writes its own folder with `stage_status.json`, `provenance.json` and `console.log`.

| Stage | What | MATLAB original | Main outputs |
|---|---|---|---|
| S0 | runtime probe: versions, BLAS, a numerical smoke test | `stage_S0_probe.m` | `probe.json` |
| S1 | gate G4a (C-05): the BdG spectrum about $\psi_0$ and about the Newton-converged TOD soliton $\phi$; the radiative loss at the operating point; Eq. (15) on the trajectory; the Gaussian layer | `run_g4a_canonical.m` 0.1.1 | `g4a_report.json`, `T04_series.csv`, `bdg_backgrounds.npz` |
| S2 | G4c (C-05b, T-10, T-12): loss scan over 9 points $(N,\delta_3)$ with convergence runs, exact-scaling check, channel census, the `[G4C-*]` token strings | `run_g4c_loss_scan.m` 0.1.2 | `g4c_loss_report.json`, `series/*.csv`, `newton_solutions.npz` |
| S3 | provenance of `data/entropy_trajectory.csv` (G-12): Ref. [25]'s solver at the archived settings (A, B) and at the settings the manuscript states (C, D) | `stage_S3_g12.m` | `g12_report.json`, `traj_*.csv` |
| S4 | oracle ledgers (T-14): v0.1.0 (29 oracles) and v0.2.0 (35) | `run_oracles.m`, `run_oracles_v020.m` | `oracle_ledger_v010.json`, `oracle_ledger_v020.json`, `physics_pass.npz` |
| S5 | the five manuscript figures (T-13, T-18) | `make_paper2_figs.m` 0.2.1 | `fig*.pdf`, `fig*.eps` (vector), `fig*.png` (600 dpi), `captions.md`/`.tex`, `figs_manifest.json`, `figure_data.npz` |
| S6 | the phase-2 data products of suite v0.1.0 | `run_phase2_all.m` | `paper2_data/` (10 CSV + `phase2_summary.json`), `paper2_figures/` (6 figures as PDF, EPS and PNG, with captions), `phase2_bundle.npz`/`.mat`, `parity_vs_matlab.json` |
| S7 | Ref. [25] robustness sweep ($5\times5\times4\times3=300$ configurations) and checks of the other Ref. [25] functions | `run_robustness_sweep.m` and SOTHE_pkg's validation protocol | `robustness_sweep.json`/`.npz`, `robustness_table.csv`, `ref25_checks.json` |

After the stages, the run folder receives:

- `SUMMARY.md`: every stage, the key numbers, the token values, the parity with MATLAB, and the open decisions;
- `summary.json`;
- `tokens_proposed.json`: the manuscript placeholder values each report supplies (nothing edits a manuscript);
- `MANIFEST.sha256`, and the tarball `SOTHE-P2_results_<RUN>.tar.gz` with its `.sha256`.

Three caveats about the outputs:

- **S6 products that are regression data, not results.** S6 reproduces suite v0.1.0 as it was. Its $\delta_3>0$ BdG rates (`radiation_rate_table.csv`, `bdg_spectrum_d3_020.csv`, `stability_map.csv`, panel (b) of fig2, fig5) linearize about $\psi_0=N\,\mathrm{sech}(N\tau)$, which is not stationary once $\delta_3\neq0$ (G3u audit, S1-1). The oracles P1.4 and P1.5 that pinned them were retired in v0.2.0.
- **The S5 figures supersede the S6 figures for the manuscript.**
- **S7 ranges are wider than Ref. [25] states.** S7 reports the sweep against the statement of Ref. [25] ($\Delta S_{\rm tot}\in[0.82,2.42]$ nats, $\eta_{\rm GSL}\in[0.52,0.62]$, every configuration with $\Delta S_{\rm tot}>0$ and $\eta_{\rm GSL}>1/2$). With `run_robustness_sweep.m` as shipped, all 300 configurations meet the two conditions, but the ranges come out wider (see the S7 section of `SUMMARY.md`).

## Figures

S5 (the five manuscript figures) and S6 (the six figures of suite v0.1.0) plot the same data with the same symbols, marker shapes and colours as the MATLAB scripts. The style follows IOP Publishing's figure guidelines (`sothe_p2/pubstyle.py`):

- **Size.** Final width 8.5 cm (one column) or 15 cm (two columns). The files are written at that size, so the lettering keeps its size in print.
- **Lettering.** Computer Modern for text and mathematics: 9 pt labels, 8 pt ticks, legends and contour labels. IOP's standard font families are Times, Helvetica, Courier and Symbol, and its graphics guide also accepts Computer Modern. `SOTHE_P2_FIG_FONT=stix` switches to the Times-like STIX fonts. Matplotlib ships both font families, so every computer draws the same glyphs.
- **Content.** No titles, text boxes or in-plot comments. Parts are labelled (a), (b), (c).
- **Colour.** Colour never carries information alone: Krein signatures also differ in marker shape, and the two curves of a panel in line style.
- **Files.** `<name>.pdf` and `<name>.eps` (vector), `<name>.png` (600 dpi).

`captions.md` and `captions.tex`, next to the figures, explain what each figure shows, with the parameters of the run. They contain no analysis; that belongs to the paper text. For S5, `figs_manifest.json` records the numbers each figure shows and `figure_data.npz` every plotted array.

## Requirements

- Python 3.9 or newer.
- NumPy 1.21 or newer and Matplotlib 3.5 or newer. Tested with NumPy 1.26.4 + Matplotlib 3.8.4 (the UHeM Anaconda module) and with NumPy 2.2.6 + Matplotlib 3.10.9.
- SciPy (optional): with it, S6 also writes `phase2_bundle.mat`.

No compiler and no MATLAB are needed.

```bash
pip install -r requirements.txt          # or: conda install numpy matplotlib scipy
```

## Quick start (any computer)

```bash
tar xzf SOTHE-P2_v1.0.0.tar.gz           # or: git clone -b p2 https://github.com/codekyha/Project_SOTHE.git SOTHE-P2
cd SOTHE-P2
python p2.py verify                      # every file against MANIFEST.sha256
python p2.py probe --local               # versions, BLAS, a smoke test (seconds)
python p2.py run                         # all stages S1..S7 with one worker per core (at most 40)
python p2.py show                        # SUMMARY.md of the newest run
```

Other ways to run:

- `python p2.py run --fast`: a smoke test. It uses short propagations, no convergence runs, and small grids. Numbers from a fast run are not the protocol values.
- `python p2.py run --stages "S6"`: one stage (or any subset).
- `python p2.py run --workers 8`: set the number of processes.

Wall time of the validation run on 2 cores (container, Python 3.12.3, NumPy 1.26.4); expect variation between machines:

| Stage | S1 | S2 | S3 | S4 | S5 | S6 | S7 | Whole run |
|---|---|---|---|---|---|---|---|---|
| Time | 22 s | 8.4 min | 16 s | 28 s | 26 s | 2.1 min | 7.3 min | 19.3 min |

With 40 workers on a SLURM node the whole run should take about 3 min.

On Windows, use `python` from an Anaconda Prompt or a venv. The commands are the same. `tar` is built into Windows 10 and later.

## Commands

`python p2.py help` prints the full list. The main commands:

| Command | What it does |
|---|---|
| `install [--force]` | checks the MANIFEST; detects the SLURM account, the Python module and its packages; writes `config/site.env` |
| `probe [--local \| --wait]` | runtime check here (`--local`) or as a SLURM job |
| `submit [--stages "S1 S2"] [--fast] [--cpus N] [--time D-HH:MM] [--partition P] [--dry-run \| --local]` | the whole run as one SLURM job |
| `run [--stages ...] [--fast] [--workers N]` | the whole run here, without SLURM |
| `stages` | the stages and their outputs |
| `list`, `status [RUN]`, `show [RUN]` | run folders, queue and stage states, the summary |
| `repack [RUN]` | rebuilds `SOTHE-P2_results_<RUN>.tar.gz` and its `.sha256` |
| `verify` | sha256 of the pack against `MANIFEST.sha256` |
| `version` | the version |

`RUN` is a RUNID, a SLURM job id, or a path.

Environment overrides:

- `SOTHE_P2_RUNS`: where run folders go (default `runs/`).
- `SOTHE_P2_RESULTS_DIR`: where results tarballs go, created if missing (default: `RESULTS_DIR` of `config/site.env` if that folder exists, else the pack folder).
- `PHASE2_STAB_NTAU`, `PHASE2_STAB_GRIDN`, `PHASE2_STAB_GRIDD`: the resolution of the S6 stability map, as in the MATLAB suite.
- `P2_G4C_REPORT`: the G4c report that Fig. 4 of S5 reads, when S2 is not part of the run.
- `SOTHE_P2_FIG_FONT=stix`: figure lettering in the Times-like STIX fonts instead of Computer Modern (see [Figures](#figures)).

A run records the ones that are set in `RUN.env`, and a SLURM job also in `job.sbatch`.

`python -m sothe_p2 <command>` is equivalent to `python p2.py <command>`. After `pip install -e .`, so is `sothe-p2 <command>`. Use an editable install only: the commands need the pack folder (`data/`, `reference/`, `config/`, `runs/`).

## Run folder

The complete list of files is in [docs/OUTPUTS.md](docs/OUTPUTS.md).

```
runs/<RUNID>/                    RUNID = UTC start time, + _FAST for --fast, + _local without SLURM
  RUN.env  JOBID  job.sbatch  job.log  slurm-<id>.out|err
  S1_g4a/  S2_g4c/  S3_g12/  S4_oracles/  S5_figures/  S6_phase2/  S7_sweep/
  SUMMARY.md  summary.json  tokens_proposed.json  MANIFEST.sha256
```

## Agreement with MATLAB

`reference/matlab_R2025b_U1/` holds the products of `run_phase2_all.m` v0.1.0 as uploaded to the project on 2026-07-18: the ten CSV files, `phase2_summary.json` (MATLAB 25.2.0.3042426, R2025b Update 1) and the six figures. S6 compares its own products with them (`parity_vs_matlab.json`). The validation run of 2026-09-27 (NumPy 1.26.4) gave:

| Quantity | Agreement |
|---|---|
| grid columns (frequencies, drives, couplings, map axes) | identical, bit for bit |
| closed-form columns (thermality, flows, grey body, entanglement, false-positive map) | $\le 3.2\times10^{-14}$ absolute |
| $\kappa_{\rm kin}$ | $8.7\times10^{-15}$ absolute ($5\times10^{-15}$ relative) |
| BdG eigenvalues and rates | $\le 2.3\times10^{-13}$ absolute; the zero sector of $\delta_3=0$ to $7\times10^{-12}$ (Jordan blocks; it varies between runtimes, up to $5\times10^{-11}$ so far) |
| oracle ledger v0.1.0 | 29/29, as in MATLAB |

The CSV and JSON writers are byte-exact ports: re-writing the MATLAB files through them reproduces the files ([tests/test_writers.py](tests/test_writers.py)). What remains is round-off of the math library, LAPACK and the FFT; [docs/PORT_MAP.md](docs/PORT_MAP.md) lists every difference. `reference/recorded_octave840/` holds the recorded Octave 8.4.0 reports of the kit (G4a 0.1.1, G4c 0.1.2), which the tests use.

## Tests

```bash
python -m unittest discover -s tests            # 41 tests, about 45 s; SOTHE_P2_QUICK=1 skips the ledgers and the suite figures (about 7 s)
```

The tests cover:

- MATLAB `linspace`, `round`, `trapz` and the other compat helpers;
- byte-identity of the writers;
- tokens, census, Gaussian layer, sub-grid peak and Kaup orders against the recorded Octave reports;
- the Ref. [25] ports, including against the validation port;
- $\kappa$, thermality and the BdG rate against MATLAB;
- both oracle ledgers;
- the figure style: scaled axes, byte-reproducible PDF, EPS and PNG files, the caption files, and the six suite figures;
- the command line (including a dry-run submission that needs no SLURM).

## Layout

```
p2.py                     command line (python p2.py help)
sothe_p2/                 the package
  cli.py                  commands, SLURM job scripts, run folders, results tarball
  stages.py               stage drivers S0..S7          runtime.py  provenance, stage records
  suite.py                suite v0.1.0 numerics          kit.py      g4a_kit numerics
  g4a.py g4c.py g12.py    stages S1 S2 S3                oracles.py  ledgers v0.1.0 / v0.2.0
  figures.py              S5 manuscript figures          suite_figs.py  S6 suite figures
  pubstyle.py             figure style (IOP guidelines), captions files
  phase2.py               S6 products and writers        postprocess.py SUMMARY.md, tokens
  compat.py               MATLAB semantics               par.py      process pool
  ref25/                  Ref. [25]: src.py (ports of the MATLAB functions), gnlse_dimensionless.py (validation port, unchanged)
data/entropy_trajectory.csv          archived GSL trajectory (input of S1, S3, S4, S6)
reference/                           MATLAB R2025b products, Octave 8.4.0 reports
config/site.env.example              site settings (python p2.py install writes config/site.env)
docs/PORT_MAP.md  docs/OUTPUTS.md    MATLAB-to-Python map; every output file
tests/                               unittest suite
tools/build_release.py               release tarball + MANIFEST (+ --doi)
tools/make_zenodo_pack.py            files and metadata of the Zenodo upload
README_UHeM_Altay.md  RELEASE.md  CHANGELOG.md  CITATION.cff  LICENSE
```

## Citation and licence

Cite the Zenodo record of this version (DOI above) and, once published, the paper it supports. `CITATION.cff` carries the metadata. The code builds on H. Oguz, *Generalized thermodynamics of solitonic event horizons in dispersive field theories*, Class. Quantum Grav. **43**, 135014 (2026), [doi:10.1088/1361-6382/ae811b](https://doi.org/10.1088/1361-6382/ae811b). Its code and data are at [doi:10.5281/zenodo.20713660](https://doi.org/10.5281/zenodo.20713660) (Ref. [25]; `sothe_p2/ref25/`, MIT).

MIT licence ([LICENSE](LICENSE)).

## Releases

[RELEASE.md](RELEASE.md) covers:

- building the release tarball;
- the Zenodo pack and its manual upload;
- the `p2` branch, the `p2-v1.0.0` tag and the GitHub release;
- where the DOI and the tag go afterwards.
