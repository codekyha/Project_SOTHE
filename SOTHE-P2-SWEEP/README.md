# SOTHE-P2-SWEEP 1.0.0

One full parameter sweep of the Paper-2 computations of the SOTHE project, in Python, run as **one SLURM job** on UHeM Altay (or on any computer). It turns the computational items of the referee-grade review `REVIEW_G5` into numbers: the Hawking window over $(N,k_{\rm op})$ (S1-1), the anchor $(N,k_{\rm op})=(1,1.5)$ (R-02A), the scattering solve behind the manuscript's D-06 sentence (S1-2), Raman and self-steepening in silica (S1-3), the Ref. [25] robustness ranges (W2, FINDING S7) and the parity tables (R-11).

**Status of every number: RECORDED Python output (HR-3; only MATLAB R2025b output is canonical).** Block B7 is diagnostic: scattering-derived $\alpha$ and $\beta$ are not Paper-2 claims (INV5). Nothing in this pack edits a manuscript.

## Quick start (Altay)

```bash
module load ANACONDA/Anaconda3-2024.06-1-python-3.12
python sweep.py install          # MANIFEST check; module, packages, account -> config/site.env
python sweep.py probe --local    # seconds on the login node: checks + wall-time estimate, "PROBE OK"
python sweep.py submit           # ONE cpu2dq job, 128 workers: the whole sweep, SUMMARY.md, results tarball
python sweep.py status           # progress per block;  python sweep.py show <jobid>  when it has ended
```

The results arrive as `/ari/users/$USER/SOTHE-P2-SWEEP_results_<RUNID>.tar.gz` with its `.sha256`. The full runbook (upload, resume, repack, download, cost, troubleshooting) is `README_UHeM_Altay.md`.

On any other computer: `python sweep.py run` (all cores) or `python sweep.py run --fast` (smoke test, seconds).

## What the sweep computes

| Block | What | Review items |
|---|---|---|
| B1 | kinematic layer over $N\in\{0.6,\dots,2.0\}$, $k_{\rm op}\in\{1,\dots,3\}$, $\delta_3\in\{0,0.02,0.05\}$, $\kappa_g\in\{0,0.3,0.6\}$: closed form, extraction, Eq. (7″), drive scan, $\kappa_g$ flow | S1-1 A, W11 |
| B2 | entanglement layer ($E_N$, $\nu_-$, $\bar n_{\max}$, loss-robust $E_N$) on the grid and on the Hawking window | S1-1, R-02, R-10 |
| B3 | radiative loss and drift of $N\,\mathrm{sech}(N\tau)$ under Eq. (1), $N\times\{0.02,0.05\}$ plus the nine points of job 530047, with $d\xi/2$ and $2n$ convergence runs | S1-1 A |
| B4 | comoving channel census with the measured drift epochs; Hawking window $[\max(\mu,0.05),\min(\omega^{(-)}_{\max},1.6)]$ | S1-1, R-02, R-08 |
| B5 | BdG about the stationary TOD soliton at every $(N,\delta_3)$: Krein census, zero-sector ranks ($\dim\ker M$, $\dim\ker M^2$), $\psi_0$ artefact | A-2, S1-1 A |
| B6 | Eq. (1) with Raman (Blow–Wood; linear $T_R=3$ fs) and self-steepening by RK4IP: SSFS against Gordon, drift, census, $T_0=50\dots1000$ fs | S1-3, D-4 |
| B7 | BdG scattering on the soliton-plus-flow background ($\delta_3=0$): S-matrix, conversion ratio, $n_H$ against the Planck form; Kaup and null checks | S1-2 / D-06 (INV5) |
| B8 | Ref. [25] robustness sweep: the shipped script, the production resolution, the super-Gaussian of the version of record, the draft $\delta_3$ grid, $\xi_{\max}=12$ | W2, OD-P8 |
| B9 | Supplement Tables 2 and 3 (extraction-bias parity) | R-11 |
| B10 | anchor $(1,1.5,0.05,0.3)$: the thermal and entanglement CSV files in the paper2 format, anchor summary | R-02A |

Equations, protocols and checks: `docs/PHYSICS.md`. Output files and columns: `docs/OUTPUTS.md`.

## How it runs

- One job, three phases in one process pool: A (478 independent tasks), B (24 census tasks, which need the drifts of B3), C (reductions, tables, figures, the oracle ledger, `SUMMARY.md`, the results tarball). About 9300 core-seconds in all; a few minutes on 128 workers.
- Longest task first; one task per worker process with single-threaded BLAS and FFT, so the numbers do not depend on the number of workers.
- Every task writes its result atomically to `tasks/<block>/<id>.json`. `python sweep.py resume <RUN>` recomputes only the missing or failed tasks (after a walltime limit or a node failure) and redoes phase C.
- Exit code: 0 (all tasks ran, every hard oracle passed), 2 (a hard oracle failed), 3 (a task failed).

## Oracle ledger (`ORACLES.json`)

- **hard**: identities, reproductions of recorded runs and exact checks. Examples: $\kappa$ at the operating point equals the v0.1.0 pin to $10^{-12}$; the nine loss points and the 300 Ref. [25] configurations of job 530047; Tables 2 and 3; Gordon's rate; the Kaup reflectionless limit; the Krein flux balance; the bitwise identity of the extended Ref. [25] solver. A failure means the pipeline is wrong.
- **claim**: a statement of the manuscript or of `REVIEW_G5` put to the test (for example, the Hawking window at the anchor, the SSFS rate of S1-3, the Ref. [25] ranges). A failure is a finding.
- **info**: tracked numbers with no pass rule.

## Contents

| Path | What |
|---|---|
| `sweep.py` | command line (`python sweep.py help`) |
| `sweep_p2/` | the sweep: `grid.py` (grids), `tasks.py` (task list), `runner.py` (pool, resume), `reduce.py` (phase C), `figs.py`, `summary_md.py`, `probe.py`, `cli.py`; physics blocks `b_kin.py`, `b_loss.py`, `b_census.py`, `b_bdg.py`, `b_raman.py` (new), `b_scatter.py` (new), `b_ref25.py` |
| `sothe_p2/` | SOTHE-P2 1.0.0, vendored unchanged (see `VENDORED.md`): `kit.py`, `suite.py`, `compat.py`, `phase2.py`, `pubstyle.py`, `ref25/` |
| `reference/job530047/` | the recorded Altay results the hard oracles compare with (`g4c_loss_report.json`, `g4a_report.json`, `robustness_table.csv`) |
| `config/site.env.example` | settings template (`install` writes `config/site.env`) |
| `docs/` | `PHYSICS.md`, `OUTPUTS.md` |
| `tests/test_sweep.py` | unit tests (`python tests/test_sweep.py`, about 10 s) |
| `tools/build_pack.py` | MANIFEST and reproducible tarball |
| `MANIFEST.sha256` | SHA-256 of every file (`python sweep.py verify`) |

## Requirements

Python 3.8 or later with NumPy, SciPy and Matplotlib (`requirements.txt`). Tested with the Altay module's versions (Python 3.12 module: NumPy 1.26.4, SciPy 1.13.1, Matplotlib 3.8.4; here with Python 3.11) and with NumPy 2.4.4, SciPy 1.17.1, Matplotlib 3.10.9.

## Licence and credit

MIT (`LICENSE`), © 2026 Hasan Oguz. The Ref. [25] functions in `sothe_p2/ref25/` are the author's own (SOTHE_pkg v1.0.0, MIT). Built on 2026-09-28 by Cowork for the PI.
