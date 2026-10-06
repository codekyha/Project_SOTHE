# Figures and numbers of the article

Scripts, plotted data and figure files behind the figures and the printed numbers of

> H. Oguz, *Linear spectrum, radiative loss and emission channels of solitons with third-order dispersion underlying optical analogue horizons* (2026).

The folder is part of the record `SOTHE-P2 1.0.0, SOTHE-P2-SWEEP 1.0.0 and SOTHE-P2-RAD 1.0.0: code and data for this article` (Zenodo) and of the
branch `p2` of <https://github.com/codekyha/Project_SOTHE>.

| Path | What it is |
|---|---|
| `make_article_figures.py` | draws Figs. 1 to 3, 5 and 6 and Fig. S1 of the article (PDF, EPS, PNG) and writes the plotted arrays to `figure_data.npz`; the file names `fig4_channels` and `fig5_entanglement` keep the numbering of an earlier version of the text, in which there was no figure of the steady emission: they are Figs. 5 and 6 |
| `make_article_numbers.py` | evaluates the numbers quoted in the article and writes each with its source to `article_numbers.json` |
| `zero_sector_tests.py` | the rank, perturbation, convergence, quartet, grid-scaling and scan tests of section S4 of the supplement (BdG discretization and convergence): recomputed from stage S1 of job 530047 with `sothe_p2` of SOTHE-P2 1.0.0, or read from block B5 of the sweep run and from closed forms (scan, census, stencil symbols), and a comparison of each numerical statement of that section with the recomputed value |
| `reproduce.py` | unpacks the archives of the record, runs the two scripts of the figures and numbers and compares the results with `output/`; with `--zero-sector` it also runs `zero_sector_tests.py` |
| `output/` | the six figures that these scripts draw (PDF, EPS, 600 dpi PNG), `figure_data.npz` (28 arrays), `article_numbers.json` (62 entries) and `zero_sector_tests.json` (result of `zero_sector_tests.py --grid`), as used in the article |
| `LICENSE` | MIT licence (the same text as in the three packages) |
| `MANIFEST.sha256` | SHA-256 of every file of this folder |

## Reproduce

Python 3.9 or newer with NumPy, SciPy and Matplotlib. Download from the record the two results archives (and, outside the git
layout, `SOTHE-P2-SWEEP_v1.0.0.tar.gz`), then

```
python reproduce.py --p2-results SOTHE-P2_results_20260927T175224Z.tar.gz --sweep-results SOTHE-P2-SWEEP_results_20260928T131100Z.tar.gz
python reproduce.py ... --sweep-pkg-tar SOTHE-P2-SWEEP_v1.0.0.tar.gz     # when the folder ../SOTHE-P2-SWEEP is not present
python reproduce.py ... --zero-sector                                    # also the zero-sector tests (about 2 minutes)
python reproduce.py ... --zero-sector --grid                             # and the grid scan (about 4 more minutes)
```

`reproduce.py` unpacks everything in a temporary folder, so the archives and `output/` stay untouched. To run a script by hand,
set the variables named in its header (`SOTHE_SWEEP_PKG`, `SOTHE_SWEEP_RUN`, `SOTHE_P2_FIGDATA`, `SOTHE_P2_RUN`, `SOTHE_ARTICLE_OUT`,
`SOTHE_ARTICLE_NUMBERS`).

**Checked.** With Python 3.13.16, NumPy 2.5.3, SciPy 1.18.1 and Matplotlib 3.11.2, on the three archives of the record: `reproduce.py` ends with `RESULT: PASS`; 0 of 18 figure files are byte-identical to `output/`, `figure_data.npz` has 28 arrays (identical) and `article_numbers.json` has 62 entries (byte-identical). `reproduce.py --zero-sector` also passes, and `zero_sector_tests.py --grid` finds that all 38 checks hold: the identity of its operator with the one of `sothe_p2` and the 37 groups of numerical statements of section S4 of the supplement.

## Where each figure comes from

| Figure | Data |
|---|---|
| 1 (a), (b) | SOTHE-P2-SWEEP run 20260928T131100Z, block B10: `kappa_drive_independence.csv`, `kinematic_kappa_flow.csv`; closed forms of the kinematic model |
| 2 (a) | SOTHE-P2 run 20260927T175224Z, stage S5: `figure_data.npz` (BdG spectrum about the integrable profile, delta3 = 0) |
| 2 (b) | SOTHE-P2-SWEEP block B5, task N = 1, delta3 = 0.05 (BdG spectrum about the stationary soliton, Krein signs, zero sector) |
| 3 (a) to (c) | SOTHE-P2-SWEEP block B3: `series/loss_N*_d*.csv`, `loss_points.csv` |
| 4 (a) to (d) | not drawn here: SOTHE-P2-RAD 1.0.0, `python rad.py figures` on its results archive `SOTHE-P2-RAD_results_20261005T215641Z.tar.gz` (`figures/fig_emission_recoil.pdf` and `figure_data.npz` there; the numbers of the steady emission and the recoil are in `analysis/` and `SUMMARY.md` of the same archive) |
| 5 (file `fig4_channels`) | SOTHE-P2-SWEEP block B4 (`summary.json`, window table) and the measured drifts of block B3; closed-form boost edge and threshold |
| 6 (file `fig5_entanglement`) | closed form of the logarithmic negativity at the closed-form kappa |
| S1 (a), (b) | SOTHE-P2-SWEEP block B7: `scatter_omega.csv` |

## Zero-sector tests

Section S4 of the supplement quotes results for the zero sector of the Bogoliubov-de Gennes operator $M$ about the stationary soliton $\phi$: $\dim\ker M=2$ and $\dim\ker M^2=4$ from singular values, the $\sqrt{\delta}$ splitting of the two $2\times2$ Jordan blocks under random perturbations, the null vectors of the phase and translation symmetries, the round-off floors and the growth of the spectral radius with $n_\tau$, the convergence at $\delta_3=0$, and the complex quartet that the same four modes form about the integrable profile $\psi_0$ for $\delta_3>0$ (its size, its convergence in $n_\tau$ and in the box, and the phase-matching wavenumber $k_{\rm RR}$); its last three paragraphs add the scan over $N$ and $\delta_3$, the census at the operating point and the dispersion of the two stencils. `zero_sector_tests.py` recomputes the first group with `sothe_p2` from the stationary soliton that job 530047 stored (`S1_g4a/bdg_backgrounds.npz`). The scan and census numbers it reads from block B5 of the archived sweep run (the table `B5_bdg/bdg_points.csv` and the spectrum `B5_bdg/spectra/bdg_N1_d0.05.csv`), and it computes the symbols of the two stencils and the zero crossings of the discrete branch in closed form, comparing the symbols with the matrices of `sothe_p2`. It then compares the numerical statements of the section with these values: 38 checks in all with `--grid`, the first being that its operator is the one of `sothe_p2` and each of the others covering one group of printed passages, which the script quotes (`TEX_QUOTES`; two checks need the grid scan and are skipped without `--grid`). A printed value is accepted when the recomputed value rounds to it; a printed bound is checked as a bound, at the round-off level of the arithmetic with a factor of 2 (2.5 for the null vector of the phase mode); an approximate law is accepted within the tolerance that its own line of output states. A check that holds only within the rounding of a printed value or bound, or only within the factor allowed at the round-off level, says so (a NOTE line, and `"literal": false` in the result file). Quantities at the round-off level (the four eigenvalues of the sector, the smallest singular values, the floors) depend on the BLAS library, and the supplement quotes them as bounds. The script checks numbers, not the qualitative statements or the definitions of the section, and it does not judge the argument built on the numbers.  `output/zero_sector_tests.json` is the result of the build of this folder, with the grid scan (`--grid`).

## Status of the numbers

The two archives used here hold Python output of SOTHE-P2 1.0.0 and SOTHE-P2-SWEEP 1.0.0 (NumPy 1.26.4 with MKL, UHeM Altay); the third results archive of the record, that of SOTHE-P2-RAD 1.0.0, has its own `rad.py analyse` and `rad.py figures`. The figure
scripts only plot recorded data and closed forms; they run no new simulation. `zero_sector_tests.py` is different: it recomputes the BdG operator, its spectra and, with `--grid`, the stationary soliton at other grids, with `sothe_p2`, so its results are new Python output (`output/zero_sector_tests.json` records the versions of Python and NumPy), except the scan and census numbers, which are read from block B5 of the archived sweep run. Licence: MIT.
