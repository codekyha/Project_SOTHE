# SOTHE-P2-RAD 1.0.0

Steady Cherenkov emission and recoil dynamics of solitons of the generalized nonlinear Schrödinger equation with
third-order dispersion: the emission computations of the article *Linear spectrum, radiative loss and emission channels of
solitons with third-order dispersion underlying optical analogue horizons* (H. Oguz, 2026).

Zenodo DOI of this version: [10.5281/zenodo.23156910](https://doi.org/10.5281/zenodo.23156910) (the record of the article, with
SOTHE-P2 1.0.0 and SOTHE-P2-SWEEP 1.0.0). Repository: <https://github.com/codekyha/Project_SOTHE>, branch `p2`.

## What it computes

* **Steady emission.** A soliton prepared adiabatically (the third-order coefficient raised slowly from zero) radiates at a
  constant rate, too small to resolve in the oscillating core norm. The package measures it from the resonant radiation
  that the soliton leaves behind, a plane wave whose intensity times its group velocity relative to the soliton is the
  emitted norm flux, and checks it against the slope of the core norm and of the norm left in the box.
* **Scaled rate function.** The rate of the local soliton depends on its amplitude $N$, its mean wavenumber $k_s$ and
  $\delta_3$ through the scaled coefficient $\tilde\delta=N\delta_3/(1-6\delta_3k_s)^{3/2}$. An interpolation formula
  $C\tilde\delta^{\,p}e^{g\tilde\delta-\pi q_r(\tilde\delta)}$, with $q_r$ the drift-corrected phase-matched wavenumber, is
  fitted to the rates of the prepared solitons, compared with the first-order (Born) rate, and tested on the launched solitons
  ($N=1$, $1.6$, $2$), which are not used in the fit.
* **Recoil.** Norm and momentum balance with that rate function, and no other parameter, predict the loss, the wavenumber
  shift and the velocity of launched solitons to $\xi=160$, and the exponent $c$ of $-d\ln Q/d\xi\simeq c/\xi$ over windows of
  $\xi$.
* **The loss table of the article** is reproduced (the runs of SOTHE-P2 1.0.0, same scheme, grid and absorber).

`docs/METHODS.md` gives the equations; `docs/OUTPUTS.md` the files of a results folder.

## Use

Python 3.9 or newer with NumPy, SciPy and Matplotlib (`requirements.txt`).

```bash
python rad.py verify                      # MANIFEST.sha256 of the package
python rad.py test                        # unit tests, about 15 seconds
python rad.py list                        # the 43 runs of the catalogue
python rad.py run --out results --workers 4          # all runs (about 45 core-minutes) -> results/<UTC stamp>/
python rad.py analyse --results results/<stamp>      # analysis/ and SUMMARY.md
python rad.py figures --results results/<stamp>      # figures/fig_emission_recoil.{pdf,eps,png} and figure_data.npz
python rad.py manifest --results results/<stamp>     # MANIFEST.sha256 of the results folder
```

`--set` restricts a run to groups (`table`, `prepared`, `recoil`, `scale`, `check`) or run identifiers.
`slurm/` has a job script for a SLURM cluster.

## The archived results

`SOTHE-P2-RAD_results_<stamp>.tar.gz` in the same Zenodo record is the run behind the article: the series of every run,
`analysis/`, `figures/`, `SUMMARY.md`, `RUN.json` (platform, library versions, wall times), `CODE_SHA256.txt` and
`MANIFEST.sha256`. `python rad.py verify --results <unpacked folder>` checks it; `python rad.py analyse` and
`python rad.py figures` regenerate its analysis and figure from its series. `CODE_SHA256.txt` lists the SHA-256 of the six
files of this package that ran (`rad.py`, `sothe_rad/__init__.py`, `solver.py`, `model.py`, `experiments.py`, `runner.py`);
with the archive unpacked next to the folder `SOTHE-P2-RAD`, `sha256sum -c ../<stamp>/CODE_SHA256.txt` in that folder checks
them. On the platform of `RUN.json`, the regenerated `analysis/` and `figures/` are byte-identical to the archived ones (the
name of the results folder, which `SUMMARY.md` and `summary.json` carry, aside).

## Licence and citation

MIT (`LICENSE`). Cite the Zenodo record above and the article (`CITATION.cff`).
