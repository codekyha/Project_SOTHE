# Changelog of SOTHE-P2-RAD

## 1.0.0 (2026-10-06)

First release, archived in the Zenodo record of the article together with SOTHE-P2 1.0.0 and SOTHE-P2-SWEEP 1.0.0.

- `sothe_rad/solver.py`: the split-step scheme, grid and absorber of the radiative-loss runs of SOTHE-P2 1.0.0
  (`kit.splitstep_loss`), extended by a ramped third-order coefficient, a moving computational frame and the diagnostics of
  the emission analysis (tail-window flux, momentum centroid, phase rotation rate, Hamiltonian).
- `sothe_rad/model.py`: dispersion, drift-corrected resonance, local-soliton scaling, first-order emission rate.
- `sothe_rad/ode.py`: norm- and momentum-balance recoil model and its asymptotic exponent.
- `sothe_rad/analysis.py`, `sothe_rad/figures.py`: rates of the prepared solitons, the scaled rate function calibrated on
  them and tested on the launched solitons, the recoil model against the runs, window exponents, momentum balance, the loss
  table of the article reproduced, numerical checks; the article figure.
- `sothe_rad/experiments.py`: the run catalogue (43 runs).
- `sothe_rad/pubstyle.py`: copy of `sothe_p2/pubstyle.py` of SOTHE-P2 1.0.0 (figure style of the article).
- `tools/build_release.py` (release tarball and `MANIFEST.sha256`; `--doi` sets the version DOI), `tools/pack_results.py`
  (results archive), `tests/test_rad.py` (15 unit tests).
- The production run of the results archive used `rad.py`, `sothe_rad/__init__.py`, `solver.py`, `model.py`,
  `experiments.py` and `runner.py` as released here (byte-identical; their SHA-256 are listed in `CODE_SHA256.txt` of the
  results archive). `analysis.py`, `ode.py` and `figures.py` were completed after the run and produced the `analysis/`,
  `figures/` and `SUMMARY.md` of the archive from its series.
