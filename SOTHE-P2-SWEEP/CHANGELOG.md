# Changelog

## 1.0.0 (2026-09-28)

First version. One full parameter sweep of the Paper-2 computations in Python, run as one SLURM job on UHeM Altay (`cpu2dq`, 128 workers) or on any computer, built from the computational items of the referee-grade review `REVIEW_G5` of `final_G5`.

- Blocks B1–B10 (`python sweep.py blocks`): kinematic layer, entanglement layer, radiative loss and drift, comoving census with the Hawking window, BdG about the stationary TOD soliton, Raman and self-steepening (new), BdG scattering with the kinematic flow (new; D-06, INV5), Ref. [25] robustness-sweep variants, parity tables, the anchor $(N,k_{\rm op})=(1,1.5)$.
- New numerics: `sweep_p2/b_raman.py` (RK4IP with the Blow–Wood response in closed form, linear Raman, self-steepening, window tracking by exact integer rolls; checked against Gordon's rate and against the split step at $f_R=s=0$) and `sweep_p2/b_scatter.py` (sparse frequency-domain solve with exact discrete exterior modes; checked by the Kaup reflectionless limit, the no-soliton null and the Krein flux balance).
- Reused unchanged: SOTHE-P2 1.0.0 (`sothe_p2/`, vendored; `VENDORED.md`).
- Command line `sweep.py`: install, probe (with a speed calibration and wall-time estimate), submit, run, resume (local or `--submit`), blocks, list, status (tasks per block), show, repack, verify.
- Oracle ledger with three classes (hard, claim, info); exit codes 0 / 2 / 3.
- Validation: see `validation/VALIDATION.md` next to the pack tarball (not part of the tarball).
