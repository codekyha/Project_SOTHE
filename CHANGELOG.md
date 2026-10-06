# Changelog

## SOTHE-P2 1.0.0 (2026-09-27)

This is the first release under the name SOTHE-P2: package `sothe_p2`, branch `p2` of `codekyha/Project_SOTHE`, tag `p2-v1.0.0`. It succeeds P2_pypack 2.2.0 and replaces the MATLAB pack P2_altay_pack 1.0.0 on machines without MATLAB.

**Coverage and new stages**

- **Complete rewrite.** All 68 MATLAB files of P2_altay_pack 1.0.0 now have a Python counterpart (`docs/PORT_MAP.md`). The LaTeX tools of that pack are not included.
- **New stage S6, phase-2 data products.** A port of `run_phase2_all.m` v0.1.0:
  - the ten CSV files, written by a byte-exact port of `p2_write_csv.m`;
  - `phase2_summary.json`, written by a byte-exact port of `p2_write_json.m`;
  - the six figures of suite v0.1.0 (`suite_figs.py`, with `p2_magma`, `p2_viridis`, `p2_label_contours`);
  - `phase2_bundle.npz`, and `phase2_bundle.mat` when SciPy is installed;
  - a parity report against the MATLAB R2025b products in `reference/matlab_R2025b_U1/`.
- **New stage S7, Ref. [25] robustness sweep.** A port of `run_robustness_sweep.m` (300 configurations, process-parallel). S7 also exercises the other Ref. [25] functions through the protocol of SOTHE_pkg's own validation script:
  - the nominal run;
  - the Cherenkov phase-matching scan with `find_rr_lobe` and `rr_phase_match_fit`;
  - one `gnlse_dimensional` run.
- **`sothe_p2/ref25/src.py`.** Statement-by-statement ports of the eight Ref. [25] MATLAB functions (MATLAB `linspace`, `round`, `trapz`, `cumtrapz`; `find_rr_lobe` with its own reflect padding and 3-sigma kernel). The unchanged validation port stays the S3 solver. S3 adds run A through `src.py` as a cross-check.

**Agreement with MATLAB**

- **MATLAB grids.** `compat.mlinspace` implements MATLAB's `linspace` formula. Every grid column of the ten products is now bit-identical to the MATLAB R2025b file; `numpy.linspace` differed at 22 of the 60 frequencies.
- **Parity with the MATLAB R2025b products** (container run, NumPy 1.26.4):
  - grids are identical;
  - closed-form columns agree to $\le 3.2\times10^{-14}$;
  - eigenvalue-derived columns agree to a few $10^{-13}$ (the $\delta_3=0$ zero sector, with its Jordan blocks, to $\le 5\times10^{-11}$);
  - all 29 v0.1.0 oracles pass on the full-grid pass.

**Figures and conventions**

- **No watermarks.** The watermark of the manuscript figures (stage S5) is removed without a switch (PI decision, 2026-09-27). S5 passes when all five figures are written. Figure status and provenance are in `figs_manifest.json`.
- **Publication-grade figures** (PI instruction, 2026-09-27). S5 and S6 render in one style that follows IOP Publishing's figure guidelines (`sothe_p2/pubstyle.py`):
  - 8.5 cm or 15 cm final width, Computer Modern at 8 and 9 pt;
  - no titles, text boxes or in-plot comments; parts labelled (a), (b), (c);
  - vector PDF and EPS, and PNG at 600 dpi.

  `captions.md` and `captions.tex` explain each figure with the parameters of the run; they contain no analysis. `docs/PORT_MAP.md` (item 5) lists what the MATLAB layouts had that the figures leave out. `SOTHE_P2_FIG_FONT=stix` gives Times-like lettering.
- **Fig. 3: Krein-neutral markers and the full vertical range** (P-4, PI decision of 2026-09-28, folded into 1.0.0 before its DOI build). The zero-mode sector (the four eigenvalues of smallest modulus) and every non-real eigenvalue are drawn as Krein-neutral crosses with their own legend entry: $\sigma_3M$ is Hermitian, so $\mathrm{Im}\,\Omega\,(\lVert u\rVert^2-\lVert v\rVert^2)=0$ and the Krein signs the solver returns for them are round-off (they differed between runtimes). The vertical range is $1.15\max|\mathrm{Im}\,\Omega|$ over both panels, so every eigenvalue of the window is shown. `figures.fig3_from_data()` re-renders the figure from a run's `figure_data.npz`; the manuscript's Fig. 3 was produced that way from the Altay run `20260927T175224Z`. No number changes.
- **No environment writes.** Stages no longer write `PHASE2_STAB_*` into the process environment: S4 and S5 pass the small stability grid explicitly, and S6 uses the full grid.

**Command line**

- `stages` and `version` commands.
- `SOTHE_P2_RUNS` and `SOTHE_P2_RESULTS_DIR` overrides. The optional settings of a run (`PHASE2_STAB_*`, `P2_G4C_REPORT`, `SOTHE_P2_FIG_FONT`) are recorded in `RUN.env` and `job.sbatch`.
- Results tarball `SOTHE-P2_results_<RUN>.tar.gz`.
- Default stages S1 to S7.
- A pack-folder check.

**Release tooling and tests**

- `tools/build_release.py`: reproducible tarball, `MANIFEST.sha256`, and the `--doi` writer.
- `tools/make_zenodo_pack.py`: the files and metadata of the manual Zenodo upload.
- `tests/`: unittest suite. It includes byte-identity of the CSV and JSON writers against the MATLAB files, and census, format and layer checks against the recorded Octave 8.4.0 reports.

**Metadata and documentation**

- `README.md`, `README_UHeM_Altay.md`, `RELEASE.md`, `CITATION.cff`, `LICENSE` (MIT), `pyproject.toml` (editable install), `requirements.txt`, `.gitattributes` (LF checkout keeps the MANIFEST valid on Windows clones).

## Lineage before SOTHE-P2

### P2_pypack 2.2.0 (2026-09-27)

- **Sub-grid peak.** The peak position is the vertex of the log-parabola through the grid maximum and its neighbours (PI instruction). Reports carry `tokens_grid_peak` and `tokens_changed_by_peak`. The MATLAB kit keeps the grid maximum (decision P-3M open).
- **Run-folder MANIFEST.** The live logs (`job.log`, `slurm-*.out/err`) are no longer in the run-folder `MANIFEST.sha256`.

### P2_pypack 2.1.0 (2026-09-26)

- **Pure Python.** The MATLAB sources were removed from the pack.
- **Raw data kept for analysis**: series CSV, Newton solutions, BdG backgrounds, `figure_data.npz`.
- **Altay job 529643** (UHeM, 2026-09-26) ran this version. It reproduced every recorded Octave token.

### P2_pypack 2.0.0 (2026-09-26)

- **First Python port** of P2_altay_pack 1.0.0 for UHeM Altay, where MATLAB R2025b is not available (probe 529374).
- **Stages S0 to S5**, with the same run-folder layout as the MATLAB pack.

### P2_altay_pack 1.0.0 (2026-09-25)

- **MATLAB/SLURM pack**, 68 `.m` files:
  - counssug `g4a_kit`;
  - `sothe_phase2_matlab` v0.1.0 and the v0.2.0 overlay;
  - `make_paper2_figs` 0.2.1;
  - SOTHE_pkg v1.0.0 `src/`;
  - stage drivers S0 to S5.
- It remains the reference for a MATLAB R2025b run.
