# SOTHE-P2-RAD on a SLURM cluster

The 43 runs of the catalogue are independent; with one worker per run the whole catalogue takes about as long as its
longest run (about 3 minutes on one core of a current CPU) plus the analysis (about a minute).

1. Copy `SOTHE-P2-RAD_v1.0.0.tar.gz` and its `.sha256` file to the cluster, check and unpack:
   `sha256sum -c SOTHE-P2-RAD_v1.0.0.tar.gz.sha256 && tar xzf SOTHE-P2-RAD_v1.0.0.tar.gz && cd SOTHE-P2-RAD`.
2. Edit `slurm/rad_altay.sbatch`: replace `ACCOUNT_NAME` and `PARTITION_NAME` by your project account and partition, and the
   module line by the Python module of your site (Python >= 3.9 with NumPy, SciPy and Matplotlib).
3. `python rad.py test` (about 15 seconds; must end with `RESULT: 15 tests, 0 failed`).
4. `sbatch slurm/rad_altay.sbatch`. The results folder `results/<UTC stamp>/` then holds the runs, `analysis/`, `figures/`,
   `SUMMARY.md` and `MANIFEST.sha256`; pack it with `python tools/pack_results.py results/<stamp>`.

The archived results were computed in one process pool of two workers on a Linux machine (`RUN.json` of the archive gives
the platform and the library versions). Values that depend on the round-off of the FFT library differ between machines in
their last digits only; the analysis reports rates to three or four significant figures.
