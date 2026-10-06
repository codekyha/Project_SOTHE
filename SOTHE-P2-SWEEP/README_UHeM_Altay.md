# SOTHE-P2-SWEEP on UHeM Altay: install, run, resume, repack

This page covers one complete cycle on Altay: upload the pack, install it, check the runtime, submit **one SLURM job** that runs the whole sweep, follow it, and fetch the results tarball.

Facts from earlier Altay runs (SOTHE-P2 1.0.0, job 530047; P2_pypack, job 529643):

- **Module.** `ANACONDA/Anaconda3-2024.06-1-python-3.12` gives Python 3.12.4, NumPy 1.26.4 (MKL 2023.1), SciPy 1.13.1 and Matplotlib 3.8.4.
- **Partition.** `cpu2dq` allocates a whole 128-core EPYC 7742 node and charges all 128 cores, so the sweep uses 128 workers by default.
- **Account.** Project `modfkt`. User `hoguz`, home `/ari/users/hoguz`.
- **Queue.** Job 530047 waited about 7.4 min in the queue.

## 1. Upload

Windows, WinSCP: connect to `altay.uhem.itu.edu.tr` as `hoguz`, set the transfer mode to **Binary**, and put `SOTHE-P2-SWEEP_v1.0.0.tar.gz` and `SOTHE-P2-SWEEP_v1.0.0.tar.gz.sha256` into `/ari/users/hoguz/`.

Linux or macOS:

```bash
scp SOTHE-P2-SWEEP_v1.0.0.tar.gz SOTHE-P2-SWEEP_v1.0.0.tar.gz.sha256 hoguz@altay.uhem.itu.edu.tr:/ari/users/hoguz/
```

## 2. Install (once per version)

```bash
ssh hoguz@altay.uhem.itu.edu.tr
cd /ari/users/$USER
sha256sum -c SOTHE-P2-SWEEP_v1.0.0.tar.gz.sha256        # must print: OK
tar xzf SOTHE-P2-SWEEP_v1.0.0.tar.gz                    # -> /ari/users/$USER/SOTHE-P2-SWEEP
cd SOTHE-P2-SWEEP
module load ANACONDA/Anaconda3-2024.06-1-python-3.12      # after every login
python sweep.py install                                 # MANIFEST check; module, packages, account -> config/site.env
python sweep.py probe --local                           # seconds: versions, BLAS, pool, 8 checks, wall-time estimate
```

`probe --local` ends with `PROBE OK` and prints the estimated wall time of the full sweep on 32, 64 and 128 workers (from a speed calibration on the node it runs on). `python sweep.py probe --wait` runs the same check as a SLURM job.

`install` writes `config/site.env` from `config/site.env.example`:

| Setting | Default |
|---|---|
| `SLURM_ACCOUNT` | `modfkt` |
| `SLURM_PARTITION` | `cpu2dq` |
| `SLURM_WALLTIME` | `0-01:00` |
| `SLURM_CPUS` (worker processes) | `128` |
| `PYTHON_MODULE` | `ANACONDA/Anaconda3-2024.06-1-python-3.12` |
| `BLOCKS` | `B1 B2 B3 B4 B5 B6 B7 B8 B9 B10` |
| `RESULTS_DIR` | `/ari/users/$USER` |

Edit `config/site.env` to change them; `python sweep.py install --force` rewrites it from the template. The pack keeps its own `runs/` folder and does not touch `SOTHE-P2/` or `P2_pypack/`.

## 3. Run

```bash
python sweep.py submit              # ONE cpu2dq job: phases A, B, C, SUMMARY.md, results tarball
python sweep.py status              # queue state (PENDING reason, SLURM start estimate) and tasks done per block
python sweep.py list                # every run folder: RUNID, job id, kind, state
tail -f runs/<jobid>/job.log        # live log (Ctrl-C stops watching, not the job)
python sweep.py show <jobid>        # SUMMARY.md, once the job has ended
ls -l /ari/users/$USER/SOTHE-P2-SWEEP_results_*.tar.gz*
```

The job needs no attention while it waits or runs; log out and come back with `python sweep.py list`.

Options:

```bash
python sweep.py submit --fast                   # smoke test (about 15 s of work): reduced grids, not protocol values
python sweep.py submit --blocks "B6 B7"         # a subset; dependencies are added (B2 <- B1 B4, B4 <- B3, B10 <- B1 B2 B3 B4 B5 B7)
python sweep.py submit --cpus 64 --time 0-00:30 # fewer workers, shorter walltime
python sweep.py submit --dry-run                # writes runs/dryrun_<UTC>/job.sbatch and prints the sbatch line; submits nothing
python sweep.py blocks                          # what each block computes, and its output folder
```

Interactive alternative (the same job script, in your shell on a compute node):

```bash
srun -A modfkt -p cpu2dq -N 1 -n 1 -c 128 -t 0-01:00 --pty bash
module load ANACONDA/Anaconda3-2024.06-1-python-3.12
cd /ari/users/$USER/SOTHE-P2-SWEEP && python sweep.py submit --local
```

Run folders are named by the UTC start time: `runs/<UTC>/` (`submit`), `runs/<UTC>_FAST/` (`submit --fast`), `runs/<UTC>_local/` (`run`), `runs/probe_<UTC>/`, `runs/dryrun_<UTC>/`. `runs/<jobid>` links to the same folder.

## 4. Resume (walltime, node failure, cancellation)

Every task writes its result as soon as it finishes, so an interrupted job loses only the tasks that were running.

```bash
python sweep.py status <jobid>                 # which blocks are incomplete
python sweep.py resume <jobid> --submit        # a new cpu2dq job that computes only the missing or failed tasks, then phase C
python sweep.py resume <jobid> --submit --time 0-02:00 --cpus 128
python sweep.py resume <RUN>                   # the same in this shell (another computer, or inside srun)
```

`resume` keeps the run folder, its grids (`config.json`) and its RUNID, and appends the new job id to `JOBID`. Failed tasks keep their traceback in `tasks/<block>/<id>.error.json` and are retried.

## 5. Results and repack

At the end of the job the results tarball and its checksum are written into `RESULTS_DIR`:

- `/ari/users/$USER/SOTHE-P2-SWEEP_results_<RUNID>.tar.gz`: the whole run folder with its `MANIFEST.sha256`;
- `/ari/users/$USER/SOTHE-P2-SWEEP_results_<RUNID>.tar.gz.sha256`.

```bash
python sweep.py repack               # the newest run
python sweep.py repack <jobid>       # a given run (RUNID, job id, or runs/<name>)
SWEEP_RESULTS_DIR=/path python sweep.py repack <RUN>    # the tarball somewhere else
```

`repack` rebuilds the run folder's `MANIFEST.sha256`, the tarball and its `.sha256`. Use it when the job ended before the tarball, after `resume` from another computer, or after adding files to the run folder.

Download both files with WinSCP (Binary) and check them on Windows:

```powershell
Get-FileHash -Algorithm SHA256 .\SOTHE-P2-SWEEP_results_<RUNID>.tar.gz    # compare with the .sha256 file
tar -xzf .\SOTHE-P2-SWEEP_results_<RUNID>.tar.gz                          # then open <RUNID>\SUMMARY.md
```

or on Linux or macOS:

```bash
scp 'hoguz@altay.uhem.itu.edu.tr:/ari/users/hoguz/SOTHE-P2-SWEEP_results_<RUNID>.tar.gz*' .
sha256sum -c SOTHE-P2-SWEEP_results_<RUNID>.tar.gz.sha256
tar xzf SOTHE-P2-SWEEP_results_<RUNID>.tar.gz && cd <RUNID> && sha256sum -c MANIFEST.sha256
```

Start with `SUMMARY.md` (oracle ledger and block summaries), then `ORACLES.json`, `summary.json`, the per-block CSV folders and `figures/`. The live logs `job.log` and `slurm-<id>.out/.err` are not in the MANIFEST: they are still written after the tarball.

## 6. Troubleshooting

| Symptom | Fix |
|---|---|
| `python: command not found`, or `ModuleNotFoundError: numpy` / `scipy` | `module load ANACONDA/Anaconda3-2024.06-1-python-3.12` (after every login) |
| `the pack folder is not at ...` | run the commands inside `SOTHE-P2-SWEEP/`, or as `python /ari/users/$USER/SOTHE-P2-SWEEP/sweep.py <command>` |
| `MANIFEST mismatch` at `install` or `verify` | re-upload in Binary mode, then `sha256sum -c SOTHE-P2-SWEEP_v1.0.0.tar.gz.sha256` |
| job stays `PENDING (Resources)` or `(Priority)` | the partition is full; the job starts by itself (`python sweep.py status` shows SLURM's start estimate) |
| `sbatch` prints a `BILGI:` line | UHeM's accounting notice; expected |
| `slurm-<id>.err` has `Anaconda icin bu module ilaveten ... conda.sh` | UHeM's advisory at module load; harmless |
| `Matplotlib is building the font cache` | once per account |
| `TIMEOUT` | `python sweep.py resume <jobid> --submit --time 0-02:00` |
| `OUT_OF_MEMORY` or `the worker pool broke` | `python sweep.py resume <jobid> --submit --cpus 64` (the finished tasks are kept) |
| exit code 3, `FAILED <task>` in `job.log` | read `tasks/<block>/<id>.error.json`; `resume` retries |
| exit code 2 | a hard oracle failed: `python sweep.py show <jobid>` lists it; send the results tarball for analysis |

## 7. Cost

`cpu2dq` charges all 128 cores for the wall time, whatever `--cpus` asks. The full sweep is about 9300 core-seconds (2.6 core-hours) of work in 502 tasks; the validation run took 1 h 18 min on 2 cores. On 128 workers the wall time is set by the longest tasks (the Raman propagations, about a minute each): about 3 minutes of computation plus start-up and phase C, so expect 4 to 6 minutes, i.e. roughly 9 to 13 charged core-hours (an estimate; `probe` prints the estimate for the node it runs on). `sacct -j <jobid> -o JobID,Elapsed,AllocCPUS,CPUTimeRAW` gives the billed figure. `probe --local`, `--dry-run` and runs on your own computer cost nothing.

## 8. Command reference

```
python sweep.py install [--force]        python sweep.py probe [--local|--wait]     python sweep.py submit [options]
python sweep.py run [options]            python sweep.py resume RUN [--submit]      python sweep.py blocks
python sweep.py list                     python sweep.py status [RUN]               python sweep.py show [RUN]
python sweep.py repack [RUN]             python sweep.py verify                     python sweep.py version | help
```
