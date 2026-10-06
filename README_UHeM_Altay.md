# SOTHE-P2 on UHeM Altay: install, run, repack

This page covers one complete cycle on Altay: upload the pack, install it, check the runtime, submit one SLURM job, follow it, and fetch the results tarball. The whole run (stages S1 to S7) is a single job on the `cpu2dq` partition and takes a few minutes of wall time once it starts.

Facts measured on Altay (probe 529374, job 529643, 2026-09-26):

- **Module.** `ANACONDA/Anaconda3-2024.06-1-python-3.12` gives Python 3.12.4, NumPy 1.26.4 (MKL 2023.1), SciPy 1.13.1 and Matplotlib 3.8.4 on the compute nodes.
- **Partition.** `cpu2dq` allocates a whole 128-core EPYC 7742 node and charges all 128 cores.
- **Account.** Project `modfkt`.
- **Queue and run time.** Job 529643 (P2_pypack 2.1.0, stages S1 to S5) waited 7 min 43 s in the queue and ran for 1 min 43 s. SOTHE-P2 adds S6 and S7, which should take about one more minute with 40 workers (an estimate from the 2-core runs; not yet measured on Altay).
- **MATLAB.** MATLAB R2025b is not available on Altay.

## 1. Upload

On Windows, WinSCP:

1. Connect to `altay.uhem.itu.edu.tr` as `hoguz`.
2. Set Transfer settings to **Binary**.
3. Put `SOTHE-P2_v1.0.0.tar.gz` and `SOTHE-P2_v1.0.0.tar.gz.sha256` into `/ari/users/hoguz/`.

From a Linux or macOS shell:

```bash
scp SOTHE-P2_v1.0.0.tar.gz SOTHE-P2_v1.0.0.tar.gz.sha256 hoguz@altay.uhem.itu.edu.tr:/ari/users/hoguz/
```

Or clone the release on Altay, if outbound HTTPS is allowed there:

```bash
git clone -b p2 https://github.com/codekyha/Project_SOTHE.git /ari/users/$USER/SOTHE-P2
cd /ari/users/$USER/SOTHE-P2 && git checkout p2-v1.0.0
```

## 2. Install (once per version)

```bash
ssh hoguz@altay.uhem.itu.edu.tr
cd /ari/users/$USER
sha256sum -c SOTHE-P2_v1.0.0.tar.gz.sha256          # must print: OK
tar xzf SOTHE-P2_v1.0.0.tar.gz                      # -> /ari/users/$USER/SOTHE-P2
cd SOTHE-P2
module load ANACONDA/Anaconda3-2024.06-1-python-3.12  # after every login
python p2.py install                                # MANIFEST check; module, packages, account -> config/site.env
python p2.py verify                                 # MANIFEST again (after edits or a re-upload)
python p2.py probe --local                          # seconds on the login node: versions, BLAS, smoke test, "PROBE OK"
```

`install` writes `config/site.env` from `config/site.env.example`:

| Setting | Default |
|---|---|
| account | `modfkt` |
| partition | `cpu2dq` |
| walltime | `0-01:00` |
| workers | 40 |
| module | `ANACONDA/Anaconda3-2024.06-1-python-3.12` |
| stages | `S1 S2 S3 S4 S5 S6 S7` |
| results folder | `/ari/users/$USER` |

Edit `config/site.env` to change any of them. `python p2.py install --force` rewrites it from the template.

The pack keeps its own run folders. If an older pack is installed (`P2_pypack/`), leave it in place: SOTHE-P2 does not touch it.

## 3. Run

```bash
python p2.py submit                  # one cpu2dq job: S1..S7, SUMMARY.md, results tarball
python p2.py status                  # queue state (PENDING reason, SLURM start estimate) and per-stage progress
python p2.py list                    # every run folder: RUNID, job id, kind, state
tail -f runs/<jobid>/job.log         # live log (Ctrl-C stops watching, not the job)
python p2.py show <jobid>            # SUMMARY.md, once the job has ended
ls -l /ari/users/$USER/SOTHE-P2_results_*.tar.gz*
```

The job needs no attention while it waits or runs. Log out, and come back with `python p2.py list`.

Options:

```bash
python p2.py submit --fast                     # smoke test: short propagations, no convergence runs, small grids (not protocol values)
python p2.py submit --stages "S6 S7"           # a subset, e.g. only the phase-2 products and the sweep
python p2.py submit --cpus 64 --time 0-00:30   # more workers, shorter walltime
python p2.py submit --partition <name>         # another partition (ask UHeM before putting CPU work on a GPU partition)
python p2.py submit --dry-run                  # writes runs/dryrun_<UTC>/job.sbatch and prints the sbatch line; submits nothing
python p2.py probe --wait                      # the runtime check as a SLURM job (queues like any job)
```

The job script (`runs/<RUN>/job.sbatch`) loads the module, sets one BLAS thread per process and runs `python p2.py _job ...`. Parallelism comes from worker processes, one per `--cpus`.

Run folders are named by the UTC start time:

| Folder | Made by |
|---|---|
| `runs/20261001T091500Z/` | `submit` |
| `runs/20261001T091500Z_FAST/` | `submit --fast` |
| `runs/20261001T091500Z_local/` | `run` |
| `runs/probe_<UTC>/` | `probe` |
| `runs/dryrun_<UTC>/` | `submit --dry-run` |

`runs/<jobid>` links to the same folder.

## 4. Results

At the end of the job, `p2.py` writes the results tarball and its checksum into `RESULTS_DIR`:

- `/ari/users/$USER/SOTHE-P2_results_<RUNID>.tar.gz` holds the whole run folder with its `MANIFEST.sha256`;
- `/ari/users/$USER/SOTHE-P2_results_<RUNID>.tar.gz.sha256` holds its SHA-256.

Download both with WinSCP (Binary). Check the checksum on Windows:

```powershell
Get-FileHash -Algorithm SHA256 .\SOTHE-P2_results_<RUNID>.tar.gz     # compare with the .sha256 file
tar -xzf .\SOTHE-P2_results_<RUNID>.tar.gz                           # then open <RUNID>\SUMMARY.md
```

Inside the run folder, `sha256sum -c MANIFEST.sha256` (Linux) verifies every file. The live logs `job.log` and `slurm-<id>.out/.err` are not in the MANIFEST: they are still written after the tarball.

## 5. Repack

```bash
python p2.py repack               # the newest run
python p2.py repack <jobid>       # a given run (RUNID, job id, or runs/<name>)
```

`repack` rebuilds the run folder's `MANIFEST.sha256`, `SOTHE-P2_results_<RUNID>.tar.gz` and its `.sha256` in `RESULTS_DIR`. Use it in three cases:

- the job ended before the tarball (walltime, cancellation);
- you added files to the run folder (for example an analysis);
- you want the tarball somewhere else: `SOTHE_P2_RESULTS_DIR=/path python p2.py repack <RUN>`.

To re-run one stage into a fresh folder, submit the subset (`--stages "S2"`). Stages never overwrite another run's folder.

## 6. Troubleshooting

| Symptom | Fix |
|---|---|
| `python: command not found`, or `ModuleNotFoundError: numpy` | `module load ANACONDA/Anaconda3-2024.06-1-python-3.12` (after every login) |
| `the pack folder is not at ...` | run the commands inside `SOTHE-P2/`, or as `python /ari/users/$USER/SOTHE-P2/p2.py <command>` |
| `MANIFEST mismatch` at `install` or `verify` | re-upload the tarball in Binary mode, then check `sha256sum -c SOTHE-P2_v1.0.0.tar.gz.sha256` |
| job stays `PENDING (Resources)` or `(Priority)` | the partition is full; the job starts by itself. `python p2.py status` shows SLURM's start estimate. Meanwhile `python p2.py probe --local` works on the login node |
| `sbatch` prints a `BILGI:` line | UHeM's accounting notice; expected |
| `slurm-<id>.err` has `Anaconda icin bu module ilaveten ... conda.sh` | UHeM's advisory at module load; harmless (the packages come from the module, see `provenance.json`) |
| `Matplotlib is building the font cache` | once per account; later runs are faster |
| a stage shows `ran=False` | `python p2.py show <jobid>`, then `runs/<RUN>/<stage>/console.log`; the other stages still ran |
| `TIMEOUT` | `python p2.py submit --time 0-02:00` (the full run needs minutes; a timeout points to a stalled node) |
| `OUT_OF_MEMORY` | fewer workers: `--cpus 32` |

## 7. Cost

`cpu2dq` charges all 128 cores of the node for the wall time, whatever `--cpus` asks. Job 529643 (S1 to S5, 1 min 43 s) cost about 3.7 core-hours. A full SOTHE-P2 run should take about 3 min, which is about 6.5 core-hours (estimate); `--fast` costs less. `sacct -j <jobid> -o JobID,Elapsed,AllocCPUS,CPUTimeRAW` gives the billed figure. `probe --local` and runs on your own computer (`python p2.py run`) cost nothing.

## 8. Command reference

```
python p2.py install [--force]      python p2.py probe [--local|--wait]      python p2.py submit [options]
python p2.py run [options]          python p2.py stages                      python p2.py list
python p2.py status [RUN]           python p2.py show [RUN]                  python p2.py repack [RUN]
python p2.py verify                 python p2.py version                     python p2.py help
```
