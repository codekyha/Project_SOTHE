"""SOTHE-P2 command line: run the Paper-2 computations and pack the results for analysis.   v1.0.0
Pure Python (NumPy + Matplotlib); UHeM Altay with SLURM, or any computer (Linux, macOS, Windows).

USAGE:  python p2.py <command> [options]          (on Altay first:  module load ANACONDA/Anaconda3-2024.06-1-python-3.12)
  install [--force]          check the MANIFEST; detect the SLURM account, the Python module and its packages;
                             write config/site.env
  probe [--local|--wait]     runtime check (versions, BLAS, numerical smoke test): here with --local (seconds, no job,
                             no core-hours), or as a SLURM job (--wait: show the job state until the verdict)
  submit [options]           the whole run as ONE SLURM job (stages, SUMMARY.md, results tarball); options:
       --stages "S1 S2"        subset of S1..S7 (default: STAGES in config/site.env)
       --fast                  smoke-test settings (short propagations, no convergence runs, small grids)
       --partition cpu2dq      SLURM partition (default SLURM_PARTITION)
       --time 0-01:00          walltime (default SLURM_WALLTIME)
       --cpus 40               cpus-per-task = parallel worker processes (default SLURM_CPUS)
       --dry-run               write the job script and print the sbatch line, submit nothing
       --local                 run the job script in this shell (inside `srun --pty bash`, or for tests)
  run [options]              the stages here, without SLURM, then the results tarball: --stages, --fast, --workers N
  stages                     list the stages S0..S7 and what each one writes
  list                       every run folder: RUN, job id, kind, state
  status [RUN]               queue (with the pending reason) + per-stage state of a run (default: the newest)
  show [RUN]                 print a run's SUMMARY.md, or a probe's probe.txt
  repack [RUN]               (re)build SOTHE-P2_results_<RUN>.tar.gz (+ .sha256) in RESULTS_DIR
  verify                     sha256 check of the pack against MANIFEST.sha256
  version | help
RUN = the RUNID, the SLURM job id, or runs/<name>.  Run folders: runs/<RUNID>, RUNID = UTC time of the start
(+ _FAST, _local); probes: runs/probe_<UTC>.  On SLURM, runs/<jobid> links to the same folder.
Every number this pack writes is RECORDED (HR-3: only MATLAB R2025b output is canonical)."""
import hashlib
import json
import os
import re
import shlex
import shutil
import socket
import subprocess
import sys
import tarfile
import time
from datetime import datetime, timezone

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")          # one BLAS thread per process; parallelism comes from worker processes
os.environ.setdefault("MPLBACKEND", "Agg")
sys.dont_write_bytecode = True

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))     # the pack folder (holds p2.py)
RUNS = os.environ.get("SOTHE_P2_RUNS") or os.path.join(ROOT, "runs")   # run folders (override: SOTHE_P2_RUNS=/scratch/...)
DEFAULTS = {"SLURM_ACCOUNT": "modfkt", "SLURM_PARTITION": "cpu2dq", "SLURM_WALLTIME": "0-01:00", "SLURM_CPUS": "40",
            "PROBE_WALLTIME": "0-00:10", "PROBE_CPUS": "2",
            "PYTHON_MODULE": "ANACONDA/Anaconda3-2024.06-1-python-3.12", "PYTHON": "python",
            "STAGES": "S1 S2 S3 S4 S5 S6 S7", "RESULTS_DIR": "/ari/users/$USER"}
# kept here (not imported from sothe_p2.stages) so that help, list, status and show need no NumPy on a login node
STAGE_DIRS = {"S1": "S1_g4a", "S2": "S2_g4c", "S3": "S3_g12", "S4": "S4_oracles", "S5": "S5_figures", "S6": "S6_phase2",
              "S7": "S7_sweep"}
STAGE_INFO = [
    ("S0", "S0_probe", "runtime probe: versions, BLAS, numerical smoke test", "probe.json"),
    ("S1", "S1_g4a", "C-05 gate G4a (run_g4a_canonical 0.1.1)", "g4a_report.json, T04_series.csv, bdg_backgrounds.npz"),
    ("S2", "S2_g4c", "C-05b T-10 + T-12, the G4C tokens (run_g4c_loss_scan 0.1.2, sub-grid peak)",
     "g4c_loss_report.json, series/*.csv, newton_solutions.npz"),
    ("S3", "S3_g12", "G-12/T-19 provenance of entropy_trajectory.csv (Ref. [25] solver, runs A-D)", "g12_report.json, traj_*.csv"),
    ("S4", "S4_oracles", "T-14 oracle ledgers v0.1.0 (29) and v0.2.0 (35)", "oracle_ledger_v010.json, oracle_ledger_v020.json, physics_pass.npz"),
    ("S5", "S5_figures", "T-13 + T-18 manuscript figures (make_paper2_figs 0.2.1)", "fig*.pdf/.eps (vector), fig*.png (600 dpi), captions.md/.tex, "
     "figs_manifest.json, figure_data.npz"),
    ("S6", "S6_phase2", "phase-2 data products (run_phase2_all v0.1.0)",
     "paper2_data/*.csv (10), phase2_summary.json, paper2_figures/ (6 figures: PDF, EPS, PNG; captions), "
     "phase2_bundle.npz/.mat, parity_vs_matlab.json"),
    ("S7", "S7_sweep", "Ref. [25] robustness sweep, 300 configurations (run_robustness_sweep)",
     "robustness_sweep.json, robustness_sweep.npz, robustness_table.csv, ref25_checks.json")]
RESULTS_PREFIX = "SOTHE-P2_results_"
FINAL_STATES = ("COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL", "OUT_OF_MEMORY", "PREEMPTED", "BOOT_FAIL", "DEADLINE")


# ============================================================ small helpers
def now(compact=False):
    t = datetime.now(timezone.utc)
    return t.strftime("%Y%m%dT%H%M%SZ" if compact else "%H:%M:%S")


def log(msg):
    print("[%s] %s" % (now(), msg), flush=True)


def warn(msg):
    print("[%s] WARNING: %s" % (now(), msg), file=sys.stderr, flush=True)


def die(msg):
    print("[%s] FATAL: %s" % (now(), msg), file=sys.stderr, flush=True)
    sys.exit(1)


def version():
    with open(os.path.join(ROOT, "VERSION")) as f:
        return f.read().strip()


def sh(cmd, timeout=120, cwd=None):
    """Run a command (list, or a string for bash -lc); return (returncode, combined output)."""
    try:
        if isinstance(cmd, str):
            cmd = ["bash", "-lc", cmd]
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True,
                           timeout=timeout, cwd=cwd)
        return p.returncode, p.stdout
    except (OSError, subprocess.TimeoutExpired) as e:
        return 127, str(e)


def load_config():
    C = dict(DEFAULTS)
    path = os.path.join(ROOT, "config", "site.env")
    if os.path.exists(path):
        with open(path, encoding="utf-8") as f:
            for ln in f:
                s = ln.strip()
                if not s or s.startswith("#") or "=" not in s:
                    continue
                k, v = s.split("=", 1)
                v = v.strip()
                if v.startswith('"') and '"' in v[1:]:
                    v = v[1:v.index('"', 1)]
                elif v.startswith("'") and "'" in v[1:]:
                    v = v[1:v.index("'", 1)]
                else:
                    v = v.split("#", 1)[0].strip()
                C[k.strip()] = v
    if os.environ.get("SOTHE_P2_RESULTS_DIR"):
        C["RESULTS_DIR"] = os.environ["SOTHE_P2_RESULTS_DIR"]
    user = os.environ.get("USER") or os.environ.get("USERNAME") or "user"
    os.environ.setdefault("USER", user)
    for k in C:
        C[k] = os.path.expanduser(os.path.expandvars(C[k]))
    if not os.path.isdir(C["RESULTS_DIR"]):
        made = False
        if os.environ.get("SOTHE_P2_RESULTS_DIR"):          # an explicit choice: create the folder
            try:
                os.makedirs(C["RESULTS_DIR"], exist_ok=True)
                made = True
            except OSError as e:
                warn("cannot create SOTHE_P2_RESULTS_DIR=%s (%s): results go to %s" % (C["RESULTS_DIR"], e, ROOT))
        if not made:                                          # e.g. the Altay default on another computer
            C["RESULTS_DIR"] = ROOT
    return C


def parse_opts(args, flags, valued):
    opts, rest = {}, []
    i = 0
    while i < len(args):
        a = args[i]
        if a in flags:
            opts[a] = True
        elif a in valued:
            if i + 1 >= len(args):
                die("option %s needs a value" % a)
            opts[a] = args[i + 1]
            i += 1
        elif a.startswith("--"):
            die("unknown option %s (see: python p2.py help)" % a)
        else:
            rest.append(a)
        i += 1
    return opts, rest


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def check_manifest(quiet=False):
    man = os.path.join(ROOT, "MANIFEST.sha256")
    if not os.path.exists(man):
        warn("no MANIFEST.sha256")
        return True
    bad = []
    with open(man, encoding="utf-8") as f:
        for ln in f:
            if not ln.strip():
                continue
            digest, name = ln.rstrip("\n").split(None, 1)
            name = name.lstrip("*").strip()
            p = os.path.join(ROOT, *name.split("/"))
            if not os.path.exists(p) or sha256(p) != digest:
                bad.append(name)
    if bad:
        for b in bad[:20]:
            warn("MANIFEST mismatch: %s" % b)
        return False
    if not quiet:
        log("MANIFEST: every file intact")
    return True


# ============================================================ runs, job ids
def run_dirs():
    if not os.path.isdir(RUNS):
        return []
    return sorted(d for d in os.listdir(RUNS) if os.path.isdir(os.path.join(RUNS, d)) and not os.path.islink(os.path.join(RUNS, d)))


def resolve_run(arg=None):
    """RUNID, job id, runs/<name> or a path -> absolute run folder; None = the newest run."""
    if not arg:
        cands = [d for d in run_dirs() if re.match(r"^20\d{6}T\d{6}Z", d)]
        return os.path.join(RUNS, cands[-1]) if cands else None
    for c in (arg, os.path.join(ROOT, arg), os.path.join(RUNS, arg)):
        if os.path.isdir(c):
            return os.path.realpath(c)
    for d in run_dirs():
        j = os.path.join(RUNS, d, "JOBID")
        if os.path.exists(j) and open(j).read().strip() == arg:
            return os.path.join(RUNS, d)
    return None


def resolve_or_die(arg):
    r = resolve_run(arg)
    if not r:
        die("no run '%s' under %s (see: python p2.py list)" % (arg or "<newest>", RUNS))
    return r


def link_jobid(job, rundir):
    link = os.path.join(RUNS, job)
    try:
        if os.path.islink(link):
            os.remove(link)
        os.symlink(os.path.basename(rundir), link)
    except OSError:
        pass                                       # job ids are also found through the JOBID files


def parse_jobid(text):
    m = re.findall(r"Submitted batch job (\d+)", text)
    return m[-1] if m else None


def squeue_state(job):
    """'STATE TIME REASON-or-NODE' of a job, '' once it has left the queue."""
    rc, out = sh(["squeue", "-h", "-j", job, "-o", "%T %M %R"], timeout=30)
    if rc != 0:
        return ""
    lines = [l for l in out.splitlines() if l.strip()]
    return lines[0].strip() if lines else ""


def start_estimate(job):
    """SLURM's expected start time of a pending job (squeue --start), or ''."""
    rc, out = sh(["squeue", "-h", "-j", job, "--start", "-o", "%S"], timeout=30)
    v = out.strip().splitlines()[-1].strip() if rc == 0 and out.strip() else ""
    return "" if v in ("", "N/A", "Unknown") else v


def queue_note(partition, cpus):
    """Before submitting: say so when the partition has no idle node (the job will then wait in the queue)."""
    if not shutil.which("sinfo"):
        return
    rc, out = sh(["sinfo", "-h", "-p", partition, "-o", "%C"], timeout=30)
    try:
        alloc, idle, other, total = [int(x) for x in out.split()[0].split("/")]
    except (ValueError, IndexError):
        return
    rc2, out2 = sh(["sinfo", "-h", "-p", partition, "-t", "idle", "-o", "%D"], timeout=30)
    idle_nodes = sum(int(x) for x in out2.split() if x.isdigit()) if rc2 == 0 else None
    rc3, out3 = sh(["squeue", "-h", "-p", partition, "-t", "PD", "-o", "%C"], timeout=30)
    pend = sum(int(x) for x in out3.split() if x.isdigit()) if rc3 == 0 else None
    if idle < int(cpus) or idle_nodes == 0:
        warn("partition %s is full now: %d of %d cores idle, %s idle nodes, %s cores already pending; "
             "the job waits in the queue until a node frees" % (partition, idle, total,
                                                               "?" if idle_nodes is None else idle_nodes,
                                                               "?" if pend is None else pend))
        print("  meanwhile, without the queue:  python p2.py probe --local   (runtime check on this node)")
        print("  or run the whole pack on any computer with Python:  python p2.py run")


def wait_job(job, file=None, pattern=None):
    """Poll every 10 s; stop when FILE matches PATTERN, the job reaches a final state or leaves the queue.
    Ctrl-C stops the waiting, not the job."""
    t0 = time.time()
    st = ""
    told = False
    try:
        while True:
            if file and pattern and os.path.exists(file) and re.search(pattern, open(file, errors="replace").read()):
                st = "finished"
                break
            st = squeue_state(job)
            if not st:
                st = "left the queue"
                break
            if st.split()[0] in FINAL_STATES:
                break
            if st.startswith("PENDING") and not told:
                est = start_estimate(job)
                print("  job %s is waiting for a free node%s" % (job, ("; SLURM expects it to start at %s" % est) if est else ""))
                told = True
            sys.stdout.write("\r  job %s: %-30s (waited %d s; Ctrl-C stops waiting, not the job)" % (job, st, time.time() - t0))
            sys.stdout.flush()
            time.sleep(10)
    except KeyboardInterrupt:
        print("\n  stopped waiting; job %s keeps running (python p2.py status)" % job)
        return
    sys.stdout.write("\r  job %s: %-30s%-60s\n" % (job, st, ""))
    time.sleep(1)


PASS_ENV = ("SOTHE_P2_FIG_FONT", "PHASE2_STAB_NTAU", "PHASE2_STAB_GRIDN", "PHASE2_STAB_GRIDD", "P2_G4C_REPORT")


def env_settings():
    """The optional settings of a run that come from the environment (recorded in RUN.env and job.sbatch)."""
    return {k: os.environ[k] for k in PASS_ENV if os.environ.get(k)}


def write_job_script(path, C, name, cpus, walltime, rundir, body):
    mod = C["PYTHON_MODULE"]
    lines = ["#!/bin/bash",
             "#SBATCH -J %s" % name,
             "#SBATCH -A %s" % C["SLURM_ACCOUNT"],
             "#SBATCH -p %s" % C["SLURM_PARTITION"],
             "#SBATCH -N 1",
             "#SBATCH -n 1",
             "#SBATCH -c %s" % cpus,
             "#SBATCH -t %s" % walltime,
             "#SBATCH -o %s/slurm-%%j.out" % rundir,
             "#SBATCH -e %s/slurm-%%j.err" % rundir,
             "# written by p2.py (SOTHE-P2 %s) at %s UTC" % (version(), now(True)),
             "set -uo pipefail",
             "if ! type module >/dev/null 2>&1; then",
             "  for f in /etc/profile.d/lmod.sh /etc/profile.d/z00_lmod.sh /usr/share/lmod/lmod/init/bash /etc/profile.d/modules.sh; do",
             "    [ -f \"$f\" ] && . \"$f\" && break",
             "  done",
             "fi",
             ("module load %s || echo \"WARNING: module load %s failed\"" % (mod, mod)) if mod and mod != "none" else "# no module",
             "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1"]
    lines += ["export %s=%s" % (k, shlex.quote(v)) for k, v in env_settings().items()]    # also when SLURM does not export
    lines += ["cd %s" % shlex.quote(ROOT), body, ""]
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines))
    os.chmod(path, 0o755)


def sbatch(script, rundir):
    if not shutil.which("sbatch"):
        die("sbatch not found: this is not a SLURM login node (use: python p2.py run)")
    rc, out = sh(["sbatch", script], timeout=60, cwd=rundir)
    for l in out.splitlines():
        print("  | " + l)
    if rc != 0:
        die("sbatch failed")
    job = parse_jobid(out)
    if not job:
        die("could not parse a job id from the sbatch output")
    with open(os.path.join(rundir, "JOBID"), "w") as f:
        f.write(job + "\n")
    link_jobid(job, rundir)
    return job


# ============================================================ commands
def cmd_install(args):
    opts, _ = parse_opts(args, ("--force",), ())
    C = load_config()
    log("SOTHE-P2 %s at %s" % (version(), ROOT))
    if not check_manifest():
        die("MANIFEST check failed: re-upload the tarball in binary mode")
    detected = {}
    log("this interpreter: %s (Python %s)" % (sys.executable, sys.version.split()[0]))
    have = {}
    for m in ("numpy", "scipy", "matplotlib"):
        try:
            have[m] = __import__(m).__version__
        except Exception:
            have[m] = None
    log("packages here: %s" % ", ".join("%s %s" % (k, v or "MISSING") for k, v in have.items()))
    slurm = bool(shutil.which("sbatch"))
    if slurm:
        rc, out = sh("module -t avail 2>&1 | grep -i '^anaconda' | sort -V", timeout=60)
        mods = [l.strip() for l in out.splitlines() if l.strip() and "/" in l]
        if mods:
            want = [m for m in mods if "python-3" in m] or mods
            detected["PYTHON_MODULE"] = want[-1] if C["PYTHON_MODULE"] not in mods else C["PYTHON_MODULE"]
            log("Anaconda modules: %s; chosen %s" % (" ".join(mods), detected["PYTHON_MODULE"]))
            rc, out = sh("module load %s >/dev/null 2>&1; python -c 'import sys, numpy, scipy, matplotlib; "
                         "print(sys.version.split()[0], numpy.__version__, scipy.__version__, matplotlib.__version__)'"
                         % detected["PYTHON_MODULE"], timeout=120)
            if rc == 0:
                v = out.strip().splitlines()[-1].split()
                log("with the module: Python %s, NumPy %s, SciPy %s, Matplotlib %s" % tuple(v[:4]))
            else:
                warn("module %s does not give numpy/scipy/matplotlib: %s" % (detected["PYTHON_MODULE"], out.strip()[-300:]))
        else:
            warn("no ANACONDA module found; set PYTHON_MODULE in config/site.env (or 'none' with PYTHON=/path/to/python)")
        rc, out = sh(["sacctmgr", "-nP", "show", "assoc", "user=%s" % os.environ["USER"], "format=account"], timeout=60)
        accts = sorted(set(l.strip() for l in out.splitlines() if l.strip())) if rc == 0 else []
        if accts:
            log("SLURM accounts of %s: %s" % (os.environ["USER"], " ".join(accts)))
            if C["SLURM_ACCOUNT"] not in accts:
                detected["SLURM_ACCOUNT"] = accts[0]
                warn("account %s not listed; using %s" % (C["SLURM_ACCOUNT"], accts[0]))
    else:
        warn("sbatch not found: no SLURM here; use  python p2.py run  (local run)")
        if sys.platform.startswith("win"):
            detected["PYTHON_MODULE"] = "none"
    src = os.path.join(ROOT, "config", "site.env.example")
    dst = os.path.join(ROOT, "config", "site.env")
    target = dst if (opts.get("--force") or not os.path.exists(dst)) else dst + ".detected"
    with open(src, encoding="utf-8") as f:
        text = f.read()
    for k, v in detected.items():
        text = re.sub(r"(?m)^(%s=)\S*" % re.escape(k), lambda m: m.group(1) + v, text)
    with open(target, "w", encoding="utf-8", newline="\n") as f:
        f.write(text)
    log("wrote %s" % os.path.relpath(target, ROOT) + ("" if target == dst else "  (config/site.env kept; --force overwrites)"))
    os.makedirs(RUNS, exist_ok=True)
    C = load_config()
    for k in ("SLURM_ACCOUNT", "SLURM_PARTITION", "SLURM_WALLTIME", "SLURM_CPUS", "PYTHON_MODULE", "STAGES", "RESULTS_DIR"):
        print("  %s=%s" % (k, C[k]))
    print("\nnext:  python p2.py probe --local    then   python p2.py submit" if slurm else "\nnext:  python p2.py run")


def _probe_body(rundir):
    """Runs inside the probe job (or here with --local): versions, BLAS, stage S0 (numerical smoke test)."""
    sys.path.insert(0, ROOT)
    log("probe %s on %s" % (os.environ.get("SLURM_JOB_ID", "local"), socket.gethostname()))
    print("--- python  %s (%s)" % (sys.version.split()[0], sys.executable))
    print("--- module  %s" % os.environ.get("LOADEDMODULES", "(none)"))
    ok = True
    try:
        from sothe_p2 import stages
        code = stages.run_stage("S0", rundir, fast=False)
        st = json.load(open(os.path.join(rundir, "S0_probe", "stage_status.json")))
        ok = code == 0 and st.get("pass")
    except Exception as e:
        warn("numerics unavailable: %s" % e)
        ok = False
    log("PROBE OK" if ok else "PROBE FAILED (see %s)" % os.path.join(rundir, "probe.txt"))


class _Tee(object):
    def __init__(self, f, stream):
        self.f, self.stream = f, stream

    def write(self, s):
        self.stream.write(s)
        self.f.write(s)
        return len(s)

    def flush(self):
        self.stream.flush()
        self.f.flush()

    def isatty(self):
        return False


def _teed(path, fn, *a):
    f = open(path, "a", encoding="utf-8")
    out, err = sys.stdout, sys.stderr
    sys.stdout, sys.stderr = _Tee(f, out), _Tee(f, err)
    try:
        return fn(*a)
    finally:
        sys.stdout.flush()
        sys.stderr.flush()
        sys.stdout, sys.stderr = out, err
        f.close()


def cmd_probe(args):
    opts, _ = parse_opts(args, ("--wait", "--local"), ("--partition",))
    C = load_config()
    C["SLURM_PARTITION"] = opts.get("--partition") or C["SLURM_PARTITION"]
    runid = "probe_" + now(True) + ("_local" if opts.get("--local") else "")
    rundir = os.path.join(RUNS, runid)
    os.makedirs(rundir, exist_ok=True)
    if opts.get("--local"):
        os.environ["P2_RUNNER"] = "p2.py probe --local"
        _teed(os.path.join(rundir, "probe.txt"), _probe_body, rundir)
        return
    script = os.path.join(rundir, "probe.sbatch")
    body = "export P2_RUNNER='p2.py probe (SLURM)'\nexec %s p2.py _probe --run %s" % (C["PYTHON"], shlex.quote(rundir))
    write_job_script(script, C, "p2py_probe", C["PROBE_CPUS"], C["PROBE_WALLTIME"], rundir, body)
    queue_note(C["SLURM_PARTITION"], C["PROBE_CPUS"])
    job = sbatch(script, rundir)
    log("probe job %s -> runs/%s/probe.txt   (also runs/%s)" % (job, runid, job))
    if opts.get("--wait"):
        pt = os.path.join(rundir, "probe.txt")
        wait_job(job, pt, r"PROBE (OK|FAILED)")
        for _ in range(6):                          # shared file systems can lag behind the job end by a few seconds
            if os.path.exists(pt) and re.search(r"PROBE (OK|FAILED)", open(pt, errors="replace").read()):
                break
            time.sleep(5)
        if os.path.exists(pt):
            for l in open(pt, errors="replace"):
                if re.search(r"PROBE|smoke|engine|versions|blas|threads|python|module|reachable|WARNING|FATAL|Error", l):
                    print("  " + l.rstrip())
        print("full output:  python p2.py show %s" % job)
    else:
        print("result:  python p2.py show %s" % job)


def _job_body(rundir, stages, fast, workers):
    sys.path.insert(0, ROOT)
    C = load_config()
    log("SOTHE-P2 run %s, job %s on %s, %d workers" % (os.path.basename(rundir), os.environ.get("SLURM_JOB_ID", "local"),
                                                    socket.gethostname(), workers))
    if fast:
        warn("FAST mode: smoke-test settings")
    from sothe_p2 import stages as ST
    codes = ST.run_all(rundir, stages, fast=fast, workers=workers)
    log("stage exit codes: %s" % codes)
    repack(rundir, C["RESULTS_DIR"])                  # the results tarball for analysis, SLURM or local
    log("done")


def cmd_job(args):
    opts, _ = parse_opts(args, ("--fast",), ("--run", "--stages", "--workers"))
    rundir = opts["--run"]
    stages = opts.get("--stages", DEFAULTS["STAGES"]).split()
    workers = int(opts.get("--workers") or os.environ.get("SLURM_CPUS_PER_TASK") or os.cpu_count() or 1)
    _teed(os.path.join(rundir, "job.log"), _job_body, rundir, stages, bool(opts.get("--fast")), workers)


def cmd_probe_job(args):
    opts, _ = parse_opts(args, (), ("--run",))
    _teed(os.path.join(opts["--run"], "probe.txt"), _probe_body, opts["--run"])


def _write_run_env(rundir, **kv):
    with open(os.path.join(rundir, "RUN.env"), "w", encoding="utf-8", newline="\n") as f:
        for k, v in kv.items():
            f.write("%s=%s\n" % (k, v))


def cmd_submit(args):
    opts, _ = parse_opts(args, ("--fast", "--dry-run", "--local"), ("--stages", "--time", "--cpus", "--partition"))
    C = load_config()
    C["SLURM_PARTITION"] = opts.get("--partition") or C["SLURM_PARTITION"]
    stages = (opts.get("--stages") or C["STAGES"]).split()
    for s in stages:
        if s not in STAGE_DIRS:
            die("unknown stage %s (S1..S7)" % s)
    fast = bool(opts.get("--fast"))
    cpus = opts.get("--cpus") or C["SLURM_CPUS"]
    walltime = opts.get("--time") or C["SLURM_WALLTIME"]
    runid = ("dryrun_" if opts.get("--dry-run") else "") + now(True) + ("_FAST" if fast else "")
    rundir = os.path.join(RUNS, runid)
    os.makedirs(rundir, exist_ok=True)
    _write_run_env(rundir, runid=runid, stages=" ".join(stages), fast=int(fast), walltime=walltime, cpus=cpus,
                   account=C["SLURM_ACCOUNT"], partition=C["SLURM_PARTITION"], python_module=C["PYTHON_MODULE"],
                   pack_version=version(), status="RECORDED (HR-3)", **env_settings())
    body = "export P2_RUNNER='p2.py submit (SLURM)'\nexec %s p2.py _job --run %s --stages %s%s --workers \"${SLURM_CPUS_PER_TASK:-%s}\"" % (
        C["PYTHON"], shlex.quote(rundir), shlex.quote(" ".join(stages)), " --fast" if fast else "", cpus)
    script = os.path.join(rundir, "job.sbatch")
    write_job_script(script, C, "p2py", cpus, walltime, rundir, body)
    if opts.get("--dry-run"):
        print("  job script: %s\n  sbatch %s   (dry run: nothing submitted)" % (script, script))
        return
    if opts.get("--local"):
        log("running the job script in this shell (no sbatch)")
        rc = subprocess.call(["bash", script], cwd=rundir)
        sys.exit(rc)
    queue_note(C["SLURM_PARTITION"], cpus)
    job = sbatch(script, rundir)
    log("submitted job %s  (run %s; stages %s; %s cpus on %s; fast %s)" % (job, runid, " ".join(stages), cpus,
                                                                         C["SLURM_PARTITION"], fast))
    st = squeue_state(job)
    if st:
        est = start_estimate(job) if st.startswith("PENDING") else ""
        print("state:      %s%s" % (st, ("   (SLURM start estimate %s)" % est) if est else ""))
    print("run folder: runs/%s   (also runs/%s)" % (runid, job))
    print("watch:      python p2.py status        log: tail -f runs/%s/job.log" % job)
    print("summary:    python p2.py show %s      (after the job ends)" % job)
    print("results:    %s/%s%s.tar.gz" % (C["RESULTS_DIR"], RESULTS_PREFIX, runid))


def cmd_run(args):
    opts, _ = parse_opts(args, ("--fast",), ("--stages", "--workers"))
    C = load_config()
    stages = (opts.get("--stages") or C["STAGES"]).split()
    fast = bool(opts.get("--fast"))
    workers = int(opts.get("--workers") or min(os.cpu_count() or 1, 40))
    runid = now(True) + ("_FAST" if fast else "") + "_local"
    rundir = os.path.join(RUNS, runid)
    os.makedirs(rundir, exist_ok=True)
    _write_run_env(rundir, runid=runid, stages=" ".join(stages), fast=int(fast), workers=workers, host=socket.gethostname(),
                   pack_version=version(), status="RECORDED (HR-3)", **env_settings())
    os.environ["P2_RUNNER"] = "p2.py run (local)"
    _teed(os.path.join(rundir, "job.log"), _job_body, rundir, stages, fast, workers)
    print("\nrun folder: %s\nsummary:    python p2.py show %s\nresults:    %s" % (
        rundir, runid, os.path.join(load_config()["RESULTS_DIR"], "%s%s.tar.gz" % (RESULTS_PREFIX, runid))))


def cmd_list(args):
    C = load_config()
    ds = run_dirs()
    if not ds:
        log("no run folders yet")
        return
    print("%-30s %-9s %-6s %s" % ("RUN (runs/<RUN>)", "JOB", "KIND", "STATE"))
    for d in ds:
        p = os.path.join(RUNS, d)
        job = open(os.path.join(p, "JOBID")).read().strip() if os.path.exists(os.path.join(p, "JOBID")) else "-"
        if d.startswith("dryrun_"):
            kind, st = "dry", "job script only (nothing submitted)"
        elif d.startswith("probe_"):
            kind = "probe"
            txt = open(os.path.join(p, "probe.txt"), errors="replace").read() if os.path.exists(os.path.join(p, "probe.txt")) else ""
            m = re.findall(r"PROBE (OK|FAILED)", txt)
            st = "PROBE " + m[-1] if m else "(no verdict yet)"
        else:
            kind = "run"
            st = "done: SUMMARY.md" if os.path.exists(os.path.join(p, "SUMMARY.md")) else "running, waiting or incomplete"
        if os.path.exists(os.path.join(C["RESULTS_DIR"], "%s%s.tar.gz" % (RESULTS_PREFIX, d))):
            st += " + tarball"
        print("%-30s %-9s %-6s %s" % (d, job, kind, st))


def cmd_status(args):
    C = load_config()
    if shutil.which("squeue"):
        rc, out = sh(["squeue", "-u", os.environ["USER"], "-o", "%.10i %.12j %.9T %.10M %.10l %R"], timeout=30)
        print(out.rstrip())
    rundir = resolve_run(args[0] if args else None)
    if not rundir:
        log("no run %s (all folders: python p2.py list)" % (args[0] if args else "yet"))
        return
    name = os.path.basename(rundir)
    job = open(os.path.join(rundir, "JOBID")).read().strip() if os.path.exists(os.path.join(rundir, "JOBID")) else "local"
    print("\nrun %s  (job %s)" % (name, job))
    if job != "local" and shutil.which("squeue"):
        qs = squeue_state(job)
        if qs:
            est = start_estimate(job) if qs.startswith("PENDING") else ""
            print("  queue: %s%s" % (qs, ("   (SLURM start estimate %s)" % est) if est else ""))
    if os.path.exists(os.path.join(rundir, "probe.txt")):
        for l in open(os.path.join(rundir, "probe.txt"), errors="replace"):
            if re.search(r"PROBE|WARNING|FATAL", l):
                print("  " + l.rstrip())
        return
    for S, d in STAGE_DIRS.items():
        p = os.path.join(rundir, d)
        stf, exf = os.path.join(p, "stage_status.json"), os.path.join(p, "stage_exit.json")
        if not os.path.isdir(p):
            print("  %-3s %-11s not started" % (S, d))
        elif os.path.exists(stf):
            st = json.load(open(stf))
            ex = json.load(open(exf)) if os.path.exists(exf) else {}
            err = (st.get("error") or "").strip().splitlines()
            print("  %-3s %-11s ran=%-5s pass=%-5s exit=%s  %.0f s  %s" % (S, d, st.get("ok"), st.get("pass"), ex.get("exit_code", "?"),
                                                                          st.get("wall_s") or 0, err[-1][:80] if err else ""))
        else:
            print("  %-3s %-11s running" % (S, d))
    if os.path.exists(os.path.join(rundir, "SUMMARY.md")):
        print("\nsummary: python p2.py show %s" % name)
    tb = os.path.join(C["RESULTS_DIR"], "%s%s.tar.gz" % (RESULTS_PREFIX, name))
    if os.path.exists(tb):
        print("tarball: %s" % tb)
    jl = os.path.join(rundir, "job.log")
    if os.path.exists(jl):
        print("")
        print("".join(open(jl, errors="replace").readlines()[-5:]).rstrip())


def cmd_show(args):
    rundir = resolve_or_die(args[0] if args else None)
    for n in ("SUMMARY.md", "probe.txt"):
        p = os.path.join(rundir, n)
        if os.path.exists(p):
            sys.stdout.write(open(p, encoding="utf-8", errors="replace").read())
            return
    log("%s: no SUMMARY.md yet" % os.path.basename(rundir))
    for n in sorted(os.listdir(rundir)):
        print("  " + n)
    jl = os.path.join(rundir, "job.log")
    if os.path.exists(jl):
        print("".join(open(jl, errors="replace").readlines()[-15:]).rstrip())


LIVE_LOGS = re.compile(r"^(job\.log|slurm-\d+\.(out|err))$")    # still written after the tarball: not in the MANIFEST


def repack(rundir, dest):
    """MANIFEST.sha256 of the run folder (every file except itself and the live logs job.log, slurm-<id>.out/.err,
    which keep growing after the tarball is written, so that the MANIFEST verifies in place and in the tarball),
    then SOTHE-P2_results_<RUN>.tar.gz of the whole folder and its .sha256."""
    rundir = os.path.realpath(rundir)
    base = os.path.basename(rundir)
    files = []
    for dp, dn, fn in os.walk(rundir):
        dn.sort()
        for n in sorted(fn):
            rel = os.path.relpath(os.path.join(dp, n), rundir).replace(os.sep, "/")
            if rel != "MANIFEST.sha256" and not LIVE_LOGS.match(rel):
                files.append(rel)
    with open(os.path.join(rundir, "MANIFEST.sha256"), "w", encoding="utf-8", newline="\n") as f:
        for rel in files:
            f.write("%s  ./%s\n" % (sha256(os.path.join(rundir, *rel.split("/"))), rel))
    name = "%s%s.tar.gz" % (RESULTS_PREFIX, base)
    out = os.path.join(dest, name)
    with tarfile.open(out, "w:gz") as tf:
        tf.add(rundir, arcname=base)
    with open(out + ".sha256", "w", encoding="utf-8", newline="\n") as f:
        f.write("%s  %s\n" % (sha256(out), name))
    log("results: %s (%.1f MB); sha256 in %s.sha256" % (out, os.path.getsize(out) / 1e6, name))
    return out


def cmd_repack(args):
    C = load_config()
    repack(resolve_or_die(args[0] if args else None), C["RESULTS_DIR"])


def cmd_verify(args):
    if not check_manifest():
        die("MANIFEST check failed")


def cmd_help(args):
    print(__doc__)


def cmd_version(args):
    print("SOTHE-P2 %s (Python %s at %s)" % (version(), sys.version.split()[0], sys.executable))


def cmd_stages(args):
    C = load_config()
    default = C["STAGES"].split()
    print("%-3s %-11s %-8s %s" % ("", "FOLDER", "DEFAULT", "WHAT / OUTPUTS"))
    for sid, d, what, outs in STAGE_INFO:
        print("%-3s %-11s %-8s %s" % (sid, d, "yes" if sid in default else ("probe" if sid == "S0" else "no"), what))
        print("%-3s %-11s %-8s   -> %s" % ("", "", "", outs))
    print("\nsubset:  python p2.py submit --stages \"S1 S2\"     (or run --stages ...)")


COMMANDS = {"install": cmd_install, "probe": cmd_probe, "submit": cmd_submit, "run": cmd_run, "list": cmd_list, "ls": cmd_list,
            "status": cmd_status, "show": cmd_show, "repack": cmd_repack, "verify": cmd_verify, "stages": cmd_stages,
            "version": cmd_version, "--version": cmd_version,
            "help": cmd_help, "-h": cmd_help, "--help": cmd_help, "_job": cmd_job, "_probe": cmd_probe_job}


NEEDS_PACK = ("install", "probe", "submit", "run", "verify", "_job", "_probe")


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    cmd = argv[0] if argv else "help"
    if cmd not in COMMANDS:
        die("unknown command %s (see: python p2.py help)" % cmd)
    if cmd in NEEDS_PACK and not os.path.exists(os.path.join(ROOT, "data", "entropy_trajectory.csv")):
        die("the pack folder is not at %s: run the commands from the unpacked SOTHE-P2 folder (python p2.py ...), "
            "or install it editable (pip install -e .)" % ROOT)
    try:
        COMMANDS[cmd](argv[1:])
    except BrokenPipeError:                         # output piped into head/less that closed early
        try:
            sys.stdout = open(os.devnull, "w")
        except OSError:
            pass


if __name__ == "__main__":
    main()
