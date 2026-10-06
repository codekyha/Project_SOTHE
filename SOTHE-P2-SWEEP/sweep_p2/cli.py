"""SOTHE-P2-SWEEP command line: one full parameter sweep of the Paper-2 computations in Python, as ONE SLURM job.   v1.0.0
Pure Python (NumPy, SciPy, Matplotlib); UHeM Altay with SLURM, or any computer (Linux, macOS, Windows).

USAGE:  python sweep.py <command> [options]      (on Altay first:  module load ANACONDA/Anaconda3-2024.06-1-python-3.12)
  install [--force]          check the MANIFEST; detect the SLURM account, the Python module and its packages;
                             write config/site.env
  probe [--local|--wait]     runtime check (versions, BLAS, process pool, numerical smoke test of every block, speed
                             calibration): here with --local (about a minute, no job), or as a SLURM job (--wait: follow it)
  submit [options]           the whole sweep as ONE SLURM job (all blocks, SUMMARY.md, results tarball); options:
       --blocks "B3 B7"        subset of B1..B10 (dependencies are added; default: BLOCKS in config/site.env)
       --fast                  smoke-test grids (seconds to minutes; nothing is compared with job 530047)
       --partition cpu2dq      SLURM partition (default SLURM_PARTITION)
       --time 0-01:00          walltime (default SLURM_WALLTIME)
       --cpus 128              cpus-per-task = worker processes (default SLURM_CPUS)
       --dry-run               write the job script and print the sbatch line, submit nothing
       --local                 run the job script in this shell (inside `srun --pty bash`, or for tests)
  run [options]              the sweep here, without SLURM, then the results tarball: --blocks, --fast, --workers N
  resume RUN [options]       finish an interrupted run (walltime, node failure): only the tasks without a result are
                             computed, then the reductions are redone; here by default, or --submit [--cpus --time
                             --partition] as a SLURM job; --workers N
  blocks                     list the blocks B1..B10 and what each one computes
  list                       every run folder: RUN, job id, kind, state
  status [RUN]               queue (with the pending reason) + per-block progress of a run (default: the newest)
  show [RUN]                 print a run's SUMMARY.md, or a probe's probe.txt
  repack [RUN]               (re)build SOTHE-P2-SWEEP_results_<RUN>.tar.gz (+ .sha256) in RESULTS_DIR
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

for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
    os.environ.setdefault(_v, "1")          # one BLAS/FFT thread per process; parallelism comes from worker processes
os.environ.setdefault("MPLBACKEND", "Agg")
sys.dont_write_bytecode = True

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))     # the pack folder (holds sweep.py)
RUNS = os.environ.get("SWEEP_RUNS") or os.path.join(ROOT, "runs")      # run folders (override: SWEEP_RUNS=/scratch/...)
DEFAULTS = {"SLURM_ACCOUNT": "modfkt", "SLURM_PARTITION": "cpu2dq", "SLURM_WALLTIME": "0-01:00", "SLURM_CPUS": "128",
            "PROBE_WALLTIME": "0-00:10", "PROBE_CPUS": "4",
            "PYTHON_MODULE": "ANACONDA/Anaconda3-2024.06-1-python-3.12", "PYTHON": "python",
            "BLOCKS": "B1 B2 B3 B4 B5 B6 B7 B8 B9 B10", "RESULTS_DIR": "/ari/users/$USER"}
# kept here (not imported from sweep_p2.grid) so that help, list, status and show need no NumPy on a login node
BLOCK_INFO = [
    ("B1", "kinematic layer over (N, k_op, delta3, kappa_g): closed form, extraction, drive scan, kappa_g flow", "B1_kinematics/"),
    ("B2", "entanglement layer on the grid and on the Hawking window (closed form)", "B2_entanglement/"),
    ("B3", "radiative loss and drift of the launched soliton (split step; dxi/2 and 2n convergence runs)", "B3_loss/"),
    ("B4", "comoving channel census with the measured drift epochs; Hawking window, escape threshold", "B4_census/"),
    ("B5", "BdG about the stationary TOD soliton: spectrum, Krein census, zero-sector ranks, psi0 artefact", "B5_bdg/"),
    ("B6", "Raman (Blow-Wood, linear T_R) and self-steepening, RK4IP: SSFS vs Gordon, drift, census", "B6_raman/"),
    ("B7", "mode-conversion scattering with the kinematic flow (delta3 = 0): S-matrix, n_H vs thermal (D-06)", "B7_scatter/"),
    ("B8", "Ref. [25] robustness sweep: shipped script, production resolution, VoR super-Gaussian, draft grid, xi 12", "B8_ref25/"),
    ("B9", "extraction-bias parity tables (Supplement Tables 2 and 3)", "B9_parity/"),
    ("B10", "anchor (N, k_op) = (1, 1.5): paper2-format CSV files, anchor summary", "B10_anchor/")]
ALL_BLOCKS = [b[0] for b in BLOCK_INFO]
RESULTS_PREFIX = "SOTHE-P2-SWEEP_results_"
FINAL_STATES = ("COMPLETED", "FAILED", "CANCELLED", "TIMEOUT", "NODE_FAIL", "OUT_OF_MEMORY", "PREEMPTED", "BOOT_FAIL", "DEADLINE")
PASS_ENV = ("SOTHE_P2_FIG_FONT",)


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
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, timeout=timeout, cwd=cwd)
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
    if os.environ.get("SWEEP_RESULTS_DIR"):
        C["RESULTS_DIR"] = os.environ["SWEEP_RESULTS_DIR"]
    user = os.environ.get("USER") or os.environ.get("USERNAME") or "user"
    os.environ.setdefault("USER", user)
    for k in C:
        C[k] = os.path.expanduser(os.path.expandvars(C[k]))
    if not os.path.isdir(C["RESULTS_DIR"]):
        made = False
        if os.environ.get("SWEEP_RESULTS_DIR"):
            try:
                os.makedirs(C["RESULTS_DIR"], exist_ok=True)
                made = True
            except OSError as e:
                warn("cannot create SWEEP_RESULTS_DIR=%s (%s): results go to %s" % (C["RESULTS_DIR"], e, ROOT))
        if not made:                                 # e.g. the Altay default on another computer
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
            die("unknown option %s (see: python sweep.py help)" % a)
        else:
            rest.append(a)
        i += 1
    return opts, rest


def parse_blocks(text):
    bl = [b.strip().upper() for b in re.split(r"[\s,]+", text or "") if b.strip()]
    for b in bl:
        if b not in ALL_BLOCKS:
            die("unknown block %s (B1..B10; see: python sweep.py blocks)" % b)
    return bl


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
    bad, n = [], 0
    with open(man, encoding="utf-8") as f:
        for ln in f:
            if not ln.strip():
                continue
            digest, name = ln.rstrip("\n").split(None, 1)
            name = name.lstrip("*").strip()
            if name.startswith("./"):
                name = name[2:]
            p = os.path.join(ROOT, *name.split("/"))
            n += 1
            if not os.path.exists(p) or sha256(p) != digest:
                bad.append(name)
    if bad:
        for b in bad[:20]:
            warn("MANIFEST mismatch: %s" % b)
        return False
    if not quiet:
        log("MANIFEST: all %d files intact" % n)
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
        if os.path.exists(j) and arg in open(j).read().split():
            return os.path.join(RUNS, d)
    return None


def resolve_or_die(arg):
    r = resolve_run(arg)
    if not r:
        die("no run '%s' under %s (see: python sweep.py list)" % (arg or "<newest>", RUNS))
    return r


def link_jobid(job, rundir):
    link = os.path.join(RUNS, job)
    try:
        if os.path.islink(link):
            os.remove(link)
        os.symlink(os.path.basename(rundir), link)
    except OSError:
        pass


def parse_jobid(text):
    m = re.findall(r"Submitted batch job (\d+)", text)
    return m[-1] if m else None


def squeue_state(job):
    rc, out = sh(["squeue", "-h", "-j", job, "-o", "%T %M %R"], timeout=30)
    if rc != 0:
        return ""
    lines = [l for l in out.splitlines() if l.strip()]
    return lines[0].strip() if lines else ""


def start_estimate(job):
    rc, out = sh(["squeue", "-h", "-j", job, "--start", "-o", "%S"], timeout=30)
    v = out.strip().splitlines()[-1].strip() if rc == 0 and out.strip() else ""
    return "" if v in ("", "N/A", "Unknown") else v


def queue_note(partition, cpus):
    if not shutil.which("sinfo"):
        return
    rc, out = sh(["sinfo", "-h", "-p", partition, "-o", "%C"], timeout=30)
    try:
        alloc, idle, other, total = [int(x) for x in out.split()[0].split("/")]
    except (ValueError, IndexError):
        return
    rc2, out2 = sh(["sinfo", "-h", "-p", partition, "-t", "idle", "-o", "%D"], timeout=30)
    idle_nodes = sum(int(x) for x in out2.split() if x.isdigit()) if rc2 == 0 else None
    if idle < int(cpus) or idle_nodes == 0:
        warn("partition %s is full now: %d of %d cores idle, %s idle nodes; the job waits in the queue until a node frees" % (
            partition, idle, total, "?" if idle_nodes is None else idle_nodes))
        print("  meanwhile, without the queue:  python sweep.py probe --local   (runtime check on this node)")


def wait_job(job, file=None, pattern=None):
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
        print("\n  stopped waiting; job %s keeps running (python sweep.py status)" % job)
        return
    sys.stdout.write("\r  job %s: %-30s%-60s\n" % (job, st, ""))
    time.sleep(1)


def env_settings():
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
             "# written by sweep.py (SOTHE-P2-SWEEP %s) at %s UTC" % (version(), now(True)),
             "set -uo pipefail",
             "if ! type module >/dev/null 2>&1; then",
             "  for f in /etc/profile.d/lmod.sh /etc/profile.d/z00_lmod.sh /usr/share/lmod/lmod/init/bash /etc/profile.d/modules.sh; do",
             "    [ -f \"$f\" ] && . \"$f\" && break",
             "  done",
             "fi",
             ("module load %s || echo \"WARNING: module load %s failed\"" % (mod, mod)) if mod and mod != "none" else "# no module",
             "export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 NUMEXPR_NUM_THREADS=1 VECLIB_MAXIMUM_THREADS=1",
             "export MPLBACKEND=Agg PYTHONDONTWRITEBYTECODE=1"]
    lines += ["export %s=%s" % (k, shlex.quote(v)) for k, v in env_settings().items()]
    lines += ["cd %s" % shlex.quote(ROOT), body, ""]
    with open(path, "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(lines))
    os.chmod(path, 0o755)


def sbatch(script, rundir):
    if not shutil.which("sbatch"):
        die("sbatch not found: this is not a SLURM login node (use: python sweep.py run)")
    rc, out = sh(["sbatch", script], timeout=60, cwd=rundir)
    for l in out.splitlines():
        print("  | " + l)
    if rc != 0:
        die("sbatch failed")
    job = parse_jobid(out)
    if not job:
        die("could not parse a job id from the sbatch output")
    with open(os.path.join(rundir, "JOBID"), "a") as f:
        f.write(job + "\n")
    link_jobid(job, rundir)
    return job


def _write_run_env(rundir, **kv):
    with open(os.path.join(rundir, "RUN.env"), "a", encoding="utf-8", newline="\n") as f:
        f.write("# %s UTC\n" % now(True))
        for k, v in kv.items():
            f.write("%s=%s\n" % (k, v))


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


# ============================================================ commands
def cmd_install(args):
    opts, _ = parse_opts(args, ("--force",), ())
    C = load_config()
    log("SOTHE-P2-SWEEP %s at %s" % (version(), ROOT))
    if not check_manifest():
        die("MANIFEST check failed: re-upload the tarball in binary mode, or unpack it again")
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
                         % detected["PYTHON_MODULE"], timeout=180)
            if rc == 0 and out.strip():
                v = out.strip().splitlines()[-1].split()
                log("with the module: Python %s, NumPy %s, SciPy %s, Matplotlib %s" % tuple((v + ["?"] * 4)[:4]))
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
        warn("sbatch not found: no SLURM here; use  python sweep.py run  (local run)")
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
    for k in ("SLURM_ACCOUNT", "SLURM_PARTITION", "SLURM_WALLTIME", "SLURM_CPUS", "PYTHON_MODULE", "BLOCKS", "RESULTS_DIR"):
        print("  %s=%s" % (k, C[k]))
    print("\nnext:  python sweep.py probe --local    then   python sweep.py submit" if slurm else "\nnext:  python sweep.py probe --local   then   python sweep.py run")


def _probe_body(rundir):
    sys.path.insert(0, ROOT)
    log("probe %s on %s" % (os.environ.get("SLURM_JOB_ID", "local"), socket.gethostname()))
    print("--- python  %s (%s)" % (sys.version.split()[0], sys.executable))
    print("--- module  %s" % os.environ.get("LOADEDMODULES", "(none)"))
    ok = False
    try:
        from sweep_p2 import probe
        ok = probe.run(rundir, workers=int(os.environ.get("SLURM_CPUS_PER_TASK") or 2))
    except Exception as e:
        import traceback
        traceback.print_exc()
        warn("numerics unavailable: %s" % e)
    log("PROBE OK" if ok else "PROBE FAILED (see %s)" % os.path.join(rundir, "probe.txt"))


def cmd_probe(args):
    opts, _ = parse_opts(args, ("--wait", "--local"), ("--partition",))
    C = load_config()
    C["SLURM_PARTITION"] = opts.get("--partition") or C["SLURM_PARTITION"]
    runid = "probe_" + now(True) + ("_local" if opts.get("--local") else "")
    rundir = os.path.join(RUNS, runid)
    os.makedirs(rundir, exist_ok=True)
    if opts.get("--local"):
        _teed(os.path.join(rundir, "probe.txt"), _probe_body, rundir)
        return
    script = os.path.join(rundir, "probe.sbatch")
    body = "exec %s sweep.py _probe --run %s" % (C["PYTHON"], shlex.quote(rundir))
    write_job_script(script, C, "sweep_probe", C["PROBE_CPUS"], C["PROBE_WALLTIME"], rundir, body)
    queue_note(C["SLURM_PARTITION"], C["PROBE_CPUS"])
    job = sbatch(script, rundir)
    log("probe job %s -> runs/%s/probe.txt   (also runs/%s)" % (job, runid, job))
    if opts.get("--wait"):
        pt = os.path.join(rundir, "probe.txt")
        wait_job(job, pt, r"PROBE (OK|FAILED)")
        for _ in range(6):
            if os.path.exists(pt) and re.search(r"PROBE (OK|FAILED)", open(pt, errors="replace").read()):
                break
            time.sleep(5)
        if os.path.exists(pt):
            for l in open(pt, errors="replace"):
                if re.search(r"PROBE|check|estimate|versions|blas|python|module|WARNING|FATAL|Error", l):
                    print("  " + l.rstrip())
        print("full output:  python sweep.py show %s" % job)
    else:
        print("result:  python sweep.py show %s" % job)


def _job_body(rundir, blocks, fast, workers, resume):
    sys.path.insert(0, ROOT)
    C = load_config()
    log("SOTHE-P2-SWEEP run %s, job %s on %s, %d workers%s" % (os.path.basename(rundir), os.environ.get("SLURM_JOB_ID", "local"),
                                                              socket.gethostname(), workers, " (resume)" if resume else ""))
    if fast:
        warn("FAST mode: smoke-test grids")
    from sweep_p2 import runner
    code = 1
    try:
        code = runner.run_sweep(rundir, fast=fast, blocks=blocks, workers=workers, resume=resume)
    except Exception:
        import traceback
        traceback.print_exc()
        warn("the sweep stopped with an exception; finished tasks are kept (python sweep.py resume %s)" % os.path.basename(rundir))
    with open(os.path.join(rundir, "EXIT_CODE"), "w") as f:
        f.write("%d\n" % code)
    repack(rundir, C["RESULTS_DIR"])
    log("done (exit code %d: %s)" % (code, {0: "all tasks ran, every hard oracle passed", 2: "a hard oracle failed or a reduction broke",
                                             3: "some tasks failed (resume retries them)"}.get(code, "stopped early")))
    return code


def cmd_job(args):
    opts, _ = parse_opts(args, ("--fast", "--resume"), ("--run", "--blocks", "--workers"))
    rundir = opts["--run"]
    blocks = parse_blocks(opts.get("--blocks", " ".join(ALL_BLOCKS)))
    workers = int(opts.get("--workers") or os.environ.get("SLURM_CPUS_PER_TASK") or os.cpu_count() or 1)
    code = _teed(os.path.join(rundir, "job.log"), _job_body, rundir, blocks, bool(opts.get("--fast")), workers, bool(opts.get("--resume")))
    sys.exit(code)


def cmd_probe_job(args):
    opts, _ = parse_opts(args, (), ("--run",))
    _teed(os.path.join(opts["--run"], "probe.txt"), _probe_body, opts["--run"])


def _submit(C, rundir, runid, blocks, fast, cpus, walltime, resume, dry, local, target=None):
    target = target or rundir
    body = "exec %s sweep.py _job --run %s --blocks %s%s%s --workers \"${SLURM_CPUS_PER_TASK:-%s}\"" % (
        C["PYTHON"], shlex.quote(target), shlex.quote(" ".join(blocks)), " --fast" if fast else "", " --resume" if resume else "", cpus)
    script = os.path.join(rundir, "resume.sbatch" if resume else "job.sbatch")
    write_job_script(script, C, "sweep_resume" if resume else "sweep", cpus, walltime, target, body)
    if dry:
        print("  job script: %s\n  sbatch %s   (dry run: nothing submitted)" % (script, script))
        return
    if local:
        log("running the job script in this shell (no sbatch)")
        rc = subprocess.call(["bash", script], cwd=rundir)
        sys.exit(rc)
    queue_note(C["SLURM_PARTITION"], cpus)
    job = sbatch(script, rundir)
    log("submitted job %s  (run %s; blocks %s; %s cpus on %s; fast %s%s)" % (job, runid, " ".join(blocks), cpus, C["SLURM_PARTITION"], fast,
                                                                           "; resume" if resume else ""))
    st = squeue_state(job)
    if st:
        est = start_estimate(job) if st.startswith("PENDING") else ""
        print("state:      %s%s" % (st, ("   (SLURM start estimate %s)" % est) if est else ""))
    print("run folder: runs/%s   (also runs/%s)" % (runid, job))
    print("watch:      python sweep.py status        log: tail -f runs/%s/job.log" % job)
    print("summary:    python sweep.py show %s      (after the job ends)" % job)
    print("results:    %s/%s%s.tar.gz" % (C["RESULTS_DIR"], RESULTS_PREFIX, runid))


def cmd_submit(args):
    opts, _ = parse_opts(args, ("--fast", "--dry-run", "--local"), ("--blocks", "--time", "--cpus", "--partition"))
    C = load_config()
    C["SLURM_PARTITION"] = opts.get("--partition") or C["SLURM_PARTITION"]
    blocks = parse_blocks(opts.get("--blocks") or C["BLOCKS"])
    fast = bool(opts.get("--fast"))
    cpus = opts.get("--cpus") or C["SLURM_CPUS"]
    walltime = opts.get("--time") or C["SLURM_WALLTIME"]
    runid = ("dryrun_" if opts.get("--dry-run") else "") + now(True) + ("_FAST" if fast else "")
    rundir = os.path.join(RUNS, runid)
    os.makedirs(rundir, exist_ok=True)
    _write_run_env(rundir, runid=runid, blocks=" ".join(blocks), fast=int(fast), walltime=walltime, cpus=cpus, account=C["SLURM_ACCOUNT"],
                   partition=C["SLURM_PARTITION"], python_module=C["PYTHON_MODULE"], pack_version=version(), status="RECORDED (HR-3)",
                   **env_settings())
    _submit(C, rundir, runid, blocks, fast, cpus, walltime, False, opts.get("--dry-run"), opts.get("--local"))


def cmd_run(args):
    opts, _ = parse_opts(args, ("--fast",), ("--blocks", "--workers"))
    C = load_config()
    blocks = parse_blocks(opts.get("--blocks") or C["BLOCKS"])
    fast = bool(opts.get("--fast"))
    workers = int(opts.get("--workers") or max(1, min(os.cpu_count() or 1, 128)))
    runid = now(True) + ("_FAST" if fast else "") + "_local"
    rundir = os.path.join(RUNS, runid)
    os.makedirs(rundir, exist_ok=True)
    _write_run_env(rundir, runid=runid, blocks=" ".join(blocks), fast=int(fast), workers=workers, host=socket.gethostname(),
                   pack_version=version(), status="RECORDED (HR-3)", **env_settings())
    code = _teed(os.path.join(rundir, "job.log"), _job_body, rundir, blocks, fast, workers, False)
    print("\nrun folder: %s\nsummary:    python sweep.py show %s\nresults:    %s" % (
        rundir, runid, os.path.join(load_config()["RESULTS_DIR"], "%s%s.tar.gz" % (RESULTS_PREFIX, runid))))
    sys.exit(code)


def cmd_resume(args):
    opts, rest = parse_opts(args, ("--submit", "--dry-run"), ("--workers", "--cpus", "--time", "--partition"))
    if not rest:
        die("resume needs a RUN (see: python sweep.py list)")
    rundir = resolve_or_die(rest[0])
    cfgp = os.path.join(rundir, "config.json")
    if not os.path.exists(cfgp):
        die("%s has no config.json: the run never started; submit a new one" % rundir)
    with open(cfgp, encoding="utf-8") as f:
        stored = json.load(f)
    blocks, fast = stored["blocks"], bool(stored["cfg"]["fast"])
    C = load_config()
    runid = os.path.basename(rundir)
    if opts.get("--submit") or opts.get("--dry-run"):
        C["SLURM_PARTITION"] = opts.get("--partition") or C["SLURM_PARTITION"]
        cpus = opts.get("--cpus") or C["SLURM_CPUS"]
        walltime = opts.get("--time") or C["SLURM_WALLTIME"]
        if opts.get("--dry-run"):                  # leaves no trace in the run folder
            import tempfile
            tmp = tempfile.mkdtemp(prefix="sweep_dryrun_")
            _submit(C, tmp, runid, blocks, fast, cpus, walltime, True, True, False, target=rundir)
            print("  (the script writes into %s when submitted for real)" % rundir)
            return
        _write_run_env(rundir, resumed="SLURM", cpus=cpus, walltime=walltime)
        _submit(C, rundir, runid, blocks, fast, cpus, walltime, True, False, False)
        return
    workers = int(opts.get("--workers") or max(1, min(os.cpu_count() or 1, 128)))
    _write_run_env(rundir, resumed="local", workers=workers, host=socket.gethostname())
    code = _teed(os.path.join(rundir, "job.log"), _job_body, rundir, blocks, fast, workers, True)
    sys.exit(code)


def _progress(rundir):
    """Per-block counts of planned tasks, results and errors (reads plan_A.json / plan_B.json and tasks/)."""
    out = {}
    for pn in ("plan_A.json", "plan_B.json"):
        p = os.path.join(rundir, pn)
        if not os.path.exists(p):
            continue
        try:
            plan = json.load(open(p, encoding="utf-8"))
        except ValueError:
            continue
        for t in plan:
            b = out.setdefault(t["block"], [0, 0, 0])
            b[0] += 1
            d = os.path.join(rundir, "tasks", t["block"])
            if os.path.exists(os.path.join(d, t["id"] + ".json")):
                b[1] += 1
            elif os.path.exists(os.path.join(d, t["id"] + ".error.json")):
                b[2] += 1
    return out


def cmd_list(args):
    C = load_config()
    ds = run_dirs()
    if not ds:
        log("no run folders yet")
        return
    print("%-32s %-9s %-6s %s" % ("RUN (runs/<RUN>)", "JOB", "KIND", "STATE"))
    for d in ds:
        p = os.path.join(RUNS, d)
        job = " ".join(open(os.path.join(p, "JOBID")).read().split()) if os.path.exists(os.path.join(p, "JOBID")) else "-"
        if d.startswith("dryrun_"):
            kind, st = "dry", "job script only (nothing submitted)"
        elif d.startswith("probe_"):
            kind = "probe"
            txt = open(os.path.join(p, "probe.txt"), errors="replace").read() if os.path.exists(os.path.join(p, "probe.txt")) else ""
            m = re.findall(r"PROBE (OK|FAILED)", txt)
            st = "PROBE " + m[-1] if m else "(no verdict yet)"
        else:
            kind = "run"
            if os.path.exists(os.path.join(p, "SUMMARY.md")):
                ec = open(os.path.join(p, "EXIT_CODE")).read().strip() if os.path.exists(os.path.join(p, "EXIT_CODE")) else "?"
                st = "done: SUMMARY.md (exit %s)" % ec
            else:
                pr = _progress(p)
                st = ("tasks %d/%d" % (sum(v[1] for v in pr.values()), sum(v[0] for v in pr.values()))) if pr else "waiting or starting"
        if os.path.exists(os.path.join(C["RESULTS_DIR"], "%s%s.tar.gz" % (RESULTS_PREFIX, d))):
            st += " + tarball"
        print("%-32s %-9s %-6s %s" % (d, job, kind, st))


def cmd_status(args):
    C = load_config()
    if shutil.which("squeue"):
        rc, out = sh(["squeue", "-u", os.environ["USER"], "-o", "%.10i %.14j %.9T %.10M %.10l %R"], timeout=30)
        print(out.rstrip())
    rundir = resolve_run(args[0] if args else None)
    if not rundir:
        log("no run %s (all folders: python sweep.py list)" % (args[0] if args else "yet"))
        return
    name = os.path.basename(rundir)
    jobs = open(os.path.join(rundir, "JOBID")).read().split() if os.path.exists(os.path.join(rundir, "JOBID")) else []
    print("\nrun %s  (job %s)" % (name, " ".join(jobs) or "local"))
    if jobs and shutil.which("squeue"):
        qs = squeue_state(jobs[-1])
        if qs:
            est = start_estimate(jobs[-1]) if qs.startswith("PENDING") else ""
            print("  queue: %s%s" % (qs, ("   (SLURM start estimate %s)" % est) if est else ""))
    if os.path.exists(os.path.join(rundir, "probe.txt")):
        for l in open(os.path.join(rundir, "probe.txt"), errors="replace"):
            if re.search(r"PROBE|WARNING|FATAL|estimate", l):
                print("  " + l.rstrip())
        return
    pr = _progress(rundir)
    if not pr:
        print("  no task plan yet (the job has not started)")
    run_blocks = []
    try:
        run_blocks = json.load(open(os.path.join(rundir, "config.json"), encoding="utf-8")).get("blocks", [])
    except (OSError, ValueError):
        pass
    for b in ALL_BLOCKS:
        if b in pr:
            n, ok, bad = pr[b]
            print("  %-4s %4d/%-4d tasks done%s" % (b, ok, n, ("   %d FAILED" % bad) if bad else ""))
        elif b in run_blocks:
            print("  %-4s reduction only (phase C)" % b)
    if os.path.exists(os.path.join(rundir, "SUMMARY.md")):
        ec = open(os.path.join(rundir, "EXIT_CODE")).read().strip() if os.path.exists(os.path.join(rundir, "EXIT_CODE")) else "?"
        print("\nfinished (exit code %s); summary: python sweep.py show %s" % (ec, name))
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
    """MANIFEST.sha256 of the run folder (every file except itself and the live logs), then
    SOTHE-P2-SWEEP_results_<RUN>.tar.gz of the whole folder (sorted entries) and its .sha256."""
    rundir = os.path.realpath(rundir)
    base = os.path.basename(rundir)
    files = []
    for dp, dn, fn in os.walk(rundir):
        dn.sort()
        for n in sorted(fn):
            if n.startswith(".tmp_"):
                continue
            rel = os.path.relpath(os.path.join(dp, n), rundir).replace(os.sep, "/")
            if rel != "MANIFEST.sha256" and not LIVE_LOGS.match(rel):
                files.append(rel)
    with open(os.path.join(rundir, "MANIFEST.sha256"), "w", encoding="utf-8", newline="\n") as f:
        for rel in files:
            f.write("%s  ./%s\n" % (sha256(os.path.join(rundir, *rel.split("/"))), rel))
    name = "%s%s.tar.gz" % (RESULTS_PREFIX, base)
    out = os.path.join(dest, name)
    tmp = out + ".part"

    def norm(ti):
        ti.uid = ti.gid = 0
        ti.uname = ti.gname = ""
        return ti

    with tarfile.open(tmp, "w:gz") as tf:
        tf.add(rundir, arcname=base, recursive=False, filter=norm)
        for dp, dn, fn in os.walk(rundir):
            dn.sort()
            for d in dn:
                tf.add(os.path.join(dp, d), arcname=os.path.join(base, os.path.relpath(os.path.join(dp, d), rundir)).replace(os.sep, "/"),
                       recursive=False, filter=norm)
            for n in sorted(fn):
                if n.startswith(".tmp_"):
                    continue
                p = os.path.join(dp, n)
                tf.add(p, arcname=os.path.join(base, os.path.relpath(p, rundir)).replace(os.sep, "/"), recursive=False, filter=norm)
    os.replace(tmp, out)
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
    print("SOTHE-P2-SWEEP %s (Python %s at %s)" % (version(), sys.version.split()[0], sys.executable))


def cmd_blocks(args):
    C = load_config()
    default = parse_blocks(C["BLOCKS"])
    print("%-4s %-8s %-18s %s" % ("", "DEFAULT", "FOLDER", "WHAT"))
    for b, what, folder in BLOCK_INFO:
        print("%-4s %-8s %-18s %s" % (b, "yes" if b in default else "no", folder, what))
    print("\ndependencies: B2 <- B1, B4;  B4 <- B3;  B10 <- B1, B2, B3, B4, B5, B7 (added automatically)")
    print("subset:  python sweep.py submit --blocks \"B6 B7\"     (or run --blocks ...)")


COMMANDS = {"install": cmd_install, "probe": cmd_probe, "submit": cmd_submit, "run": cmd_run, "resume": cmd_resume, "list": cmd_list,
            "ls": cmd_list, "status": cmd_status, "show": cmd_show, "repack": cmd_repack, "verify": cmd_verify, "blocks": cmd_blocks,
            "version": cmd_version, "--version": cmd_version, "help": cmd_help, "-h": cmd_help, "--help": cmd_help,
            "_job": cmd_job, "_probe": cmd_probe_job}
NEEDS_PACK = ("install", "probe", "submit", "run", "resume", "verify", "_job", "_probe")


def main(argv=None):
    argv = list(sys.argv[1:] if argv is None else argv)
    cmd = argv[0] if argv else "help"
    if cmd not in COMMANDS:
        die("unknown command %s (see: python sweep.py help)" % cmd)
    if cmd in NEEDS_PACK and not os.path.exists(os.path.join(ROOT, "sothe_p2", "kit.py")):
        die("the pack folder is not at %s: run the commands from the unpacked SOTHE-P2-SWEEP folder (python sweep.py ...)" % ROOT)
    try:
        COMMANDS[cmd](argv[1:])
    except BrokenPipeError:
        try:
            sys.stdout = open(os.devnull, "w")
        except OSError:
            pass
