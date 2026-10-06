#!/usr/bin/env python3
"""reproduce.py -- rebuild the figures and the printed numbers of the article from the archived runs and compare them with output/.

Usage
    python reproduce.py --p2-results SOTHE-P2_results_20260927T175224Z.tar.gz \
                        --sweep-results SOTHE-P2-SWEEP_results_20260928T131100Z.tar.gz \
                        [--sweep-pkg SOTHE-P2-SWEEP | --sweep-pkg-tar SOTHE-P2-SWEEP_v1.0.0.tar.gz] [--workdir DIR]
                        [--zero-sector [--grid]]

The two results archives are files of the Zenodo record that carries this folder. The SOTHE-P2-SWEEP package is taken from
../SOTHE-P2-SWEEP (the layout of the GitHub branch), from --sweep-pkg (an unpacked folder) or from --sweep-pkg-tar.
Nothing in the archives or in output/ is modified: everything is unpacked and rebuilt in a work folder (default: a temporary
one, removed at the end).

Requires Python 3.9 or newer with NumPy, SciPy and Matplotlib. The comparison has three levels: bytes (figure files, when
the Matplotlib build is the one that wrote output/), arrays (figure_data.npz: exactly equal, or equal to a relative 1e-12)
and numbers (article_numbers.json: relative 1e-12). With --zero-sector the zero-sector tests of the supplement are run as well
(zero_sector_tests.py, about 2 minutes; --grid adds the grid scan, about 4 minutes) and every statement they check must hold.
Exit status 0 when the arrays and the numbers agree (and, if requested, the zero-sector statements hold)."""
import argparse
import hashlib
import json
import os
import shutil
import subprocess
import sys
import tarfile
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
RUNID_P2 = "20260927T175224Z"
RUNID_SWEEP = "20260928T131100Z"
EXPECTED = {   # sha256 of the archives of the Zenodo record (a mismatch is reported, not fatal)
    "SOTHE-P2_results_%s.tar.gz" % RUNID_P2: "cdb4a29f199251c7a749dadf118c5be1b77d3e946844badfece14c2dc59e6786",
    "SOTHE-P2-SWEEP_results_%s.tar.gz" % RUNID_SWEEP: "af959da19a8c48bb88369b65d1d9ec462fbf00823db987fe0de5e979fcc7170f",
    "SOTHE-P2-SWEEP_v1.0.0.tar.gz": "4a79e4e9c28a29210264b79b0246876941fafd2a2314dcb66755a31b6d042b3e",
}


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def unpack(tar_path, dest):
    os.makedirs(dest, exist_ok=True)
    with tarfile.open(tar_path, "r:*") as tf:
        for m in tf.getmembers():
            parts = [p for p in m.name.split("/") if p not in ("", ".")]
            if m.name.startswith("/") or ".." in parts or m.issym() or m.islnk():
                sys.exit("unsafe entry in %s: %s" % (tar_path, m.name))
        try:
            tf.extractall(dest, filter="data")     # Python >= 3.12
        except TypeError:
            tf.extractall(dest)


def check_archive(path):
    name = os.path.basename(path)
    if name in EXPECTED:
        ok = sha256(path) == EXPECTED[name]
        print("  %-52s sha256 %s" % (name, "as recorded" if ok else "DIFFERS from the recorded value"))
        return ok
    return True


def same_json(a, b, rtol=1e-12, path=""):
    if isinstance(a, dict) and isinstance(b, dict):
        if sorted(a) != sorted(b):
            return ["keys differ at %s" % (path or "/")]
        out = []
        for k in a:
            out += same_json(a[k], b[k], rtol, path + "/" + k)
        return out
    if isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            return ["length differs at %s" % path]
        out = []
        for i, (u, v) in enumerate(zip(a, b)):
            out += same_json(u, v, rtol, "%s[%d]" % (path, i))
        return out
    if isinstance(a, bool) or isinstance(b, bool) or isinstance(a, str) or isinstance(b, str) or a is None or b is None:
        return [] if a == b else ["%s: %r != %r" % (path, a, b)]
    x, y = float(a), float(b)
    return [] if abs(x - y) <= rtol * max(abs(x), abs(y), 1e-300) else ["%s: %r != %r" % (path, a, b)]


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--p2-results", required=True, help="SOTHE-P2_results_%s.tar.gz" % RUNID_P2)
    ap.add_argument("--sweep-results", required=True, help="SOTHE-P2-SWEEP_results_%s.tar.gz" % RUNID_SWEEP)
    ap.add_argument("--sweep-pkg", help="unpacked SOTHE-P2-SWEEP folder (default: ../SOTHE-P2-SWEEP)")
    ap.add_argument("--sweep-pkg-tar", help="SOTHE-P2-SWEEP_v1.0.0.tar.gz")
    ap.add_argument("--workdir", help="work folder (default: temporary, removed at the end)")
    ap.add_argument("--zero-sector", action="store_true", help="also run zero_sector_tests.py (about 2 minutes)")
    ap.add_argument("--grid", action="store_true", help="with --zero-sector: also its grid scan (about 4 minutes)")
    a = ap.parse_args()

    tmp = a.workdir or tempfile.mkdtemp(prefix="sothe_article_")
    os.makedirs(tmp, exist_ok=True)
    try:
        print("archives")
        allok = check_archive(a.p2_results) & check_archive(a.sweep_results)
        pkg = a.sweep_pkg
        if a.sweep_pkg_tar:
            allok &= check_archive(a.sweep_pkg_tar)
            unpack(a.sweep_pkg_tar, os.path.join(tmp, "pkg"))
            pkg = os.path.join(tmp, "pkg", "SOTHE-P2-SWEEP")
        if not pkg:
            pkg = os.path.join(HERE, "..", "SOTHE-P2-SWEEP")
        if not os.path.isfile(os.path.join(pkg, "sweep.py")):
            sys.exit("SOTHE-P2-SWEEP package not found at %s (use --sweep-pkg or --sweep-pkg-tar)" % pkg)
        unpack(a.p2_results, os.path.join(tmp, "p2"))
        unpack(a.sweep_results, os.path.join(tmp, "sweep"))
        run_p2 = os.path.join(tmp, "p2", RUNID_P2)
        run_sw = os.path.join(tmp, "sweep", RUNID_SWEEP)
        figdata = os.path.join(run_p2, "S5_figures", "figure_data.npz")
        for p in (run_p2, run_sw, figdata):
            if not os.path.exists(p):
                sys.exit("missing after unpacking: %s" % p)
        out = os.path.join(tmp, "out")
        os.makedirs(out, exist_ok=True)
        env = dict(os.environ, PYTHONDONTWRITEBYTECODE="1", SOTHE_SWEEP_PKG=os.path.abspath(pkg), SOTHE_SWEEP_RUN=run_sw,
                   SOTHE_P2_FIGDATA=figdata, SOTHE_ARTICLE_OUT=out, SOTHE_ARTICLE_NUMBERS=os.path.join(out, "article_numbers.json"))
        for script in ("make_article_figures.py", "make_article_numbers.py"):
            print("running %s" % script)
            p = subprocess.run([sys.executable, os.path.join(HERE, script)], env=env, cwd=tmp,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            if p.returncode != 0:
                sys.stdout.write(p.stdout.decode("utf-8", "replace")[-3000:])
                sys.exit("%s failed" % script)

        import numpy as np
        ref = os.path.join(HERE, "output")
        print("\ncomparison with output/")
        n_bytes = n_files = 0
        for name in sorted(os.listdir(ref)):
            if name in ("figure_data.npz", "article_numbers.json", "zero_sector_tests.json"):
                continue
            n_files += 1
            same = os.path.exists(os.path.join(out, name)) and sha256(os.path.join(out, name)) == sha256(os.path.join(ref, name))
            n_bytes += int(same)
            if not same:
                print("  %-30s bytes differ (another Matplotlib build?)" % name)
        print("  figure files: %d of %d byte-identical" % (n_bytes, n_files))
        za, zb = np.load(os.path.join(ref, "figure_data.npz")), np.load(os.path.join(out, "figure_data.npz"))
        bad = []
        if sorted(za.files) != sorted(zb.files):
            bad.append("array names differ")
        for k in za.files:
            if k in zb.files:
                x, y = za[k], zb[k]
                if x.shape != y.shape or not (np.array_equal(x, y) or np.allclose(x, y, rtol=1e-12, atol=0)):
                    bad.append("array %s differs" % k)
        ident = all(np.array_equal(za[k], zb[k]) for k in za.files if k in zb.files) and not bad
        print("  figure_data.npz: %d arrays, %s" % (len(za.files), "identical" if ident else ("equal to 1e-12" if not bad else "DIFFERENT")))
        ja = json.load(open(os.path.join(ref, "article_numbers.json"), encoding="utf-8"))
        jb = json.load(open(os.path.join(out, "article_numbers.json"), encoding="utf-8"))
        jd = same_json(ja, jb)
        bytes_json = sha256(os.path.join(ref, "article_numbers.json")) == sha256(os.path.join(out, "article_numbers.json"))
        print("  article_numbers.json: %d entries, %s" % (len(ja), "byte-identical" if bytes_json else ("equal to 1e-12" if not jd else "DIFFERENT")))
        for b in bad + jd[:10]:
            print("   !", b)
        ok = not bad and not jd
        if a.zero_sector:
            print("\nzero-sector tests (zero_sector_tests.py%s)" % (" --grid" if a.grid else ""))
            zenv = dict(env, SOTHE_P2_RUN=run_p2)
            p = subprocess.run([sys.executable, os.path.join(HERE, "zero_sector_tests.py")] + (["--grid"] if a.grid else []), env=zenv, cwd=tmp,
                               stdout=subprocess.PIPE, stderr=subprocess.STDOUT)
            txt = p.stdout.decode("utf-8", "replace")
            for ln in txt.splitlines():
                if ln.startswith(("  [PASS", "  [FAIL")) or ln.startswith(("RESULT:", "NOTE:")):      # "[PASS, within the rounding ...]" lines and the NOTE line included
                    print("  " + ln.strip() if ln.startswith(("RESULT:", "NOTE:")) else ln)
            zs_ok = p.returncode == 0
            if p.returncode != 0 and "RESULT:" not in txt:
                sys.stdout.write(txt[-3000:])
            zref = os.path.join(ref, "zero_sector_tests.json")
            zout = os.path.join(out, "zero_sector_tests.json")
            if zs_ok and os.path.exists(zref) and os.path.exists(zout):
                ra = {x["id"]: x["holds"] for x in json.load(open(zref, encoding="utf-8"))["statements"]}
                rb = {x["id"]: x["holds"] for x in json.load(open(zout, encoding="utf-8"))["statements"]}
                common = sorted(set(ra) & set(rb))
                print("  statements also in output/zero_sector_tests.json: %d, same verdict: %s" % (len(common), "yes" if all(ra[k] == rb[k] for k in common) else "NO"))
            ok = ok and zs_ok
        print("\nRESULT: %s" % ("PASS" if ok else "FAIL"))
        return 0 if ok else 1
    finally:
        if not a.workdir:
            shutil.rmtree(tmp, ignore_errors=True)


if __name__ == "__main__":
    sys.exit(main())
