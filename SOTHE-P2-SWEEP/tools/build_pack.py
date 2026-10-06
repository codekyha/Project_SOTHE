#!/usr/bin/env python3
"""Build the SOTHE-P2-SWEEP pack tarball and its MANIFEST.   Standard library only.

  python tools/build_pack.py                  vendored-code check, MANIFEST.sha256, dist/SOTHE-P2-SWEEP_v<VERSION>.tar.gz + .sha256
  python tools/build_pack.py --check          verify MANIFEST.sha256 against the files (and the vendored sothe_p2/), build nothing
  python tools/build_pack.py --bundle         also write dist/SOURCE_BUNDLE_SOTHE-P2-SWEEP_v<VERSION>.txt (every text file, with sha256)
  python tools/build_pack.py --out DIR        write into DIR (default: dist/ next to the pack)

The tarball is reproducible: entries sorted, owner root/root, permissions 0644 (0755 for sweep.py, tools/*.py and
directories), every mtime = SOURCE_DATE_EPOCH if set, else RELEASE_DATE below (00:00 UTC), and a gzip header without
name or time.  It unpacks into one folder, SOTHE-P2-SWEEP/.  Excluded: runs/, dist/, validation/, .git/, __pycache__/,
config/site.env*, results tarballs, editor and OS files.
The vendored package sothe_p2/ must be byte-identical to SOTHE-P2 1.0.0 (reference/SOTHE-P2_1.0.0_MANIFEST.sha256)."""
import argparse
import calendar
import fnmatch
import gzip
import hashlib
import io
import os
import re
import sys
import tarfile

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
TOP = "SOTHE-P2-SWEEP"
RELEASE_DATE = (2026, 9, 28)
EXCL_DIRS = {"runs", "dist", "validation", ".git", "__pycache__", ".pytest_cache", ".mypy_cache", ".idea", ".vscode", "build"}
EXCL_FILES = ["MANIFEST.sha256", "config/site.env", "config/site.env.detected", "*.pyc", "*.pyo", "SOTHE-P2-SWEEP_results_*",
              "SOTHE-P2-SWEEP_v*.tar.gz", "SOTHE-P2-SWEEP_v*.tar.gz.sha256", "*.part", ".tmp_*", ".DS_Store", "Thumbs.db", "desktop.ini",
              "*~", "*.swp"]
EXEC = ["sweep.py", "tools/*.py", "tests/*.py"]
TEXT = (".py", ".md", ".txt", ".env", ".example", ".cff", ".toml", ".csv", ".json", ".sha256", "LICENSE", "VERSION", ".gitignore")


def version():
    with open(os.path.join(ROOT, "VERSION"), encoding="utf-8") as f:
        return f.read().strip()


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def pack_files():
    out = []
    for dp, dn, fn in os.walk(ROOT):
        dn[:] = sorted(d for d in dn if d not in EXCL_DIRS)
        for n in sorted(fn):
            rel = os.path.relpath(os.path.join(dp, n), ROOT).replace(os.sep, "/")
            if any(fnmatch.fnmatch(rel, p) or fnmatch.fnmatch(n, p) for p in EXCL_FILES):
                continue
            out.append(rel)
    return sorted(out)


def check_vendored():
    """sothe_p2/ against the SOTHE-P2 1.0.0 MANIFEST: every listed file identical, no file added."""
    ref = os.path.join(ROOT, "reference", "SOTHE-P2_1.0.0_MANIFEST.sha256")
    listed, bad = set(), []
    with open(ref, encoding="utf-8") as f:
        for ln in f:
            if not ln.strip():
                continue
            d, n = ln.split(None, 1)
            n = n.strip()
            n = n[2:] if n.startswith("./") else n
            if not n.startswith("sothe_p2/"):
                continue
            listed.add(n)
            p = os.path.join(ROOT, *n.split("/"))
            if not os.path.exists(p) or sha256(p) != d:
                bad.append(n)
    extra = [r for r in pack_files() if r.startswith("sothe_p2/") and r not in listed]
    print("vendored sothe_p2/: %d files identical to SOTHE-P2 1.0.0, %d differ, %d added" % (len(listed) - len(bad), len(bad), len(extra)))
    for b in bad + extra:
        print("   " + b)
    return not bad and not extra


def check_versions():
    ver = version()
    with open(os.path.join(ROOT, "sweep_p2", "__init__.py"), encoding="utf-8") as f:
        m = re.search(r'(?m)^__version__\s*=\s*"([^"]+)"', f.read())
    if not m or m.group(1) != ver:
        sys.exit("version mismatch: VERSION %s, sweep_p2/__init__.py %s" % (ver, m.group(1) if m else None))
    return ver


def write_manifest(files):
    with open(os.path.join(ROOT, "MANIFEST.sha256"), "w", encoding="utf-8", newline="\n") as f:
        for rel in files:
            f.write("%s  ./%s\n" % (sha256(os.path.join(ROOT, *rel.split("/"))), rel))


def check_manifest():
    man = os.path.join(ROOT, "MANIFEST.sha256")
    if not os.path.exists(man):
        print("no MANIFEST.sha256")
        return False
    listed, bad = set(), []
    with open(man, encoding="utf-8") as f:
        for ln in f:
            if not ln.strip():
                continue
            digest, name = ln.rstrip("\n").split(None, 1)
            n = name.strip().lstrip("*")
            rel = n[2:] if n.startswith("./") else n
            listed.add(rel)
            p = os.path.join(ROOT, *rel.split("/"))
            if not os.path.exists(p):
                bad.append("missing  " + rel)
            elif sha256(p) != digest:
                bad.append("changed  " + rel)
    extra = [r for r in pack_files() if r not in listed]
    for b in bad:
        print(b)
    for e in extra:
        print("unlisted " + e)
    print("MANIFEST: %d files, %d changed or missing, %d not listed" % (len(listed), len(bad), len(extra)))
    return not bad and not extra


def release_epoch():
    if os.environ.get("SOURCE_DATE_EPOCH"):
        return int(os.environ["SOURCE_DATE_EPOCH"])
    return calendar.timegm(RELEASE_DATE + (0, 0, 0))


def build(outdir):
    ver = check_versions()
    if not check_vendored():
        sys.exit("the vendored sothe_p2/ differs from SOTHE-P2 1.0.0: restore it before building")
    files = pack_files()
    write_manifest(files)
    files = sorted(files + ["MANIFEST.sha256"])
    mtime = release_epoch()
    dirs = set()
    for rel in files:
        parts = rel.split("/")[:-1]
        for i in range(1, len(parts) + 1):
            dirs.add("/".join(parts[:i]))
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w", format=tarfile.PAX_FORMAT) as tf:
        def info(name, isdir, size=0, mode=0o644):
            ti = tarfile.TarInfo(name)
            ti.type = tarfile.DIRTYPE if isdir else tarfile.REGTYPE
            ti.mode = 0o755 if isdir else mode
            ti.size = size
            ti.mtime = mtime
            ti.uid = ti.gid = 0
            ti.uname = ti.gname = "root"
            return ti
        tf.addfile(info(TOP, True))
        for rel, isdir in sorted([(d, True) for d in dirs] + [(f, False) for f in files]):
            if isdir:
                tf.addfile(info(TOP + "/" + rel, True))
                continue
            with open(os.path.join(ROOT, *rel.split("/")), "rb") as f:
                data = f.read()
            mode = 0o755 if any(fnmatch.fnmatch(rel, e) for e in EXEC) else 0o644
            tf.addfile(info(TOP + "/" + rel, False, len(data), mode), io.BytesIO(data))
    os.makedirs(outdir, exist_ok=True)
    name = "SOTHE-P2-SWEEP_v%s.tar.gz" % ver
    out = os.path.join(outdir, name)
    with open(out, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as gz:
            gz.write(buf.getvalue())
    with open(out + ".sha256", "w", encoding="utf-8", newline="\n") as f:
        f.write("%s  %s\n" % (sha256(out), name))
    print("%s: %d files, %.2f MB, sha256 %s" % (out, len(files), os.path.getsize(out) / 1e6, sha256(out)))
    return out, files


def bundle(outdir, files):
    """All text files of the pack in one file (for reading without unpacking), each with its sha256 and size."""
    ver = version()
    out = os.path.join(outdir, "SOURCE_BUNDLE_SOTHE-P2-SWEEP_v%s.txt" % ver)
    bar = "=" * 100
    with open(out, "w", encoding="utf-8", newline="\n") as f:
        f.write("SOURCE BUNDLE of SOTHE-P2-SWEEP %s: every text file of the pack, verbatim (binary files are listed only).\n" % ver)
        f.write("Pack tarball: SOTHE-P2-SWEEP_v%s.tar.gz (MANIFEST.sha256 inside).\n\n" % ver)
        for rel in files:
            p = os.path.join(ROOT, *rel.split("/"))
            is_text = rel.endswith(TEXT) or os.path.basename(rel) in TEXT
            f.write("%s\nFILE: %s   sha256 %s   %d bytes%s\n%s\n" % (bar, rel, sha256(p), os.path.getsize(p), "" if is_text else "   (binary, not shown)", bar))
            if is_text:
                with open(p, encoding="utf-8", errors="replace") as g:
                    f.write(g.read().rstrip("\n") + "\n\n")
    print("%s (%.2f MB)" % (out, os.path.getsize(out) / 1e6))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=os.path.join(ROOT, "dist"))
    ap.add_argument("--check", action="store_true")
    ap.add_argument("--bundle", action="store_true")
    a = ap.parse_args(argv)
    if a.check:
        sys.exit(0 if (check_vendored() and check_manifest()) else 1)
    out, files = build(a.out)
    if a.bundle:
        bundle(a.out, files)


if __name__ == "__main__":
    main()
