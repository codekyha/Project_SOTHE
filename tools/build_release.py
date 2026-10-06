#!/usr/bin/env python3
"""Build the SOTHE-P2 release tarball (the Altay/GitHub/Zenodo pack) and its MANIFEST.   Standard library only.

  python tools/build_release.py                 MANIFEST.sha256 of the pack, then dist/SOTHE-P2_v<VERSION>.tar.gz + .sha256
  python tools/build_release.py --check         verify MANIFEST.sha256 against the files, build nothing
  python tools/build_release.py --doi 10.5281/zenodo.1234567
                                                write the reserved Zenodo version DOI into CITATION.cff and README.md, then build
  python tools/build_release.py --out DIR       write the tarball into DIR (default: dist/ next to the pack)

The tarball is reproducible: entries sorted, owner root/root, permissions 0644 (0755 for p2.py, tools/*.py and
directories), every mtime = SOURCE_DATE_EPOCH if set, else the date-released of CITATION.cff (00:00 UTC), and a
gzip header without name or time.  It unpacks into one folder, SOTHE-P2/.  Excluded: runs/, dist/, .git/,
__pycache__/, config/site.env*, results tarballs, editor and OS files."""
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
TOP = "SOTHE-P2"
EXCL_DIRS = {"runs", "dist", ".git", "__pycache__", ".pytest_cache", ".mypy_cache", ".idea", ".vscode", ".ipynb_checkpoints",
             "build", "sothe_p2.egg-info"}
EXCL_FILES = ["MANIFEST.sha256", "config/site.env", "config/site.env.detected", "*.pyc", "*.pyo", "SOTHE-P2_results_*",
              "SOTHE-P2_v*.tar.gz", "SOTHE-P2_v*.tar.gz.sha256", ".DS_Store", "Thumbs.db", "desktop.ini", "*~", "*.swp"]
EXEC = ["p2.py", "tools/*.py"]
DOI_RE = r"10\.5281/zenodo\.\d+"


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
    with open(os.path.join(ROOT, "CITATION.cff"), encoding="utf-8") as f:
        m = re.search(r'(?m)^date-released:\s*"?(\d{4})-(\d{2})-(\d{2})', f.read())
    if not m:
        sys.exit("CITATION.cff has no date-released; set SOURCE_DATE_EPOCH")
    return calendar.timegm((int(m.group(1)), int(m.group(2)), int(m.group(3)), 0, 0, 0))


def set_doi(doi):
    """Put the reserved version DOI into CITATION.cff (doi: line) and README.md (the DOI line of the citation block)."""
    if not re.fullmatch(DOI_RE, doi):
        sys.exit("--doi must look like 10.5281/zenodo.1234567")
    p = os.path.join(ROOT, "CITATION.cff")
    with open(p, encoding="utf-8") as f:
        s = f.read()
    s2, n = re.subn(r'(?m)^# doi: "10\.5281/zenodo\.RESERVED".*$|^doi: "%s"$' % DOI_RE, 'doi: "%s"' % doi, s)
    if n != 1:
        sys.exit("CITATION.cff: expected one doi line (or its placeholder), found %d" % n)
    with open(p, "w", encoding="utf-8", newline="\n") as f:
        f.write(s2)
    p = os.path.join(ROOT, "README.md")
    with open(p, encoding="utf-8") as f:
        s = f.read()
    s2, n = re.subn(r"(?m)^(Zenodo DOI of this version: ).*$", r"\g<1>[%s](https://doi.org/%s)" % (doi, doi), s)
    if n != 1:
        sys.exit("README.md: expected one 'Zenodo DOI of this version:' line, found %d" % n)
    with open(p, "w", encoding="utf-8", newline="\n") as f:
        f.write(s2)
    print("DOI %s written into CITATION.cff and README.md" % doi)


def check_versions():
    """VERSION, sothe_p2/__init__.py, pyproject.toml and CITATION.cff must carry the same version."""
    ver = version()
    found = {"VERSION": ver}
    for rel, pat in (("sothe_p2/__init__.py", r'(?m)^__version__\s*=\s*"([^"]+)"'), ("pyproject.toml", r'(?m)^version\s*=\s*"([^"]+)"'),
                     ("CITATION.cff", r'(?m)^version:\s*"?([^"\s]+)"?')):
        with open(os.path.join(ROOT, *rel.split("/")), encoding="utf-8") as f:
            m = re.search(pat, f.read())
        found[rel] = m.group(1) if m else None
    bad = {k: v for k, v in found.items() if v != ver}
    if bad:
        sys.exit("version mismatch (VERSION = %s): %s" % (ver, ", ".join("%s = %s" % kv for kv in bad.items())))
    return ver


def build(outdir):
    ver = check_versions()
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
        entries = sorted([(d, True) for d in dirs] + [(f, False) for f in files])
        for rel, isdir in entries:
            if isdir:
                tf.addfile(info(TOP + "/" + rel, True))
                continue
            p = os.path.join(ROOT, *rel.split("/"))
            mode = 0o755 if any(fnmatch.fnmatch(rel, e) for e in EXEC) else 0o644
            with open(p, "rb") as f:
                data = f.read()
            tf.addfile(info(TOP + "/" + rel, False, len(data), mode), io.BytesIO(data))
    os.makedirs(outdir, exist_ok=True)
    name = "SOTHE-P2_v%s.tar.gz" % ver
    out = os.path.join(outdir, name)
    with open(out, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as gz:
            gz.write(buf.getvalue())
    with open(out + ".sha256", "w", encoding="utf-8", newline="\n") as f:
        f.write("%s  %s\n" % (sha256(out), name))
    print("%s: %d files, %.2f MB, sha256 %s" % (out, len(files), os.path.getsize(out) / 1e6, sha256(out)))
    return out


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--out", default=os.path.join(ROOT, "dist"))
    ap.add_argument("--doi", default="")
    ap.add_argument("--check", action="store_true")
    a = ap.parse_args(argv)
    if a.check:
        sys.exit(0 if check_manifest() else 1)
    if a.doi:
        set_doi(a.doi)
    return build(a.out)


if __name__ == "__main__":
    main()
