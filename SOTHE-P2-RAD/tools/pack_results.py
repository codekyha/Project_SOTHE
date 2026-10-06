#!/usr/bin/env python3
"""Pack a results folder of SOTHE-P2-RAD into SOTHE-P2-RAD_results_<stamp>.tar.gz (+ .sha256), reproducibly.

  python tools/pack_results.py RESULTS/<stamp> [--out DIR]

The folder must carry MANIFEST.sha256 (python rad.py manifest --results RESULTS/<stamp>), which is checked first.  Entries are
sorted, owner root/root, permissions 0644 (directories 0755), every mtime = the end of the run (finished_utc of RUN.json),
and the gzip header has no name or time.  The archive unpacks into one folder named <stamp>.  Standard library only."""
import argparse
import calendar
import gzip
import hashlib
import io
import json
import os
import sys
import tarfile
import time


def sha256(path):
    h = hashlib.sha256()
    with open(path, "rb") as f:
        for b in iter(lambda: f.read(1 << 20), b""):
            h.update(b)
    return h.hexdigest()


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("folder")
    ap.add_argument("--out", default=None)
    a = ap.parse_args(argv)
    root = os.path.abspath(a.folder)
    stamp = os.path.basename(root.rstrip(os.sep))
    man = os.path.join(root, "MANIFEST.sha256")
    if not os.path.isfile(man):
        sys.exit("%s has no MANIFEST.sha256 (python rad.py manifest --results %s)" % (root, a.folder))
    listed = []
    with open(man, encoding="utf-8") as f:
        for ln in f:
            if ln.strip():
                d, rel = ln.rstrip("\n").split("  ", 1)
                p = os.path.join(root, *rel.split("/"))
                if not os.path.isfile(p) or sha256(p) != d:
                    sys.exit("MANIFEST.sha256: %s is missing or changed" % rel)
                listed.append(rel)
    files = []
    for dp, dn, fn in os.walk(root):
        dn.sort()
        for n in fn:
            files.append(os.path.relpath(os.path.join(dp, n), root).replace(os.sep, "/"))
    extra = sorted(set(files) - set(listed) - {"MANIFEST.sha256"})
    if extra:
        sys.exit("files not listed in MANIFEST.sha256: %s" % ", ".join(extra[:5]))
    with open(os.path.join(root, "RUN.json"), encoding="utf-8") as f:
        fin = json.load(f)["finished_utc"]
    mtime = calendar.timegm(time.strptime(fin[:19], "%Y-%m-%dT%H:%M:%S"))
    files = sorted(files)
    dirs = sorted({"/".join(r.split("/")[:i]) for r in files for i in range(1, len(r.split("/")))})
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w", format=tarfile.PAX_FORMAT) as tf:
        def info(name, isdir, size=0):
            ti = tarfile.TarInfo(name)
            ti.type = tarfile.DIRTYPE if isdir else tarfile.REGTYPE
            ti.mode = 0o755 if isdir else 0o644
            ti.size = size
            ti.mtime = mtime
            ti.uid = ti.gid = 0
            ti.uname = ti.gname = "root"
            return ti
        tf.addfile(info(stamp, True))
        for rel, isdir in sorted([(d, True) for d in dirs] + [(f, False) for f in files]):
            if isdir:
                tf.addfile(info(stamp + "/" + rel, True))
                continue
            with open(os.path.join(root, *rel.split("/")), "rb") as f:
                data = f.read()
            tf.addfile(info(stamp + "/" + rel, False, len(data)), io.BytesIO(data))
    outdir = os.path.abspath(a.out) if a.out else os.path.dirname(root)
    os.makedirs(outdir, exist_ok=True)
    name = "SOTHE-P2-RAD_results_%s.tar.gz" % stamp
    out = os.path.join(outdir, name)
    with open(out, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as gz:
            gz.write(buf.getvalue())
    with open(out + ".sha256", "w", encoding="utf-8", newline="\n") as f:
        f.write("%s  %s\n" % (sha256(out), name))
    print("%s: %d files, %.2f MB, sha256 %s" % (out, len(files), os.path.getsize(out) / 1e6, sha256(out)))
    return 0


if __name__ == "__main__":
    sys.exit(main())
