#!/usr/bin/env python3
"""Assemble the files for the manual Zenodo upload of a SOTHE-P2 release (own record, own DOI).   Standard library only.

  python tools/make_zenodo_pack.py
        builds (or rebuilds) dist/SOTHE-P2_v<VERSION>.tar.gz with tools/build_release.py, then writes
        dist/zenodo_SOTHE-P2_v<VERSION>/ with:
          SOTHE-P2_v<VERSION>.tar.gz (+ .sha256)   the code, reference outputs and archived trajectory
          README_ZENODO.md                         the record description (paste it, or upload it as a file)
          zenodo_metadata.json                     every field of the upload form (also valid for the Zenodo REST API)
          CHECKSUMS.sha256                         sha256 of every file above and below
          UPLOAD_STEPS.txt                         the manual steps
  options
    --results PATH   add a run: a SOTHE-P2_results_<RUN>.tar.gz (its .sha256 is checked if present) or a run folder
                     (packed reproducibly as SOTHE-P2_results_<RUN>.tar.gz); repeatable
    --extra PATH     add any other file (e.g. an archive of canonical MATLAB outputs); repeatable
    --doi DOI        the version DOI reserved on Zenodo (also written into CITATION.cff and README.md before the build)
    --date DATE      publication date for the metadata, YYYY-MM-DD (default: today, UTC)
    --out DIR        default dist/ next to the pack
Nothing is uploaded: the upload and the publication are manual steps of the PI (RELEASE.md)."""
import argparse
import datetime
import gzip
import io
import json
import os
import shutil
import sys
import tarfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import build_release as B  # noqa: E402

REPO = "https://github.com/codekyha/Project_SOTHE"
REF25_DOI = "10.5281/zenodo.20713660"
REF25_PAPER_DOI = "10.1088/1361-6382/ae811b"
CREATOR = {"name": "Oguz, Hasan", "affiliation": "Istanbul Okan University; Pamukkale University", "orcid": "0000-0001-7484-4415"}
KEYWORDS = ["analog gravity", "optical event horizon", "soliton", "Bogoliubov-de Gennes", "resonant radiation",
            "third-order dispersion", "generalized nonlinear Schroedinger equation", "surface gravity", "Python"]


def sha256(p):
    return B.sha256(p)


def pack_run_folder(rundir, outdir):
    """A run folder -> SOTHE-P2_results_<RUN>.tar.gz, reproducible (sorted, fixed owner, mtime of the release)."""
    rundir = os.path.realpath(rundir)
    base = os.path.basename(rundir)
    name = "SOTHE-P2_results_%s.tar.gz" % base
    out = os.path.join(outdir, name)
    mtime = B.release_epoch()
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w", format=tarfile.PAX_FORMAT) as tf:
        for dp, dn, fn in os.walk(rundir):
            dn.sort()
            rel_d = os.path.relpath(dp, rundir).replace(os.sep, "/")
            arc_d = base if rel_d == "." else base + "/" + rel_d
            ti = tarfile.TarInfo(arc_d)
            ti.type, ti.mode, ti.mtime, ti.uid, ti.gid, ti.uname, ti.gname = tarfile.DIRTYPE, 0o755, mtime, 0, 0, "root", "root"
            tf.addfile(ti)
            for n in sorted(fn):
                p = os.path.join(dp, n)
                if os.path.islink(p):
                    continue
                with open(p, "rb") as f:
                    data = f.read()
                ti = tarfile.TarInfo(arc_d + "/" + n)
                ti.size, ti.mode, ti.mtime, ti.uid, ti.gid, ti.uname, ti.gname = len(data), 0o644, mtime, 0, 0, "root", "root"
                tf.addfile(ti, io.BytesIO(data))
    with open(out, "wb") as raw:
        with gzip.GzipFile(filename="", mode="wb", fileobj=raw, mtime=0, compresslevel=9) as gz:
            gz.write(buf.getvalue())
    with open(out + ".sha256", "w", encoding="utf-8", newline="\n") as f:
        f.write("%s  %s\n" % (sha256(out), name))
    return out


def add_results(path, outdir):
    if os.path.isdir(path):
        return pack_run_folder(path, outdir)
    if not os.path.isfile(path):
        sys.exit("--results: no such file or folder: %s" % path)
    dst = os.path.join(outdir, os.path.basename(path))
    shutil.copyfile(path, dst)
    side = path + ".sha256"
    if os.path.exists(side):
        want = open(side, encoding="utf-8").read().split()[0]
        if want != sha256(dst):
            sys.exit("--results: %s does not match its .sha256 (re-transfer it in binary mode)" % path)
        shutil.copyfile(side, dst + ".sha256")
    return dst


def description(ver, results, extras):
    items = ["<li><code>SOTHE-P2_v%s.tar.gz</code>: the Python package (source, tests, documentation), the archived entropy "
             "trajectory <code>data/entropy_trajectory.csv</code>, and reference outputs of the MATLAB R2025b and Octave 8.4.0 "
             "implementations (<code>reference/</code>).</li>" % ver]
    for r in results:
        items.append("<li><code>%s</code>: the output folder of one complete run (reports, data products, figures, "
                     "<code>SUMMARY.md</code>, <code>MANIFEST.sha256</code>).</li>" % os.path.basename(r))
    for e in extras:
        items.append("<li><code>%s</code></li>" % os.path.basename(e))
    return ("<p>SOTHE-P2 %s computes the numerical results of the second paper of the SOTHE project (solitonic event horizons "
            "in dispersive media): Bogoliubov-de Gennes spectra about the stationary third-order-dispersion soliton, the radiative "
            "loss of the soliton core from direct integration of the generalized nonlinear Schroedinger equation, the kinematic "
            "surface gravity and its thermal and entanglement layers, and the data products and figures. It is a Python "
            "(NumPy, Matplotlib) rewrite of the project's MATLAB programs and runs on a laptop or on a SLURM cluster.</p>"
            "<p>Files:</p><ul>%s</ul>"
            "<p>Unpack the tarball and read <code>README.md</code>: <code>python p2.py run</code> runs every stage on one "
            "computer; <code>README_UHeM_Altay.md</code> describes a SLURM run. <code>CHECKSUMS.sha256</code> lists the SHA-256 "
            "of every file of this record.</p>"
            "<p>Source repository: <a href=\"%s/tree/p2-v%s\">%s</a>, tag <code>p2-v%s</code>. The code builds on the package of "
            "H. Oguz, Class. Quantum Grav. 43, 135014 (2026), doi:%s, archived at doi:%s.</p>"
            % (ver, "".join(items), REPO, ver, REPO, ver, REF25_PAPER_DOI, REF25_DOI))


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--results", action="append", default=[])
    ap.add_argument("--extra", action="append", default=[])
    ap.add_argument("--doi", default="")
    ap.add_argument("--date", default=datetime.datetime.now(datetime.timezone.utc).strftime("%Y-%m-%d"))
    ap.add_argument("--out", default=os.path.join(B.ROOT, "dist"))
    a = ap.parse_args(argv)
    ver = B.version()
    if a.doi:
        B.set_doi(a.doi)
    tarball = B.build(a.out)
    zdir = os.path.join(a.out, "zenodo_SOTHE-P2_v%s" % ver)
    if os.path.isdir(zdir):
        shutil.rmtree(zdir)
    os.makedirs(zdir)
    shutil.copyfile(tarball, os.path.join(zdir, os.path.basename(tarball)))
    shutil.copyfile(tarball + ".sha256", os.path.join(zdir, os.path.basename(tarball) + ".sha256"))
    results = [add_results(r, zdir) for r in a.results]
    extras = []
    for e in a.extra:
        if not os.path.isfile(e):
            sys.exit("--extra: no such file: %s" % e)
        extras.append(shutil.copyfile(e, os.path.join(zdir, os.path.basename(e))))
    meta = {"title": "SOTHE-P2 %s: code and data for the Paper-2 computations of the SOTHE project" % ver,
            "upload_type": "software", "publication_date": a.date, "version": ver, "language": "eng",
            "creators": [CREATOR], "description": description(ver, results, extras), "access_right": "open",
            "license": "MIT", "keywords": KEYWORDS,
            "related_identifiers": [
                {"identifier": "%s/tree/p2-v%s" % (REPO, ver), "relation": "isSupplementTo", "scheme": "url",
                 "resource_type": "software"},
                {"identifier": REF25_DOI, "relation": "references", "scheme": "doi", "resource_type": "software"},
                {"identifier": REF25_PAPER_DOI, "relation": "references", "scheme": "doi", "resource_type": "publication-article"}],
            "notes": "Add the Paper-2 article DOI (relation: isSupplementTo) with 'New version' or 'Edit' once it is published."}
    if a.doi:
        meta["doi"] = a.doi
    with open(os.path.join(zdir, "zenodo_metadata.json"), "w", encoding="utf-8", newline="\n") as f:
        json.dump({"metadata": meta}, f, indent=1, ensure_ascii=False)
        f.write("\n")
    readme = ["# SOTHE-P2 %s (Zenodo record)" % ver, "",
              "Author: Hasan Oguz (Istanbul Okan University; Pamukkale University), ORCID 0000-0001-7484-4415.", "",
              "This record holds the files listed below. The description of the record is the same text.", "",
              "| File | Content |", "|---|---|",
              "| `%s` | Python package SOTHE-P2 %s: source, tests, documentation, `data/entropy_trajectory.csv`, reference outputs of the MATLAB R2025b and Octave 8.4.0 implementations |" % (os.path.basename(tarball), ver)]
    for r in results:
        readme.append("| `%s` | output folder of one complete run (reports, CSV and JSON products, figures, `SUMMARY.md`, `MANIFEST.sha256`) |" % os.path.basename(r))
    for e in extras:
        readme.append("| `%s` | (added with --extra) |" % os.path.basename(e))
    readme += ["| `zenodo_metadata.json` | the metadata of this record |", "| `CHECKSUMS.sha256` | SHA-256 of every file |", "",
               "## Use", "", "```", "tar xzf %s" % os.path.basename(tarball), "cd SOTHE-P2", "python p2.py verify        # checks every file against MANIFEST.sha256",
               "python p2.py run           # all stages on this computer (Python >= 3.9, NumPy, Matplotlib)", "```", "",
               "SLURM clusters: `README_UHeM_Altay.md` inside the tarball.", "",
               "## Source and related records", "",
               "- Repository: %s, tag `p2-v%s` (branch `p2`)." % (REPO, ver),
               "- Builds on H. Oguz, *Generalized thermodynamics of solitonic event horizons in dispersive field theories*, "
               "Class. Quantum Grav. **43**, 135014 (2026), doi:%s; code and data doi:%s." % (REF25_PAPER_DOI, REF25_DOI), "",
               "## Licence", "", "MIT (`LICENSE` in the tarball).", ""]
    with open(os.path.join(zdir, "README_ZENODO.md"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(readme))
    steps = ["Manual Zenodo upload of SOTHE-P2 %s (own record)" % ver, "",
             "1. https://zenodo.org -> log in -> New upload (if the DOI was reserved earlier, open that draft instead).",
             "2. Files: drag in every file of this folder except UPLOAD_STEPS.txt.",
             "3. Basic information: copy the fields of zenodo_metadata.json (resource type Software, title, publication date,",
             "   creators with ORCID and affiliation, description, version %s, licence MIT, keywords)." % ver,
             "4. Related works: the three entries of related_identifiers (GitHub tag URL: 'Is supplement to';",
             "   the two DOIs: 'References').",
             "5. DOI: 'Do you already have a DOI?' -> No -> Get a DOI now! (if not reserved yet). Copy it.",
             "6. Preview, then Publish. The record shows two DOIs: this version's DOI and the concept DOI (all versions).",
             "   The manuscript token [ZENODO-VERSION-DOI] takes this version's DOI; [RELEASE-TAG] takes p2-v%s." % ver,
             "7. Later versions: open the record -> New version -> replace the files -> Publish (new version DOI, same concept DOI).",
             "", "sha256 of the files: CHECKSUMS.sha256"]
    with open(os.path.join(zdir, "UPLOAD_STEPS.txt"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(steps) + "\n")
    names = sorted(n for n in os.listdir(zdir) if n not in ("CHECKSUMS.sha256", "UPLOAD_STEPS.txt"))
    with open(os.path.join(zdir, "CHECKSUMS.sha256"), "w", encoding="utf-8", newline="\n") as f:
        for n in names:
            f.write("%s  %s\n" % (sha256(os.path.join(zdir, n)), n))
    print("\nZenodo pack: %s" % zdir)
    for n in sorted(os.listdir(zdir)):
        print("  %-48s %10.1f kB" % (n, os.path.getsize(os.path.join(zdir, n)) / 1024))
    print("\n" + "\n".join(steps))
    return zdir


if __name__ == "__main__":
    main()
