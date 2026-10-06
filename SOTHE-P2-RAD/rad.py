#!/usr/bin/env python3
"""SOTHE-P2-RAD command line.

  python rad.py list                                  the run catalogue
  python rad.py run --out RESULTS [--set S ...] [--workers W]
                                                      integrate the catalogue (or the groups/ids S) into RESULTS/<UTC stamp>/
  python rad.py analyse --results RESULTS/<stamp>     rates, universal rate function, recoil model, tables, SUMMARY.md
  python rad.py figures --results RESULTS/<stamp>     the article figure and the diagnostic figures
  python rad.py manifest --results RESULTS/<stamp>    write MANIFEST.sha256 of a results folder
  python rad.py verify [--results RESULTS/<stamp>]    check MANIFEST.sha256 of this package (and of a results folder)
  python rad.py test                                  unit tests (tests/test_rad.py)

Groups of the catalogue: table, prepared, recoil, scale, check (sothe_rad/experiments.py)."""
import argparse
import os
import subprocess
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.dont_write_bytecode = True


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list")
    r = sub.add_parser("run")
    r.add_argument("--out", required=True)
    r.add_argument("--set", nargs="*", default=["all"])
    r.add_argument("--workers", type=int, default=1)
    r.add_argument("--stamp", default=None)
    a = sub.add_parser("analyse")
    a.add_argument("--results", required=True)
    f = sub.add_parser("figures")
    f.add_argument("--results", required=True)
    m = sub.add_parser("manifest")
    m.add_argument("--results", required=True)
    v = sub.add_parser("verify")
    v.add_argument("--results", default=None)
    sub.add_parser("test")
    args = ap.parse_args()

    if args.cmd == "list":
        from sothe_rad import experiments
        for run in experiments.catalogue():
            extra = {k: v for k, v in run.items() if k not in ("id", "group", "N", "d3")}
            print("%-34s %-9s N=%-4g d3=%-6g %s" % (run["id"], run["group"], run["N"], run["d3"], extra))
        return 0
    if args.cmd == "run":
        from sothe_rad import runner
        folder = runner.run_catalogue(os.path.abspath(args.out), args.set, args.workers, args.stamp, " ".join(["rad.py"] + sys.argv[1:]))
        print("results: %s" % folder)
        return 0
    if args.cmd == "analyse":
        from sothe_rad import analysis
        analysis.analyse(os.path.abspath(args.results))
        return 0
    if args.cmd == "figures":
        from sothe_rad import figures
        figures.make_all(os.path.abspath(args.results))
        return 0
    if args.cmd == "manifest":
        from sothe_rad import runner
        n = runner.write_manifest(os.path.abspath(args.results))
        print("MANIFEST.sha256: %d files" % n)
        return 0
    if args.cmd == "verify":
        from sothe_rad import runner
        ok = True
        n, bad, extra = runner.check_manifest(HERE)
        print("package: MANIFEST: %d files, %d changed or missing, %d not listed" % (n, len(bad), len(extra)))
        ok &= not bad and not extra
        if args.results:
            n, bad, extra = runner.check_manifest(os.path.abspath(args.results))
            print("results: MANIFEST: %d files, %d changed or missing, %d not listed" % (n, len(bad), len(extra)))
            ok &= not bad and not extra
        print("intact" if ok else "NOT intact")
        return 0 if ok else 1
    if args.cmd == "test":
        return subprocess.call([sys.executable, "-B", os.path.join(HERE, "tests", "test_rad.py")], cwd=HERE)
    return 2


if __name__ == "__main__":
    sys.exit(main())
