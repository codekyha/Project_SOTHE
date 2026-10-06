"""The command line: help, stages, version, a dry-run submission (no SLURM needed), the MANIFEST check."""
import glob
import os
import subprocess
import sys
import tempfile
import unittest

import _util


def p2(*args, env=None):
    e = dict(os.environ)
    e.update(env or {})
    p = subprocess.run([sys.executable, os.path.join(_util.ROOT, "p2.py")] + list(args), stdout=subprocess.PIPE,
                       stderr=subprocess.STDOUT, universal_newlines=True, env=e, cwd=_util.ROOT, timeout=120)
    return p.returncode, p.stdout


class TestCli(unittest.TestCase):
    def test_help_version_stages(self):
        rc, out = p2("help")
        self.assertEqual(rc, 0)
        self.assertIn("SOTHE-P2", out)
        rc, out = p2("version")
        self.assertEqual(rc, 0)
        self.assertIn(_util.text(os.path.join(_util.ROOT, "VERSION")).strip(), out)
        rc, out = p2("stages")
        self.assertEqual(rc, 0)
        for s in ("S0", "S1", "S2", "S3", "S4", "S5", "S6", "S7"):
            self.assertIn(s, out)

    def test_dry_run_submit(self):
        with tempfile.TemporaryDirectory() as d:
            rc, out = p2("submit", "--dry-run", "--stages", "S1 S6", env={"SOTHE_P2_RUNS": d})
            self.assertEqual(rc, 0, out)
            scripts = glob.glob(os.path.join(d, "dryrun_*", "job.sbatch"))
            self.assertEqual(len(scripts), 1)
            text = _util.text(scripts[0])
            self.assertIn("#SBATCH -A", text)
            self.assertIn("p2.py _job", text)
            self.assertIn("S1 S6", text)

    def test_settings_recorded(self):
        with tempfile.TemporaryDirectory() as d:
            rc, out = p2("submit", "--dry-run", "--stages", "S5", env={"SOTHE_P2_RUNS": d, "SOTHE_P2_FIG_FONT": "stix",
                                                                     "PHASE2_STAB_NTAU": "96"})
            self.assertEqual(rc, 0, out)
            run = glob.glob(os.path.join(d, "dryrun_*"))[0]
            job = _util.text(os.path.join(run, "job.sbatch"))
            self.assertIn("export SOTHE_P2_FIG_FONT=stix", job)
            self.assertIn("export PHASE2_STAB_NTAU=96", job)
            env = _util.text(os.path.join(run, "RUN.env"))
            self.assertIn("SOTHE_P2_FIG_FONT=stix\n", env)

    def test_results_dir_created(self):
        from sothe_p2 import cli
        with tempfile.TemporaryDirectory() as d:
            target = os.path.join(d, "new", "results")
            old = os.environ.get("SOTHE_P2_RESULTS_DIR")
            os.environ["SOTHE_P2_RESULTS_DIR"] = target
            try:
                C = cli.load_config()
            finally:
                if old is None:
                    os.environ.pop("SOTHE_P2_RESULTS_DIR", None)
                else:
                    os.environ["SOTHE_P2_RESULTS_DIR"] = old
            self.assertEqual(C["RESULTS_DIR"], target)
            self.assertTrue(os.path.isdir(target))

    def test_unknown_stage(self):
        with tempfile.TemporaryDirectory() as d:
            rc, out = p2("submit", "--dry-run", "--stages", "S9", env={"SOTHE_P2_RUNS": d})
            self.assertNotEqual(rc, 0)

    def test_verify(self):
        if not os.path.exists(os.path.join(_util.ROOT, "MANIFEST.sha256")):
            self.skipTest("no MANIFEST.sha256 (built by tools/build_release.py)")
        rc, out = p2("verify")
        self.assertEqual(rc, 0, out)


if __name__ == "__main__":
    unittest.main()
