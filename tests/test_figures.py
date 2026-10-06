"""The publication style (pubstyle.py): scaled small-value axes, reproducible PDF/EPS/PNG files, the caption files,
and the six suite figures rendered from a small physics pass."""
import hashlib
import os
import tempfile
import unittest

import _util

try:
    import matplotlib  # noqa: F401
    HAVE_MPL = True
except ImportError:
    HAVE_MPL = False


def sha(path):
    with open(path, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


@unittest.skipUnless(HAVE_MPL, "matplotlib is not installed")
class TestPubstyle(unittest.TestCase):
    def setUp(self):
        from sothe_p2 import pubstyle
        self.S = pubstyle
        self.plt = pubstyle.setup()

    def test_scaled_axis(self):
        S, plt = self.S, self.plt
        f = S.figure(plt, S.W1, 6.0)
        a = f.add_subplot(1, 1, 1)
        a.plot([-1, 1], [-1e-3, 1e-3])
        a.set_ylim(-1.15e-3, 1.15e-3)
        S.frame(a)
        S.scaled_axis(a, "y", r"$\mathrm{Im}\,\Omega$")
        self.assertEqual(a.get_ylabel(), r"$10^{3}\,$$\mathrm{Im}\,\Omega$")
        f.canvas.draw()
        labels = [t.get_text() for t in a.get_yticklabels() if t.get_text()]
        self.assertEqual(labels, [r"$\mathdefault{-1}$", r"$\mathdefault{0}$", r"$\mathdefault{1}$"])
        self.assertEqual(a.yaxis.get_offset_text().get_text(), "")
        b = f.add_subplot(2, 1, 2)
        b.set_ylim(-0.11, 0.11)
        S.scaled_axis(b, "y", r"$\mathrm{Im}\,\omega$")           # |values| >= 0.01: numbers kept, plain label
        self.assertEqual(b.get_ylabel(), r"$\mathrm{Im}\,\omega$")
        plt.close(f)

    def test_nice_ticks(self):
        from matplotlib.ticker import MaxNLocator
        S, plt = self.S, self.plt
        f = S.figure(plt, S.W1, 6.0)
        a = f.add_subplot(1, 1, 1)
        a.plot([0, 1], [0, 1])
        S.frame(a)
        self.assertIsInstance(a.xaxis.get_major_locator(), MaxNLocator)
        steps = [t for t in a.get_yticks()]
        d = sorted({round(steps[i + 1] - steps[i], 12) for i in range(len(steps) - 1)})
        self.assertEqual(len(d), 1)
        m = d[0] / 10 ** int(__import__("math").floor(__import__("math").log10(d[0])))
        self.assertIn(round(m, 9), (1.0, 2.0, 5.0))
        plt.close(f)

    def test_save_is_reproducible(self):
        import numpy as np
        S = self.S
        hashes = []
        with tempfile.TemporaryDirectory() as d:
            for k in range(2):
                plt = S.setup()
                f = S.figure(plt, S.W1, 5.6)
                a = f.add_subplot(1, 1, 1)
                x = np.linspace(0, 1.6, 60)
                a.plot(x, np.exp(-x), "-", color=S.ORANGE)
                a.set_xlabel(r"$\omega$")
                a.set_ylabel(r"$E_N$ (nats)")
                S.frame(a)
                S.label(a, "(a)")
                sub = os.path.join(d, str(k))
                os.makedirs(sub)
                names = S.save(plt, f, sub, "t")
                self.assertEqual(names, ["t.png", "t.pdf", "t.eps"])
                hashes.append([sha(os.path.join(sub, n)) for n in names])
                with open(os.path.join(sub, "t.pdf"), "rb") as fh:
                    self.assertNotIn(b"/CreationDate", fh.read())
        self.assertEqual(hashes[0], hashes[1])

    def test_font_choice(self):
        import os as _os
        old = _os.environ.pop("SOTHE_P2_FIG_FONT", None)
        try:
            self.assertEqual(self.S.font_choice(), "cm")
            _os.environ["SOTHE_P2_FIG_FONT"] = "STIX"
            self.assertEqual(self.S.font_choice(), "stix")
            self.S.setup()
            from matplotlib import rcParams
            self.assertEqual(rcParams["mathtext.fontset"], "stix")
            _os.environ["SOTHE_P2_FIG_FONT"] = "comic"
            self.assertEqual(self.S.font_choice(), "cm")               # unknown values fall back to Computer Modern
        finally:
            _os.environ.pop("SOTHE_P2_FIG_FONT", None)
            if old is not None:
                _os.environ["SOTHE_P2_FIG_FONT"] = old
            self.S.setup()

    def test_captions(self):
        with tempfile.TemporaryDirectory() as d:
            files = self.S.write_captions(d, [("fig1_x", "Figure 1", r"(a) $\kappa$ against $\omega$.")], "Test")
            self.assertEqual(files, ["captions.md", "captions.tex"])
            tex = _util.text(os.path.join(d, "captions.tex"))
            self.assertIn("% Figure 1: fig1_x\n\\caption{(a) $\\kappa$ against $\\omega$.}\n", tex)
            md = _util.text(os.path.join(d, "captions.md"))
            self.assertIn("## Figure 1 (`fig1_x`)", md)


@unittest.skipUnless(HAVE_MPL, "matplotlib is not installed")
@unittest.skipIf(os.environ.get("SOTHE_P2_QUICK"), "SOTHE_P2_QUICK set: no physics pass")
class TestSuiteFigures(unittest.TestCase):
    def test_render_all(self):
        from sothe_p2 import suite, suite_figs
        P = suite.phase2_params(suite.SMALL_STAB)
        R = suite.compute_all(P, _util.TRAJ)
        with tempfile.TemporaryDirectory() as d:
            out = suite_figs.render_all(R, P, d)
            self.assertEqual(out["errors"], {})
            self.assertEqual(out["figures"], [n + ".png" for n in suite_figs.FIGURES])
            for n in suite_figs.FIGURES:
                for ext in (".pdf", ".eps", ".png"):
                    self.assertTrue(os.path.getsize(os.path.join(d, n + ext)) > 1000, n + ext)
            tex = _util.text(os.path.join(d, "captions.tex"))
            self.assertEqual(sum(1 for line in tex.splitlines() if line.startswith("\\caption{")), 6)


if __name__ == "__main__":
    unittest.main()
