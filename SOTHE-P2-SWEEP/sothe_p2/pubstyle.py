"""Publication style of every figure of the package (stages S5 and S6), after IOP Publishing's figure guidelines.

  size      final width 8.5 cm (one column) or 15 cm (two columns); the figure is saved at exactly that size
  text      Computer Modern (text and mathematics), 9 pt labels, 8 pt ticks, legends and contour labels;
            SOTHE_P2_FIG_FONT=stix selects the Times-like STIX fonts instead (IOP's standard families are Times,
            Helvetica, Courier and Symbol; its graphics guide also accepts Computer Modern)
  content   no title, no text boxes, no annotation: what a figure shows is explained in its caption
            (captions.md / captions.tex next to the figures); parts labelled (a), (b), (c) above their top-left corner
  colour    never the only carrier of information (marker shape or line style too); no transparency
            (tints are pre-blended, so the EPS files are exact)
  files     <name>.pdf and <name>.eps (vector), <name>.png (600 dpi)
Sources: IOP Publishing, "Figures" (publishingsupport.iopscience.iop.org/questions/figures-journal-articles) and
"Preparing graphics for IOP journals"."""
import os

FONTS = {   # Matplotlib bundles both families, so every machine renders the same glyphs
    "cm": {"font.family": "serif", "font.serif": ["cmr10", "Computer Modern Roman", "DejaVu Serif"], "mathtext.fontset": "cm"},
    "stix": {"font.family": "serif", "font.serif": ["STIXGeneral", "STIX Two Text", "Times New Roman", "DejaVu Serif"],
             "mathtext.fontset": "stix"},
}


def font_choice():
    """The lettering: 'cm' (default) or 'stix' (Times-like), from SOTHE_P2_FIG_FONT."""
    v = os.environ.get("SOTHE_P2_FIG_FONT", "cm").strip().lower()
    return v if v in FONTS else "cm"


CM = 1 / 2.54
W1 = 8.5 * CM                  # one-column figure width [in]
W2 = 15.0 * CM                 # two-column figure width [in]
FS = 9                         # labels
FS_SMALL = 8                   # ticks, legends, contour labels
LW = 1.0                       # data lines [pt]
LW_AX = 0.6                    # axes [pt]
MS = 3.5                       # markers [pt]
ORANGE = (0.90, 0.45, 0.10)
BLUE = (0.15, 0.35, 0.75)
RED = (0.80, 0.10, 0.10)
KPOS = (0.80, 0.12, 0.12)      # Krein +  (circles)
KNEG = (0.12, 0.28, 0.78)      # Krein -  (diamonds)
GREY = (0.45, 0.45, 0.45)


def tint(color, a):
    """The colour of COLOR at opacity A over white, as an opaque colour."""
    return tuple(1 - a * (1 - c) for c in color)


BAND = tint((0.62, 0.64, 0.78), 0.35)    # the |Re Omega| < mu band
FILL = tint(ORANGE, 0.22)                # area under E_N


def setup():
    """Matplotlib with the publication rcParams; returns pyplot."""
    import logging
    import matplotlib
    logging.getLogger("fontTools").setLevel(logging.ERROR)          # cmr10 has old head-table timestamps
    logging.getLogger("matplotlib.font_manager").setLevel(logging.ERROR)
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    matplotlib.rcParams.update(FONTS[font_choice()])
    matplotlib.rcParams.update({
        "axes.formatter.use_mathtext": True, "axes.unicode_minus": True, "font.size": FS, "axes.labelsize": FS,
        "axes.titlesize": FS, "xtick.labelsize": FS_SMALL, "ytick.labelsize": FS_SMALL, "legend.fontsize": FS_SMALL,
        "axes.linewidth": LW_AX, "lines.linewidth": LW, "lines.markersize": MS, "lines.markeredgewidth": 0.6,
        "xtick.direction": "in", "ytick.direction": "in", "xtick.top": True, "ytick.right": True,
        "xtick.major.size": 3.0, "ytick.major.size": 3.0, "xtick.minor.size": 1.6, "ytick.minor.size": 1.6,
        "xtick.major.width": LW_AX, "ytick.major.width": LW_AX, "xtick.minor.width": 0.5, "ytick.minor.width": 0.5,
        "legend.frameon": False, "legend.handlelength": 1.6, "legend.borderaxespad": 0.4,
        "axes.autolimit_mode": "round_numbers", "axes.xmargin": 0.0, "axes.ymargin": 0.0,
        "axes.formatter.limits": [-3, 4], "savefig.dpi": 600, "figure.dpi": 100,
        "pdf.fonttype": 42, "ps.fonttype": 42, "svg.fonttype": "none",
        "figure.constrained_layout.h_pad": 0.02, "figure.constrained_layout.w_pad": 0.02,
        "figure.constrained_layout.hspace": 0.04, "figure.constrained_layout.wspace": 0.04})
    return plt


def figure(plt, width, height_cm):
    """A figure of the final size (width in inches: W1 or W2), laid out without clipping (constrained layout)."""
    return plt.figure(figsize=(width, height_cm * CM), facecolor="w", layout="constrained")


def label(ax, s):
    """Part label (a), (b), ... above the top-left corner of the axes."""
    ax.text(0.0, 1.015, s, transform=ax.transAxes, ha="left", va="bottom", fontsize=FS)


def nice_ticks(ax):
    """Linear axes with automatic ticks: steps 1, 2 or 5 times a power of ten (never 0.25-type steps)."""
    from matplotlib.ticker import AutoLocator, MaxNLocator
    for axis, scale in ((ax.xaxis, ax.get_xscale()), (ax.yaxis, ax.get_yscale())):
        if scale == "linear" and isinstance(axis.get_major_locator(), AutoLocator):
            axis.set_major_locator(MaxNLocator(nbins="auto", steps=[1, 2, 5, 10]))


def frame(ax, mirror=True):
    for sp in ax.spines.values():
        sp.set_linewidth(LW_AX)
    ax.tick_params(which="both", direction="in", top=mirror, right=mirror)
    nice_ticks(ax)


def scaled_axis(ax, which, symbol, nbins=5):
    """Label a small-valued axis as 10^k symbol with plain tick numbers, instead of an offset text such as x10^-3
    (an offset text would collide with the part label).  Axes with |values| >= 0.01 keep their numbers."""
    import math
    from matplotlib.ticker import FixedFormatter, FixedLocator, MaxNLocator
    lo, hi = ax.get_ylim() if which == "y" else ax.get_xlim()
    setlabel = ax.set_ylabel if which == "y" else ax.set_xlabel
    m = max(abs(lo), abs(hi))
    e = int(math.floor(math.log10(m))) if m > 0 else 0
    if e > -3:
        setlabel(symbol)
        return
    f = 10.0 ** (-e)
    ticks = [t for t in MaxNLocator(nbins=nbins, steps=[1, 2, 5, 10]).tick_values(lo * f, hi * f) if lo * f - 1e-9 <= t <= hi * f + 1e-9]
    txt = ["%.6g" % round(t, 9) for t in ticks]
    d = max((len(s.split(".")[1]) if "." in s else 0) for s in txt)
    labels = [r"$\mathdefault{%s}$" % ("%.*f" % (d, t)).replace("-0.0", "0.0").replace("-0", "0") if abs(t) < 1e-12 else
              r"$\mathdefault{%s}$" % ("%.*f" % (d, t)) for t in ticks]
    axis = ax.yaxis if which == "y" else ax.xaxis
    axis.set_major_locator(FixedLocator([t / f for t in ticks]))
    axis.set_major_formatter(FixedFormatter(labels))
    setlabel(r"$10^{%d}\,$%s" % (-e, symbol))


def save(plt, fig, outdir, name):
    """<name>.pdf, <name>.eps (vector) and <name>.png (600 dpi), at the figure's final size; returns the names."""
    base = os.path.join(outdir, name)
    fig.savefig(base + ".pdf", metadata={"CreationDate": None, "ModDate": None, "Producer": None})
    old = os.environ.get("SOURCE_DATE_EPOCH")
    os.environ["SOURCE_DATE_EPOCH"] = "0"                     # a reproducible %%CreationDate in the EPS header
    try:
        fig.savefig(base + ".eps")
    finally:
        if old is None:
            os.environ.pop("SOURCE_DATE_EPOCH", None)
        else:
            os.environ["SOURCE_DATE_EPOCH"] = old
    fig.savefig(base + ".png", dpi=600, facecolor="w")
    plt.close(fig)
    return [name + ".png", name + ".pdf", name + ".eps"]


def write_captions(outdir, caps, heading):
    """captions.md (Markdown with LaTeX mathematics) and captions.tex (\\caption{...} lines) for [(file, number, text)]."""
    md = ["# %s" % heading, "",
          "Captions explain what each figure shows; the analysis is in the paper text. Mathematics is LaTeX.", ""]
    tex = ["%% %s: \\caption{} texts (explanation only; the analysis is in the paper text)" % heading]
    for name, num, text in caps:
        md += ["## %s (`%s`)" % (num, name), "", text, ""]
        tex += ["", "%% %s: %s" % (num, name), "\\caption{%s}" % text]
    with open(os.path.join(outdir, "captions.md"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(md))
    with open(os.path.join(outdir, "captions.tex"), "w", encoding="utf-8", newline="\n") as f:
        f.write("\n".join(tex) + "\n")
    return ["captions.md", "captions.tex"]
