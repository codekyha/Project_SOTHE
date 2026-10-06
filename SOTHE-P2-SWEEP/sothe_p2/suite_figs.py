"""The six figures of sothe_phase2_matlab v0.1.0 (stage S6), rendered for publication (pubstyle.py, IOP guidelines).

  fig1_thermality            fig1_thermality.m 0.1.0            (a) ln|beta/alpha|^2 and its fit; (b) kappa against drive
  fig2_bdg_krein             fig2_bdg_krein.m 0.1.0             Krein-signed BdG spectra about psi0, delta3 = 0 and 0.02
  fig3_gapless_flow          fig3_gapless_flow.m 0.1.0          kappa against kappa_g, kappa = N at kappa_g = 0, T_H axis
  fig4_false_positive        fig4_false_positive.m 0.1.0        |beta| = |tanh(kappa_g d)| with contours
  fig5_stability_map         fig5_stability_map.m 0.1.0         max|Im omega| over (N, delta3) with contours
  fig6_partner_entanglement  fig6_partner_entanglement.m 0.1.0  E_N(omega) and its closed-form check
  magma, viridis             p2_magma.m, p2_viridis.m           the polynomial colormap fits (public domain)
  label_contours             p2_label_contours.m                one label at the middle vertex of every contour piece
Same data, symbols, marker shapes and colours as the MATLAB scripts.  Following the PI's instruction (2026-09-27; also
DELTA G3i), the figures carry no title and no annotation: the titles, text boxes and in-plot comments of the MATLAB
scripts, and the cyan null-locus lines of fig4 (on the frame, at kappa_g = 0 and d = 0), are left out, and
captions.md / captions.tex explain each figure.  Contours are labelled inline.  Where a suite script leaves the axis
limits automatic (fig1, fig3, fig6), the figure takes those of make_paper2_figs.m 0.2.1 (in fig3, the star at kappa = N
then sits inside the frame); limits that a MATLAB script sets are kept.
fig2 (b) and fig5 linearize about psi0 (G3u audit, S1-1): regression products of suite v0.1.0."""
import math
import traceback

import numpy as np

from . import compat as C
from . import pubstyle as S


def magma(n=256):
    """p2_magma.m: magma colormap, public-domain polynomial fit (Matt DesLauriers), clipped to [0, 1]."""
    t = C.mlinspace(0, 1, n).reshape(-1, 1)
    c0 = np.array([-0.002136485053939, -0.000749655052795, -0.005386127855323])
    c1 = np.array([0.251364492735802, 0.677988825427858, 2.494026927261818])
    c2 = np.array([8.353717279216625, -3.577719514958484, 0.319563623025944])
    c3 = np.array([-27.66873308576866, 14.26473078096533, -13.64921318813922])
    c4 = np.array([52.17613981234068, -27.94360607168351, 12.94416944238394])
    c5 = np.array([-50.76852536473588, 29.04658282127291, 4.23415299384598])
    c6 = np.array([18.65570506591883, -11.48977351997711, -5.601961508734096])
    Cm = c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6)))))
    return np.minimum(np.maximum(Cm, 0), 1)


def viridis(n=256):
    """p2_viridis.m: viridis colormap, public-domain polynomial fit (Matt DesLauriers), clipped to [0, 1]."""
    t = C.mlinspace(0, 1, n).reshape(-1, 1)
    c0 = np.array([0.2777273272234177, 0.005407344544966578, 0.3340998053353061])
    c1 = np.array([0.1050930431085774, 1.404613529898575, 1.384590162594685])
    c2 = np.array([-0.3308618287255563, 0.214847559468213, 0.09509516302823659])
    c3 = np.array([-4.634230498983486, -5.799100973351585, -19.33244095627987])
    c4 = np.array([6.228269936347081, 14.17993336680509, 56.69055260068105])
    c5 = np.array([4.776384997670288, -13.74514537774601, -65.35303263337234])
    c6 = np.array([-5.435455855934631, 4.645852612178535, 26.3124352495832])
    Cm = c0 + t * (c1 + t * (c2 + t * (c3 + t * (c4 + t * (c5 + t * c6)))))
    return np.minimum(np.maximum(Cm, 0), 1)


def _cmap(rgb, name):
    from matplotlib.colors import ListedColormap
    return ListedColormap(rgb, name=name)


def contour_pieces(cs):
    """[(level, [vertex arrays of its pieces])], the content of MATLAB's contour matrix."""
    import warnings
    try:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            segs = cs.allsegs
        return list(zip(cs.levels, segs))
    except AttributeError:                            # Matplotlib versions without allsegs: split the level paths
        from matplotlib.path import Path
        out = []
        for lvl, path in zip(cs.levels, cs.get_paths()):
            v = path.vertices
            if path.codes is None:
                out.append((lvl, [v] if len(v) else []))
                continue
            starts = list(np.flatnonzero(path.codes == Path.MOVETO)) + [len(v)]
            out.append((lvl, [v[a:b] for a, b in zip(starts[:-1], starts[1:]) if b > a]))
        return out


def label_contours(ax, cs, fmt="%.2f", col="w"):
    """p2_label_contours.m: one label per contour piece, at its middle vertex (j = max(1, round(np/2))).
    Kept as the port of the MATLAB helper; the publication figures label contours inline (clabel)."""
    for lvl, segs in contour_pieces(cs):
        for seg in segs:
            npt = len(seg)
            if npt == 0:
                continue
            j = max(1, int(C.mround(npt / 2.0)))
            ax.text(seg[j - 1][0], seg[j - 1][1], fmt % lvl, color=col, fontsize=S.FS_SMALL, ha="center", va="center")


# ---------------------------------------------------------------- the six figures
def fig1_thermality(R, outdir):
    plt = S.setup()
    f = S.figure(plt, S.W2, 5.4)
    a, b = f.subplots(1, 2)
    a.plot(R["omega"], np.log(R["ratio"]), "o", ms=3.0, mfc=S.ORANGE, mec=S.ORANGE, ls="none", clip_on=False, zorder=3)
    xl = np.array([np.min(R["omega"]), np.max(R["omega"])])
    a.plot(xl, R["slope"] * xl + R["intercept"], "k-", lw=0.8, zorder=4)   # the line over the 60 markers
    a.set_xlim(0, 1.65)
    a.set_xlabel(r"$\omega$")
    a.set_ylabel(r"$\ln|\beta_\omega/\alpha_\omega|^2$")
    S.label(a, "(a)")
    S.frame(a)
    b.plot(R["drive"], R["kappa_drive"], "s", mfc=S.BLUE, mec=S.BLUE, ls="none", clip_on=False, zorder=3)
    m = float(np.mean(R["kappa_drive"]))
    b.set_ylim(m - 0.06, m + 0.06)
    b.set_xlim(0.65, 1.35)
    b.set_xlabel(r"drive amplitude $d$")
    b.set_ylabel(r"$\kappa$")
    S.label(b, "(b)")
    S.frame(b)
    return S.save(plt, f, outdir, "fig1_thermality")


def fig2_bdg_krein(R, outdir):
    plt = S.setup()
    f = S.figure(plt, S.W2, 5.8)
    axes = f.subplots(1, 2)
    for p, (a, tag) in enumerate(zip(axes, ("(a)", "(b)"))):
        o = R["bdg"][p]["out"]
        mu = o["mu"]
        w = o["w"][o["win"]]
        kr = o["krein"][o["win"]]
        ym = 1.15 * max(1e-3, float(np.max(np.abs(w.imag))) if w.size else 1e-3)
        hb = a.fill([-mu, -mu, mu, mu], [-ym, ym, ym, -ym], color=S.BAND, ec="none", zorder=0)[0]
        hp, = a.plot(w[kr > 0].real, w[kr > 0].imag, "o", mec=S.KPOS, mfc="none", ls="none", zorder=3)
        hn, = a.plot(w[kr < 0].real, w[kr < 0].imag, "D", ms=S.MS - 0.5, mec=S.KNEG, mfc="none", ls="none", zorder=3)
        a.set_xlim(-3 * mu, 3 * mu)
        a.set_ylim(-ym, ym)
        a.set_xlabel(r"$\mathrm{Re}\,\omega$")
        S.frame(a)
        S.scaled_axis(a, "y", r"$\mathrm{Im}\,\omega$")
        S.label(a, tag)
        if p == 0:
            a.legend([hp, hn, hb], [r"Krein $+$", r"Krein $-$", r"$|\mathrm{Re}\,\omega|<\mu$"], loc="lower left")
    return S.save(plt, f, outdir, "fig2_bdg_krein")


def fig3_gapless_flow(R, P, outdir):
    plt = S.setup()
    f = S.figure(plt, S.W1, 6.2)
    a = f.add_subplot(1, 1, 1)
    a.plot(R["kg_flow"], R["kappa_flow"], "k-", zorder=2)
    a.plot(R["kg_flow"], R["kappa_flow"], "o", mfc=S.ORANGE, mec=S.ORANGE, ls="none", clip_on=False, zorder=3)
    xr = [float(np.min(R["kg_flow"])), float(np.max(R["kg_flow"]))]
    a.plot(xr, [P["N"], P["N"]], "k:", lw=0.8)
    a.plot([0], [P["N"]], "*", ms=7.5, mfc=S.BLUE, mec="k", mew=0.5, ls="none", zorder=4)   # MATLAB 'p' (pentagram)
    pad = 0.03 * (xr[1] - xr[0])                        # limits of make_paper2_figs.m 0.2.1: the star and the line
    a.set_xlim(xr[0] - pad, xr[1] + pad)                # at kappa = N stay inside the frame, off the corner tick
    yl = (0.9 * float(np.min(R["kappa_flow"])), 1.04 * P["N"])
    a.set_ylim(*yl)
    a.set_xlabel(r"$\kappa_g$")
    a.set_ylabel(r"$\kappa$")
    S.frame(a, mirror=False)
    a.tick_params(which="both", top=True)
    a2 = a.twinx()
    a2.set_ylim(yl[0] / (2 * math.pi), yl[1] / (2 * math.pi))
    a2.set_ylabel(r"$T_H=\kappa/2\pi$")
    a2.tick_params(which="both", direction="in", right=True, labelsize=S.FS_SMALL)
    S.nice_ticks(a2)
    for sp in a2.spines.values():
        sp.set_linewidth(S.LW_AX)
    return S.save(plt, f, outdir, "fig3_gapless_flow")


def _map(x, y, Z, cmap, levels, outdir, name, xlabel, ylabel, cblabel):
    plt = S.setup()
    f = S.figure(plt, S.W1, 6.4)
    a = f.add_subplot(1, 1, 1)
    Z = np.asarray(Z, dtype=float)
    pc = a.pcolormesh(x, y, Z, shading="gouraud", cmap=cmap, vmin=0, vmax=float(np.max(Z)))       # vector Gouraud shading
    cs = a.contour(x, y, Z, levels=levels, colors="w", linewidths=0.6)
    a.clabel(cs, inline=True, fmt="%.2f", fontsize=S.FS_SMALL, inline_spacing=3)
    a.set_xlim(float(np.min(x)), float(np.max(x)))
    a.set_ylim(float(np.min(y)), float(np.max(y)))
    a.set_xlabel(xlabel)
    a.set_ylabel(ylabel)
    S.frame(a)
    cb = f.colorbar(pc, ax=a, fraction=0.06, pad=0.03)
    cb.set_label(cblabel)
    cb.outline.set_linewidth(S.LW_AX)
    cb.ax.tick_params(direction="in", labelsize=S.FS_SMALL, width=S.LW_AX)
    return S.save(plt, f, outdir, name)


def fig4_false_positive(R, outdir):
    return _map(np.asarray(R["fp_kg"], dtype=float), np.asarray(R["fp_drive"], dtype=float), R["fp_beta"],
                _cmap(magma(256), "p2_magma"), [0.05, 0.10, 0.20, 0.30, 0.40], outdir, "fig4_false_positive",
                r"$\kappa_g$", r"drive amplitude $d$", r"$|\beta|$")


def fig5_stability_map(R, outdir):
    return _map(np.asarray(R["stab_N"], dtype=float), np.asarray(R["stab_d3"], dtype=float), R["stab_map"],
                _cmap(viridis(256), "p2_viridis"), [0.05, 0.10, 0.20, 0.40, 0.80], outdir, "fig5_stability_map",
                r"$N$", r"$\delta_3$", r"$\max|\mathrm{Im}\,\omega|$")


def fig6_partner_entanglement(R, outdir):
    plt = S.setup()
    f = S.figure(plt, S.W1, 5.6)
    a = f.add_subplot(1, 1, 1)
    om = np.asarray(R["omega"], dtype=float).ravel()
    EN = np.asarray(R["E_N"], dtype=float).ravel()
    a.fill(np.concatenate([om, om[::-1]]), np.concatenate([EN, np.zeros(om.size)]), color=S.FILL, ec="none", zorder=0)
    h1, = a.plot(om, EN, "-", color=S.ORANGE, lw=1.4)
    h2, = a.plot(om, R["E_N_cross"], "k--", lw=0.9)
    a.set_xlim(0, 1.65)
    a.set_ylim(0, 1.08 * float(np.max(EN)))
    a.set_xlabel(r"$\omega$")
    a.set_ylabel(r"$E_N$ (nats)")
    a.legend([h1, h2], [r"$E_N=-\ln\nu_-$", r"$2\ln(|\alpha_\omega|+|\beta_\omega|)$"], loc="upper right")
    S.frame(a)
    return S.save(plt, f, outdir, "fig6_partner_entanglement")


FIGURES = ("fig1_thermality", "fig2_bdg_krein", "fig3_gapless_flow", "fig4_false_positive", "fig5_stability_map",
           "fig6_partner_entanglement")


def captions(R, P):
    """Caption texts (explanation only, no analysis), LaTeX mathematics."""
    op = r"$N=%g$, $\delta_3=%g$, $\kappa_g=%g$, $k_{\rm op}=%g$" % (P["N"], P["d3"], P["kg"], P["k_op"])
    nw, w0, w1 = len(R["omega"]), float(R["omega"][0]), float(R["omega"][-1])
    grid = r"a periodic grid of $n_\tau=%d$ points, $\tau\in[-%g,%g)$, fourth-order finite differences" % (
        P["bdg_n_tau"], P["bdg_tau_max"], P["bdg_tau_max"])
    return [
        ("fig1_thermality", "Figure 1",
         r"(a) $\ln|\beta_\omega/\alpha_\omega|^2$ for $|\beta_\omega/\alpha_\omega|^2=e^{-2\pi\omega/\kappa}$ at the kinematic "
         r"surface gravity $\kappa=%.4f$, at %d frequencies $\omega\in[%.2f,%.2f]$ (circles), and its unweighted least-squares "
         r"straight line (solid). (b) Kinematic surface gravity $\kappa$ at %d drive amplitudes $d\in[%.1f,%.1f]$. "
         r"Parameters: %s." % (R["kappa_kin"], nw, w0, w1, len(R["drive"]), float(R["drive"][0]), float(R["drive"][-1]), op)),
        ("fig2_bdg_krein", "Figure 2",
         r"Bogoliubov–de Gennes eigenvalues $\omega$ for the linearization about $\psi_0=N\,\mathrm{sech}(N\tau)$, $N=%g$, on %s, "
         r"in the window $|\mathrm{Re}\,\omega|\le3\mu$, $\mu=N^2/2$: (a) $\delta_3=%g$; (b) $\delta_3=%g$. Open circles: "
         r"positive Krein signature; open diamonds: negative Krein signature; shaded band: $|\mathrm{Re}\,\omega|<\mu$."
         % (P["bdg_N"], grid, R["bdg"][0]["d3"], R["bdg"][1]["d3"])),
        ("fig3_gapless_flow", "Figure 3",
         r"Kinematic surface gravity $\kappa$ against the slow-light coupling gap $\kappa_g$ at %d values "
         r"$\kappa_g\in[%g,%g]$ (circles joined by a line), for $N=%g$, $\delta_3=%g$, $k_{\rm op}=%g$, $d=%g$. Star and dotted "
         r"line: $\kappa=N$. Right-hand axis: $T_H=\kappa/2\pi$."
         % (len(R["kg_flow"]), float(R["kg_flow"][0]), float(R["kg_flow"][-1]), P["N"], P["d3"], P["k_op"], P["drive0"])),
        ("fig4_false_positive", "Figure 4",
         r"$|\beta|=|\tanh(\kappa_g d)|$ on a $%d\times%d$ grid, $\kappa_g\in[%g,%g]$, drive amplitude $d\in[%g,%g]$, "
         r"with contours at $|\beta|=0.05$, 0.10, 0.20, 0.30 and 0.40."
         % (len(R["fp_kg"]), len(R["fp_drive"]), float(R["fp_kg"][0]), float(R["fp_kg"][-1]), float(R["fp_drive"][0]),
            float(R["fp_drive"][-1]))),
        ("fig5_stability_map", "Figure 5",
         r"Largest imaginary part $\max|\mathrm{Im}\,\omega|$ of the Bogoliubov–de Gennes eigenvalues about "
         r"$\psi_0=N\,\mathrm{sech}(N\tau)$ in the annulus $0.1\mu<|\mathrm{Re}\,\omega|\le3\mu$, on a $%d\times%d$ grid "
         r"$N\in[%g,%g]$, $\delta_3\in[%g,%g]$ (periodic grid of $n_\tau=%d$ points, $\tau\in[-%g,%g)$), with contours at "
         r"0.05, 0.10, 0.20, 0.40 and 0.80."
         % (len(R["stab_N"]), len(R["stab_d3"]), float(R["stab_N"][0]), float(R["stab_N"][-1]), float(R["stab_d3"][0]),
            float(R["stab_d3"][-1]), P["stab_n_tau"], P["stab_tau_max"], P["stab_tau_max"])),
        ("fig6_partner_entanglement", "Figure 6",
         r"Logarithmic negativity of the Hawking partners, $E_N=2\,\mathrm{arcsinh}|\beta_\omega|=-\ln\nu_-$ (solid), and "
         r"$2\ln(|\alpha_\omega|+|\beta_\omega|)$ (dashed), with $|\beta_\omega|^2=(e^{2\pi\omega/\kappa}-1)^{-1}$ and "
         r"$|\alpha_\omega|^2=1+|\beta_\omega|^2$, at %d frequencies $\omega\in[%.2f,%.2f]$ for $\kappa=%.4f$."
         % (nw, w0, w1, R["kappa_kin"]))]


def render_all(R, P, fdir):
    """run_phase2_all.m figure block: each figure guarded (a renderer fault must not lose data), then the captions.
    Returns {"figures": [png], "files": [every file], "errors": {figure: message}}."""
    jobs = {"fig1_thermality": lambda: fig1_thermality(R, fdir), "fig2_bdg_krein": lambda: fig2_bdg_krein(R, fdir),
            "fig3_gapless_flow": lambda: fig3_gapless_flow(R, P, fdir), "fig4_false_positive": lambda: fig4_false_positive(R, fdir),
            "fig5_stability_map": lambda: fig5_stability_map(R, fdir),
            "fig6_partner_entanglement": lambda: fig6_partner_entanglement(R, fdir)}
    pngs, files, errs = [], [], {}
    for name in FIGURES:
        try:
            out = jobs[name]()
            pngs.append(name + ".png")
            files += out
        except Exception:
            errs[name] = traceback.format_exc().strip().splitlines()[-1]
    try:
        files += S.write_captions(fdir, captions(R, P), "Phase-2 data products (suite v0.1.0): figure captions")
    except Exception:
        errs["captions"] = traceback.format_exc().strip().splitlines()[-1]
    return {"figures": pngs, "files": files, "errors": errs}
