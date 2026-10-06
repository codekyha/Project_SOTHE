"""Stage S5 (T-13 + T-18): the five manuscript figures, port of make_paper2_figs.m 0.2.1 (matlab/pack) in the
publication style of pubstyle.py (IOP guidelines: 8.5 cm or 15 cm wide, Computer Modern 8-9 pt, no titles, no
annotation text, parts (a), (b), (c), vector PDF and EPS plus 600-dpi PNG).

  fig1_thermality            Fig. 1  (a) ln|beta/alpha|^2 on 60 frequencies + least-squares line; (b) kappa at nine drives
  fig3_gapless_flow          Fig. 2  kappa over 21 kappa_g in [0, 1]; kappa = N at kappa_g = 0; right axis T_H = kappa/2pi
  fig2_bdg_krein             Fig. 3  Krein-signed BdG spectra: (a) delta3 = 0 about psi0; (b) delta3 = 0.02 about the
                                     stationary TOD soliton phi (T-13), inset: the zero-mode sector about phi and psi0;
                                     Krein-neutral markers for the zero sector and non-real eigenvalues (P-4)
  fig5_radiative_loss        Fig. 4  (a) Q(xi)/Q(0) with xi^-c; (b) tau_c(xi); (c) 1 - Q(40)/Q(0) against N delta3
  fig6_partner_entanglement  Fig. 5  E_N(omega) at the closed-form kappa = N S(k_op) (Eq. (7'), OD-P1 = A)
Same data, symbols and colours as the MATLAB script; captions.md and captions.tex explain each figure (explanation
only; the analysis belongs to the paper text).  No watermark: provenance and status are in figs_manifest.json;
figure_data.npz holds every plotted array.  Returns the manifest."""
import math
import os

import numpy as np

from . import jsonio, kit
from . import compat as C
from . import pubstyle as S
from .par import pmap
from .suite import SMALL_STAB, compute_all, phase2_params, slow_light


def _task(name):
    if name == "phi":
        T = kit.tod_soliton(1.0, 0.02, 500, 16.0)
        B = kit.bdg_spectrum_bg(T["phi"], 1.0, 0.02, 500, 16.0)
        return {"newton_res": T["res"][-1], "w": B["w"], "krein": B["krein"], "mu": B["mu"], "maxIm": B["maxIm"],
                "w_small": B["w_small"]}
    if name == "psi0":
        tau = -16.0 + np.arange(500) * (2.0 * 16.0 / 500)
        B = kit.bdg_spectrum_bg(1.0 * C.sech(1.0 * tau), 1.0, 0.02, 500, 16.0)
        return {"w_small": B["w_small"], "maxIm": B["maxIm"]}
    if name == "series":
        return kit.splitstep_loss(2.0, 0.05, 40.0, 0.002, 8192, 400.0, 0.0)
    raise ValueError(name)


def _entry(files, number, caption, data):
    print("wrote %s (%s)" % (", ".join(files), number))
    return {"file": files[0], "vector": files[1], "eps": files[2], "rendered": number, "caption": caption, "data": data}


IM_NEUTRAL = 1e-8          # |Im Omega| above this is a non-real eigenvalue (round-off of real ones is <= 1e-12 here)
N_ZERO_SECTOR = 4          # phase and translation modes with their generalized partners (two 2 x 2 Jordan blocks)


def krein_classes(w, krein):
    """Masks (positive, negative, neutral) over the whole spectrum w.  Krein-neutral: the zero-mode sector (the
    N_ZERO_SECTOR eigenvalues of smallest modulus) and every non-real eigenvalue.  With K = sigma_3 M Hermitian,
    M v = Omega v gives Im(Omega) (||u||^2 - ||v||^2) = 0, so the sign returned for these eigenvectors is round-off
    (finding A-11); for the zero sector ||u||^2 = ||v||^2 holds exactly for the phase mode (finding A-2)."""
    w = np.asarray(w)
    krein = np.asarray(krein)
    neutral = np.zeros(w.size, dtype=bool)
    neutral[np.argsort(np.abs(w), kind="stable")[:N_ZERO_SECTOR]] = True
    neutral |= np.abs(w.imag) > IM_NEUTRAL
    return (krein > 0) & ~neutral, (krein < 0) & ~neutral, neutral


def plot_fig3_bdg_krein(plt, w0, k0, mu0, wf, kf, muf, ws_psi0, ws_phi):
    """Fig. 3: (a) the spectrum about psi0 at delta3 = 0 (w0, krein k0, mu0); (b) about the stationary TOD soliton phi
    at delta3 = 0.02 (wf, kf, muf); inset of (b): the four eigenvalues of smallest modulus about psi0 (ws_psi0) and
    about phi (ws_phi) at delta3 = 0.02.  Returns (figure, info)."""
    panels = []
    for w, kr, mu in ((w0, k0, mu0), (wf, kf, muf)):
        w, kr = np.asarray(w), np.asarray(kr)
        pos, neg, neu = krein_classes(w, kr)
        win = np.abs(w.real) <= 3 * mu
        panels.append((w, pos & win, neg & win, neu & win, mu, win))
    mx = max(float(np.max(np.abs(w[win].imag))) for w, _, _, _, _, win in panels)
    ym = 1.15 * mx                          # every eigenvalue of the window inside the frame (A-10)
    f = S.figure(plt, S.W2, 5.8)
    axes = f.subplots(1, 2)
    outside, nneu = [], []
    for p, (a, (w, pos, neg, neu, mu, win), tag) in enumerate(zip(axes, panels, ("(a)", "(b)"))):
        outside.append(int(np.sum(win & (np.abs(w.imag) > ym))))
        nneu.append(int(np.sum(neu)))
        hb = a.fill([-mu, -mu, mu, mu], [-ym, ym, ym, -ym], color=S.BAND, ec="none", zorder=0)[0]
        hp, = a.plot(w[pos].real, w[pos].imag, "o", mec=S.KPOS, mfc="none", ls="none", zorder=3)
        hn, = a.plot(w[neg].real, w[neg].imag, "D", ms=S.MS - 0.5, mec=S.KNEG, mfc="none", ls="none", zorder=3)
        h0, = a.plot(w[neu].real, w[neu].imag, "x", ms=S.MS + 0.5, mec="k", mew=0.8, ls="none", zorder=4)
        a.set_xlim(-3 * mu, 3 * mu)
        a.set_ylim(-ym, ym)
        a.set_xlabel(r"$\mathrm{Re}\,\Omega$")
        S.frame(a)
        S.scaled_axis(a, "y", r"$\mathrm{Im}\,\Omega$")
        S.label(a, tag)
        if p == 0:
            a.legend([hp, hn, h0, hb], [r"Krein $+$", r"Krein $-$", r"Krein-neutral", r"$|\mathrm{Re}\,\Omega|<\mu$"],
                     loc="lower left")
    ai = axes[1].inset_axes([0.66, 0.54, 0.31, 0.42])
    ai.set_xlim(-0.16, 0.16)              # limits and ticks before the data: a small inset cannot autoscale with
    ai.set_ylim(-0.16, 0.16)              # round_numbers limits (Matplotlib's MaxNLocator gets no tick space)
    ai.set_xticks([-0.1, 0.0, 0.1])
    ai.set_yticks([-0.1, 0.0, 0.1])
    ai.plot(np.asarray(ws_psi0).real, np.asarray(ws_psi0).imag, "s", ms=3.0, mec=S.GREY, mfc="none", ls="none")
    ai.plot(np.asarray(ws_phi).real, np.asarray(ws_phi).imag, "o", ms=3.0, mec="k", mfc="k", ls="none")
    ai.tick_params(which="both", direction="in", top=True, right=True, length=2.0, labelbottom=False, labelleft=False)
    ai.set_facecolor("w")
    for sp in ai.spines.values():
        sp.set_linewidth(S.LW_AX)
    info = {"ylim": ym, "n_outside_ylim": outside, "n_neutral": nneu,
            "max_abs_im_in_window": [float(np.max(np.abs(w[win].imag))) for w, _, _, _, _, win in panels]}
    return f, info


def fig3_caption(N, n_tau, tau_max, info):
    """Caption of Fig. 3 (explanation only)."""
    grid = r"a periodic grid of $n_\tau=%d$ points, $\tau\in[-%g,%g)$, with fourth-order finite differences" % (
        n_tau, tau_max, tau_max)
    cap = (r"Bogoliubov–de Gennes eigenvalues $\Omega$ for $N=%g$ on %s, in the window $|\mathrm{Re}\,\Omega|\le3\mu$, "
           r"$\mu=N^2/2$. (a) $\delta_3=0$, linearization about $\psi_0=N\,\mathrm{sech}(N\tau)$. (b) $\delta_3=0.02$, "
           r"linearization about the stationary third-order-dispersion soliton $\phi$. Open circles: positive Krein signature; "
           r"open diamonds: negative Krein signature; crosses: Krein-neutral eigenvalues (the zero-mode sector and every "
           r"non-real eigenvalue); shaded band: $|\mathrm{Re}\,\Omega|<\mu$. Inset of (b), "
           r"$|\mathrm{Re}\,\Omega|,|\mathrm{Im}\,\Omega|\le0.16$ with ticks at 0 and $\pm0.1$: the four eigenvalues of smallest "
           r"modulus about $\phi$ (filled circles) and about $\psi_0$ at $\delta_3=0.02$ (open squares)." % (N, grid))
    for tag, n in zip(("(a)", "(b)"), info["n_outside_ylim"]):
        if n:
            cap += r" %d eigenvalue%s of %s lie%s outside the vertical range shown." % (n, "" if n == 1 else "s", tag,
                                                                                    "s" if n == 1 else "")
    return cap


def fig3_from_data(npz_path, outdir, N=1.0, n_tau=500, tau_max=16.0):
    """Re-render Fig. 3 from the figure_data.npz of a run (no computation).  Returns (files, caption, info)."""
    Z = np.load(npz_path, allow_pickle=False)
    plt = S.setup()
    os.makedirs(outdir, exist_ok=True)
    mu = float(Z["fig2_mu"])
    f, info = plot_fig3_bdg_krein(plt, Z["fig2_w_psi0_d3_0"], Z["fig2_krein_psi0_d3_0"], mu, Z["fig2_w_phi_d3_0p02"],
                                  Z["fig2_krein_phi_d3_0p02"], mu, Z["fig2_w_small_psi0"], Z["fig2_w_small_phi"])
    files = S.save(plt, f, outdir, "fig2_bdg_krein")
    return files, fig3_caption(N, n_tau, tau_max, info), info


def make_paper2_figs(outdir, pack_root, g4c_report, prov):
    plt = S.setup()
    os.makedirs(outdir, exist_ok=True)
    P = phase2_params(SMALL_STAB)          # the stability map feeds no figure here: 2 x 2 at n_tau = 48
    R = compute_all(P, os.path.join(pack_root, "data", "entropy_trajectory.csv"))
    have_rep = bool(g4c_report) and os.path.exists(g4c_report)
    names = ["series", "phi", "psi0"] if have_rep else ["phi", "psi0"]
    T = dict(zip(names, pmap(_task, names)))
    from . import PACKAGE, __version__
    M = {"version": "0.2.1 (Python port, %s %s; publication style after IOP guidelines)" % (PACKAGE, __version__),
         "provenance": prov, "status": "RECORDED (HR-3: only MATLAB R2025b output is canonical)",
         "style": {"guidelines": "IOP Publishing", "widths_cm": [8.5, 15.0], "font": S.font_choice(),
                   "files": ["pdf", "eps", "png 600 dpi"]}, "figures": []}
    FD = {}                                              # every plotted array, saved as figure_data.npz
    op = r"$N=%g$, $\delta_3=%g$, $\kappa_g=%g$, $k_{\rm op}=%g$" % (P["N"], P["d3"], P["kg"], P["k_op"])

    # ---------------- Fig. 1: fig1_thermality
    f = S.figure(plt, S.W2, 5.4)
    a, b = f.subplots(1, 2)
    a.plot(R["omega"], np.log(R["ratio"]), "o", ms=3.0, mec=S.ORANGE, mfc=S.ORANGE, ls="none", clip_on=False, zorder=3)
    xl = np.array([np.min(R["omega"]), np.max(R["omega"])])
    a.plot(xl, R["slope"] * xl + R["intercept"], "k-", lw=0.8, zorder=4)   # the line over the 60 markers
    a.set_xlim(0, 1.65)
    a.set_xlabel(r"$\omega$")
    a.set_ylabel(r"$\ln|\beta_\omega/\alpha_\omega|^2$")
    S.label(a, "(a)")
    S.frame(a)
    b.plot(R["drive"], R["kappa_drive"], "s", mec=S.BLUE, mfc=S.BLUE, ls="none", clip_on=False, zorder=3)
    m = float(np.mean(R["kappa_drive"]))
    b.set_ylim(m - 0.06, m + 0.06)
    b.set_xlim(0.65, 1.35)
    b.set_xlabel(r"drive amplitude $d$")
    b.set_ylabel(r"$\kappa$")
    S.label(b, "(b)")
    S.frame(b)
    FD.update(fig1_omega=R["omega"], fig1_ln_ratio=np.log(R["ratio"]), fig1_slope=R["slope"], fig1_intercept=R["intercept"],
              fig1_drive=R["drive"], fig1_kappa_drive=R["kappa_drive"])
    cap1 = (r"(a) $\ln|\beta_\omega/\alpha_\omega|^2$ for $|\beta_\omega/\alpha_\omega|^2=e^{-2\pi\omega/\kappa}$ at the kinematic "
            r"surface gravity $\kappa=%.4f$, at %d frequencies $\omega\in[%.2f,%.2f]$ (circles), and its unweighted "
            r"least-squares straight line (solid). (b) Kinematic surface gravity $\kappa$ at %d drive amplitudes "
            r"$d\in[%.1f,%.1f]$. Parameters: %s."
            % (R["kappa_kin"], len(R["omega"]), float(R["omega"][0]), float(R["omega"][-1]), len(R["drive"]),
               float(R["drive"][0]), float(R["drive"][-1]), op))
    M["figures"].append(_entry(S.save(plt, f, outdir, "fig1_thermality"), "Fig. 1", cap1,
                               {"kappa_extracted": R["kappa_kin"], "slope": R["slope"], "R2": R["R2"],
                                "kappa_drive_spread": R["kappa_drive_spread"]}))

    # ---------------- Fig. 2: fig3_gapless_flow
    f = S.figure(plt, S.W1, 6.2)
    a = f.add_subplot(1, 1, 1)
    a.plot(R["kg_flow"], R["kappa_flow"], "k-", zorder=2)
    a.plot(R["kg_flow"], R["kappa_flow"], "o", mec=S.ORANGE, mfc=S.ORANGE, ls="none", clip_on=False, zorder=3)
    a.plot([0, 1], [P["N"], P["N"]], "k:", lw=0.8)
    a.plot([0], [P["N"]], "s", ms=5.5, mec=S.RED, mfc=S.RED, ls="none", clip_on=False, zorder=4)
    a.set_xlim(-0.03, 1.03)
    yl = (0.9 * float(np.min(R["kappa_flow"])), 1.04 * P["N"])
    a.set_ylim(*yl)
    a.set_xlabel(r"$\kappa_g$")
    a.set_ylabel(r"$\kappa$")
    S.frame(a, mirror=False)
    a.tick_params(which="both", top=True)
    a2 = a.twinx()
    a2.set_ylim(yl[0] / (2 * math.pi), yl[1] / (2 * math.pi))
    a2.set_ylabel(r"$T_H=\kappa/2\pi$ (local)")
    a2.tick_params(which="both", direction="in", right=True, labelsize=S.FS_SMALL)
    S.nice_ticks(a2)
    for sp in a2.spines.values():
        sp.set_linewidth(S.LW_AX)
    FD.update(fig3_kg=R["kg_flow"], fig3_kappa=R["kappa_flow"], fig3_N=P["N"])
    cap2 = (r"Kinematic surface gravity $\kappa$ against the slow-light coupling gap $\kappa_g$ at %d values "
            r"$\kappa_g\in[%g,%g]$ (circles joined by a line), for $N=%g$, $\delta_3=%g$, $k_{\rm op}=%g$, $d=%g$. Square: "
            r"$\kappa=N$ at $\kappa_g=0$; dotted line: the level $\kappa=N$. Right-hand axis: the local temperature "
            r"$T_H=\kappa/2\pi$."
            % (len(R["kg_flow"]), float(R["kg_flow"][0]), float(R["kg_flow"][-1]), P["N"], P["d3"], P["k_op"], P["drive0"]))
    M["figures"].append(_entry(S.save(plt, f, outdir, "fig3_gapless_flow"), "Fig. 2", cap2,
                               {"kappa_kg0": float(R["kappa_flow"][0]), "N": P["N"], "n_kg": int(R["kg_flow"].size)}))

    # ---------------- Fig. 3: fig2_bdg_krein (T-13 in panel b; P-4: Krein-neutral markers, full vertical range)
    o0 = [x["out"] for x in R["bdg"] if x["d3"] == 0.0][0]
    Bf, Bp = T["phi"], T["psi0"]
    f, info = plot_fig3_bdg_krein(plt, o0["w"], o0["krein"], o0["mu"], Bf["w"], Bf["krein"], Bf["mu"], Bp["w_small"],
                                  Bf["w_small"])
    FD.update(fig2_w_psi0_d3_0=o0["w"], fig2_krein_psi0_d3_0=o0["krein"], fig2_w_phi_d3_0p02=Bf["w"],
              fig2_krein_phi_d3_0p02=Bf["krein"], fig2_w_small_phi=Bf["w_small"], fig2_w_small_psi0=Bp["w_small"], fig2_mu=o0["mu"],
              fig2_ylim=info["ylim"])
    cap3 = fig3_caption(P["bdg_N"], P["bdg_n_tau"], P["bdg_tau_max"], info)
    M["figures"].append(_entry(S.save(plt, f, outdir, "fig2_bdg_krein"), "Fig. 3", cap3, {
        "maxIm_d3_0": o0["maxIm"], "maxIm_about_phi": Bf["maxIm"], "zero_sector_phi": float(np.max(np.abs(Bf["w_small"]))),
        "zero_sector_psi0": float(np.max(np.abs(Bp["w_small"]))), "newton_res": Bf["newton_res"],
        "ylim": info["ylim"], "n_outside_ylim": info["n_outside_ylim"], "n_krein_neutral": info["n_neutral"],
        "max_abs_im_in_window": info["max_abs_im_in_window"]}))

    # ---------------- Fig. 4: fig5_radiative_loss
    if have_rep:
        G = jsonio.read(g4c_report)
        pts = G["points"]
        opt = [q for q in pts if q["N"] == 2 and abs(q["d3"] - 0.05) < 1e-12][0]
        Y = T["series"]
        f = S.figure(plt, S.W1, 13.0)
        a, b, c = f.subplots(3, 1)
        s = Y["xi"] >= 0.5
        qq = Y["Q"][s] / Y["Q"][0]
        a.loglog(Y["xi"][s], qq, "k-")
        i10 = int(np.argmin(np.abs(Y["xi"] - 10)))
        xx = C.mlinspace(10, 40, 50)
        a.loglog(xx, (Y["Q"][i10] / Y["Q"][0]) * (xx / 10) ** (-opt["c_over_xi"]), "--", color=S.RED)
        a.set_xlim(0.5, 40)
        ylo = math.floor(float(np.min(qq)) * 200) / 200
        a.set_ylim(ylo, 1.002)
        yt = np.arange(math.ceil(ylo * 100), 101) / 100
        xt = [0.5, 1, 2, 5, 10, 20, 40]
        from matplotlib.ticker import FixedFormatter, FixedLocator, NullFormatter, NullLocator
        a.yaxis.set_major_locator(FixedLocator(yt))
        a.yaxis.set_major_formatter(FixedFormatter(["%.2f" % v for v in yt]))
        a.yaxis.set_minor_locator(NullLocator())
        a.yaxis.set_minor_formatter(NullFormatter())
        a.xaxis.set_major_locator(FixedLocator(xt))
        a.xaxis.set_major_formatter(FixedFormatter(["%g" % v for v in xt]))
        a.xaxis.set_minor_formatter(NullFormatter())
        a.set_xlabel(r"$\xi$")
        a.set_ylabel(r"$Q(\xi)/Q(0)$")
        S.label(a, "(a)")
        S.frame(a)
        b.plot(Y["xi"], Y["tau_c"], "k-")
        b.set_xlim(0, 40)
        b.set_xlabel(r"$\xi$")
        b.set_ylabel(r"$\tau_c$")
        S.label(b, "(b)")
        S.frame(b)
        n1 = [q for q in pts if q["N"] == 1]
        n2 = [q for q in pts if q["N"] == 2]
        # N = 2 (filled squares) under N = 1 (larger open circles): the pair that coincides at N delta3 = 0.05 shows both
        h2, = c.semilogy([q["Nd3"] for q in n2], [q["loss_0_40"] for q in n2], "s", ms=3.5, mec="k", mfc="k", ls="none",
                         clip_on=False, zorder=3)
        h1, = c.semilogy([q["Nd3"] for q in n1], [q["loss_0_40"] for q in n1], "o", ms=5.5, mec="k", mfc="none", mew=0.7,
                         ls="none", clip_on=False, zorder=4)
        c.set_xlim(0, 0.11)
        c.set_xlabel(r"$N\delta_3$")
        c.set_ylabel(r"$1-Q(40)/Q(0)$")
        S.label(c, "(c)")
        S.frame(c)
        c.legend([h1, h2], [r"$N=1$", r"$N=2$"], loc="lower right")
        FD.update(fig5_xi=Y["xi"], fig5_Q=Y["Q"], fig5_tau_c=Y["tau_c"], fig5_tau_c_grid=Y["tau_c_grid"], fig5_k_s=Y["k_s"],
                  fig5_c=opt["c_over_xi"], fig5_points_N=np.array([q["N"] for q in pts]),
                  fig5_points_Nd3=np.array([q["Nd3"] for q in pts]), fig5_points_loss_0_40=np.array([q["loss_0_40"] for q in pts]))
        cap4 = (r"Split-step Fourier integration of the generalized nonlinear Schrödinger equation with third-order dispersion, "
                r"launched from $\psi_0=N\,\mathrm{sech}(N\tau)$, on $\tau\in[-L/2,L/2)$ with $L=400$ and $n=8192$ points, step "
                r"$\Delta\xi=0.002$, and an absorber at $|\tau|>0.35L$. $Q$: the norm within $6/N$ of the intensity maximum; "
                r"$\tau_c$: the position of the intensity maximum (vertex of the parabola through $\ln|\psi|^2$ at the grid "
                r"maximum and its two neighbours). (a) $Q(\xi)/Q(0)$ at $N=2$, $\delta_3=0.05$ (solid), and the power law "
                r"$[Q(10)/Q(0)]\,(\xi/10)^{-c}$ (dashed), where $c=%s$ is the mean, over the windows $[10,20]$, $[20,30]$ and "
                r"$[30,40]$, of the window midpoint times the loss rate $-\mathrm{d}\ln Q/\mathrm{d}\xi$ fitted in the window. "
                r"(b) $\tau_c(\xi)$ of the same propagation. (c) $1-Q(40)/Q(0)$ against $N\delta_3$ for $N=1$ (open circles) "
                r"and $N=2$ (filled squares), from propagations with the same settings." % kit.fmt_sci(opt["c_over_xi"], 1))
        M["figures"].append(_entry(S.save(plt, f, outdir, "fig5_radiative_loss"), "Fig. 4", cap4, {
            "c": opt["c_over_xi"], "loss_0_40_op": opt["loss_0_40"], "series_loss_rate": Y["loss_rate"],
            "series_matches_report": abs(Y["loss_rate"] - opt["loss_rate"]) <= 1e-12 * max(1, abs(opt["loss_rate"])),
            "g4c_report": g4c_report}))
    else:
        print("fig5_radiative_loss skipped: no g4c_loss_report.json (run stage S2 first or set P2_G4C_REPORT)")

    # ---------------- Fig. 5: fig6_partner_entanglement
    kcf = P["N"] * float(slow_light(P["k_op"], P["d3"], P["kg"]))
    Gl = kit.gauss_layer(kcf)
    f = S.figure(plt, S.W1, 5.6)
    a = f.add_subplot(1, 1, 1)
    a.fill(np.concatenate([Gl["w"], Gl["w"][::-1]]), np.concatenate([Gl["EN"], np.zeros(Gl["w"].size)]),
           color=S.FILL, ec="none", zorder=0)
    a.plot(Gl["w"], Gl["EN"], "-", color=S.ORANGE, lw=1.4)
    a.set_xlim(0, 1.65)
    a.set_ylim(0, 1.08 * float(np.max(Gl["EN"])))
    a.set_xlabel(r"$\omega$")
    a.set_ylabel(r"$E_N$ (nats)")
    S.frame(a)
    FD.update(fig6_omega=Gl["w"], fig6_E_N=Gl["EN"], fig6_nu_minus=Gl["nu"], fig6_nbar_max=Gl["nbar_max_w"], fig6_kappa=kcf)
    cap5 = (r"Logarithmic negativity of the Hawking partners, $E_N(\omega)=2\,\mathrm{arcsinh}|\beta_\omega|$ with "
            r"$|\beta_\omega|^2=(e^{2\pi\omega/\kappa}-1)^{-1}$, at %d frequencies $\omega\in[%.2f,%.2f]$, for the closed-form "
            r"surface gravity $\kappa=NS(k_{\rm op})=%.4f$ (%s)." % (Gl["w"].size, float(Gl["w"][0]), float(Gl["w"][-1]), kcf, op))
    M["figures"].append(_entry(S.save(plt, f, outdir, "fig6_partner_entanglement"), "Fig. 5", cap5,
                               {"kappa_closed_form": kcf, "EN_max": Gl["EN_max"], "tokens": Gl["tokens"]}))
    np.savez_compressed(os.path.join(outdir, "figure_data.npz"), **FD)
    order = {"Fig. 1": 1, "Fig. 2": 2, "Fig. 3": 3, "Fig. 4": 4, "Fig. 5": 5}
    caps = sorted(((e["file"][:-4], "Figure %d" % order[e["rendered"]], e["caption"]) for e in M["figures"]), key=lambda x: x[1])
    M["captions"] = S.write_captions(outdir, caps, "Paper 2, manuscript figures: captions")
    jsonio.write(os.path.join(outdir, "figs_manifest.json"), M)
    return M
