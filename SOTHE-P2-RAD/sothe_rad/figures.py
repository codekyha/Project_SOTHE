"""The article figure of SOTHE-P2-RAD (steady emission and recoil) and its plotted data, from an analysed results folder.

  fig_emission_recoil.{pdf,eps,png}  (a) scaled rate Rs = R/N^2 against 1/delta_tilde: steady prepared solitons (mean
                                     over the hold phase), every hold-phase rate of the prepared solitons (the calibration
                                     set), the rates of the launched solitons after the launch transient (the test set), the
                                     rate function and the first-order rate;
                                     (b) norm lost by three launched solitons, direct integration and recoil model;
                                     (c) their mean wavenumber k_s; (d) the exponent c of -d ln Q/dxi ~ c/xi over the windows
                                     [10,20], [20,40], [40,80], [80,160] (at the geometric centre of each window), run and model.
  figure_data.npz                    every plotted array
Style: pubstyle.py (the style of the other figures of the article)."""
import csv
import json
import math
import os

import numpy as np

from . import model, pubstyle as S

FIG_RUNS = ["recoil_N1_d0p080", "recoil_N1_d0p100", "recoil_N1_d0p120"]


def _csv(path):
    with open(path, encoding="utf-8") as f:
        return list(csv.DictReader(f))


def make_all(folder):
    an = os.path.join(folder, "analysis")
    with open(os.path.join(an, "summary.json"), encoding="utf-8") as f:
        Sm = json.load(f)
    out = os.path.join(folder, "figures")
    os.makedirs(out, exist_ok=True)
    plt = S.setup()
    F = Sm["rate_function"]
    C, p, g = F["C"], F["p"], F["g"]
    up = _csv(os.path.join(an, "universal_points.csv"))
    wex = _csv(os.path.join(an, "window_exponents.csv"))
    T = np.load(os.path.join(an, "model_tracks.npz"))
    FD = {}

    fig = S.figure(plt, S.W2, 11.0)
    axs = fig.subplots(2, 2)
    a, b, c, d = axs[0, 0], axs[0, 1], axs[1, 0], axs[1, 1]
    # (a) scaled rate against 1/delta_tilde
    lo, hi = F["delta_tilde_range"]
    dts = np.linspace(lo, hi, 300)
    fit = np.array([C * x ** p * math.exp(g * x - math.pi * model.q_res(x)) for x in dts])
    dtb = np.linspace(0.03, 0.105, 300)
    born = np.array([model.rate_born_scaled(x) for x in dtb])
    tst = [u for u in up if u["role"] == "test"]
    cal = [u for u in up if u["role"] == "calibration"]
    xs_t, ys_t = np.array([1 / float(u["delta_tilde"]) for u in tst]), np.array([float(u["Rs"]) for u in tst])
    xs_c, ys_c = np.array([1 / float(u["delta_tilde"]) for u in cal]), np.array([float(u["Rs"]) for u in cal])
    st = [x for x in Sm["steady"] if x["steady"]]
    xs_s = np.array([1 / x["delta_tilde"] for x in st])
    ys_s = np.array([x["R_flux"] / x["N_local"] ** 2 for x in st])
    a.semilogy(xs_t, ys_t, '.', ms=1.2, color=S.GREY, ls='none', zorder=1)
    a.semilogy(xs_c, ys_c, '.', ms=1.2, color=S.ORANGE, ls='none', zorder=2)
    a.semilogy(1 / dtb, born, '--', color=S.BLUE, zorder=3, dashes=(4, 2))
    a.semilogy(1 / dts, fit, '-', color='k', zorder=4)
    a.semilogy(xs_s, ys_s, 'o', ms=4.0, mfc='none', mec=S.ORANGE, mew=0.9, ls='none', zorder=5, clip_on=False)
    a.set_xlim(9.5, 30.0)
    a.set_ylim(1e-17, 1e-3)
    a.set_xlabel(r'$1/\tilde\delta$')
    a.set_ylabel(r'$R/N^2$')
    S.label(a, '(a)')
    S.frame(a)
    FD.update(a_dt_fit=dts, a_fit=fit, a_dt_born=dtb, a_born=born, a_test_inv_dt=xs_t, a_test_Rs=ys_t, a_calibration_inv_dt=xs_c,
              a_calibration_Rs=ys_c, a_steady_inv_dt=xs_s, a_steady_Rs=ys_s)
    # (b)-(d) launched runs against the model
    styles = [('-', S.BLUE, 's'), ('-', 'k', 'o'), ('-', S.RED, '^')]
    for rid, (ls, col, mk) in zip(FIG_RUNS, styles):
        key = lambda part, q: "%s__%s__%s" % (rid, part, q)
        if key("xi", "v") not in T.files:
            continue
        xi, Q, ks = T[key("xi", "v")], T[key("Q", "v")], T[key("ks", "v")]
        Ql = float(T[key("Qlaunch", "v")])
        mx, mQ, mks = T[key("A", "xi")], T[key("A", "Q")], T[key("A", "ks")]
        sel = xi >= 1.0
        b.plot(xi[sel], 100 * (1 - Q[sel] / Ql), ls, color=col, zorder=3)
        b.plot(mx, 100 * (1 - mQ / Ql), '--', color=col, lw=1.4, dashes=(3, 2), zorder=4)
        c.plot(xi, ks, ls, color=col, zorder=3)
        c.plot(mx, mks, '--', color=col, lw=1.4, dashes=(3, 2), zorder=4)
        w = [x for x in wex if x["id"] == rid]
        mids = np.array([math.sqrt(float(x["window"].strip("[]").split(",")[0]) * float(x["window"].strip("[]").split(",")[1])) for x in w])
        cg, cm = np.array([float(x["c_gnlse"]) for x in w]), np.array([float(x["c_model"]) for x in w])
        d.plot(mids, 100 * cg, mk, ms=4.0, mfc='none', mec=col, mew=0.9, ls='none', zorder=4)
        d.plot(mids, 100 * cm, '--', color=col, lw=1.4, dashes=(3, 2), zorder=3)
        FD.update({rid + "_xi": xi, rid + "_Q": Q, rid + "_ks": ks, rid + "_model_xi": mx, rid + "_model_Q": mQ,
                   rid + "_model_ks": mks, rid + "_Qlaunch": Ql, rid + "_window_mid": mids, rid + "_c_gnlse": cg, rid + "_c_model": cm})
    for ax, lab, yl in ((b, '(b)', r'$100\,[1-Q(\xi)/Q(0)]$'), (c, '(c)', r'$k_s$')):
        ax.set_xlim(0, 160)
        ax.set_xlabel(r'$\xi$')
        ax.set_ylabel(yl, rotation=0 if lab == '(c)' else 90, labelpad=8 if lab == '(c)' else 4)
        S.label(ax, lab)
        S.frame(ax)
    c.set_ylim(-0.3, 0.02)
    d.set_xscale('log')
    d.set_xlim(10, 160)
    d.set_xticks([10, 20, 40, 80, 160])
    d.set_xticklabels(['10', '20', '40', '80', '160'])
    d.minorticks_off()
    d.set_xlabel(r'$\xi$')
    d.set_ylabel(r'$100\,c$')
    S.label(d, '(d)')
    S.frame(d)
    S.save(plt, fig, out, "fig_emission_recoil")
    np.savez_compressed(os.path.join(out, "figure_data.npz"), **{k: np.asarray(v) for k, v in FD.items()})
    print("figures written to %s" % out)
