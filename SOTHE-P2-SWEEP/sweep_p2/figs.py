"""Figures of the sweep, in the publication style of SOTHE-P2 (sothe_p2/pubstyle.py: IOP sizes, Computer Modern,
no titles or annotations, PDF + EPS + PNG at 600 dpi, captions in captions.md / captions.tex).  Diagnostic figures
of RECORDED numbers (HR-3), not manuscript figures."""
import math
import os

import numpy as np

from sothe_p2 import pubstyle as ps

from . import b_census
from .tasks import _f


def _panel_pair(plt, h_cm=6.2):
    fig = ps.figure(plt, ps.W2, h_cm)
    ax = fig.subplots(1, 2)
    for a, lab in zip(ax, ("(a)", "(b)")):
        ps.frame(a)
        ps.label(a, lab)
    return fig, ax


def fig_window_map(plt, outdir, cfg, S):
    tab = (S.get("B4") or {}).get("window_table") or []
    if not tab:
        return None
    Ns = sorted(set(t["N"] for t in tab))
    Ks = sorted(set(t["k_op"] for t in tab))
    if len(Ns) < 2 or len(Ks) < 2:
        return None
    fig, ax = _panel_pair(plt)
    for a, key in zip(ax, ("lab_frac", "late_frac")):
        Z = np.full((len(Ns), len(Ks)), np.nan)
        for t in tab:
            v = t.get(key)
            Z[Ns.index(t["N"]), Ks.index(t["k_op"])] = np.nan if v is None else 100.0 * v
        def edges(v):
            v = np.asarray(v, dtype=float)
            m = 0.5 * (v[1:] + v[:-1])
            return np.concatenate([[v[0] - (m[0] - v[0])], m, [v[-1] + (v[-1] - m[-1])]])
        pc = a.pcolormesh(edges(Ks), edges(Ns), Z, cmap="Greys", vmin=0, vmax=100, shading="flat", edgecolors="w", linewidth=0.3)
        a.axhline(math.sqrt(3.2), color=ps.RED, lw=ps.LW, ls="--")
        a.plot([cfg["op"]["k_op"]], [cfg["op"]["N"]], ls="none", marker="x", color=ps.ORANGE, ms=6, mew=1.2)
        a.plot([cfg["anchor"]["k_op"]], [cfg["anchor"]["N"]], ls="none", marker="o", mfc="none", mec=ps.BLUE, ms=6, mew=1.2)
        a.set_xlabel(r"$k_{\rm op}$")
        a.set_ylabel(r"$N$")
    cb = fig.colorbar(pc, ax=ax, shrink=0.9, pad=0.02)
    cb.set_label(r"Hawking window (% of band)")
    cb.outline.set_linewidth(ps.LW_AX)
    return ps.save(plt, fig, outdir, "sweep_fig1_hawking_window"), (
        "sweep_fig1_hawking_window", "Figure S1",
        r"Hawking window $[\max(\mu,0.05),\min(\omega^{(-)}_{\max},1.6)]$ as a percentage of the band $[0.05,1.6]$ over $(N,k_{\rm op})$ at "
        r"$\delta_3=0.05$, $\kappa_g=0.3$: (a) laboratory frame, (b) frame comoving with the drift measured over $\xi\in[30,40]$. "
        r"Dashed: $N=\sqrt{2\cdot1.6}$, above which $\mu=N^2/2$ lies outside the band. Cross: operating point $(2,1)$; circle: anchor $(1,1.5)$.")


def fig_raman(plt, outdir, cfg, S):
    runs = S.get("_B6_runs") or []
    tab = (S.get("B6") or {}).get("S1_3_table") or []
    if not runs:
        return None
    fig, ax = _panel_pair(plt)
    N0 = 2.0 if any(abs(r["N"] - 2.0) < 1e-9 for r in runs) else runs[0]["N"]
    T0s = sorted(set(r["T0_fs"] for r in runs if r["model"] != "none"))
    T0 = 50.0 if 50.0 in T0s else (T0s[0] if T0s else None)
    styles = {"none": ("k", "-"), "ss": (ps.GREY, ":"), "lin3": (ps.RED, "--"), "bw": (ps.BLUE, "-."), "bw_ss": (ps.ORANGE, "-")}
    names = {"none": "Kerr + TOD", "ss": "+ self-steepening", "lin3": r"+ Raman, $T_R=3$ fs", "bw": "+ Raman (Blow-Wood)",
             "bw_ss": "+ Raman (BW) + s.-s."}
    for r in runs:
        if abs(r["N"] - N0) > 1e-9 or (r["model"] != "none" and T0 is not None and abs(r["T0_fs"] - T0) > 1e-9):
            continue
        c, ls = styles.get(r["model"], ("k", "-"))
        ax[0].plot(r["series"]["xi"], r["series"]["k_s"], color=c, ls=ls, label=names.get(r["model"], r["model"]))
    ax[0].set_xlabel(r"$\xi$")
    ax[0].set_ylabel(r"soliton centroid $k_s$")
    ax[0].legend(loc="lower left")
    sel = [t for t in tab if abs(t["N"] - N0) < 1e-9]
    for model, mk in (("lin3", "s"), ("bw", "o"), ("bw_ss", "^")):
        pts = sorted([t for t in sel if t["model"] == model], key=lambda t: t["T0_fs"])
        if pts:
            c, _ = styles[model]
            ax[1].loglog([t["T0_fs"] for t in pts], [abs(t["initial_rate_kfull"]) for t in pts], ls="none", marker=mk, color=c, mfc="none",
                         label=names[model])
            ax[1].loglog([t["T0_fs"] for t in pts], [abs(t["early_rate_ks"]) for t in pts], ls="none", marker=mk, color=c, ms=2.5)
            ax[1].loglog([t["T0_fs"] for t in pts], [abs(t["gordon_rate"]) for t in pts], color=c, ls=":", lw=0.8)
    rec = [t["recoil_dVd_dxi"] for t in sel if t.get("recoil_dVd_dxi")]
    if rec:
        ax[1].axhline(abs(rec[0]), color="k", ls="--", lw=ps.LW)
    from matplotlib.ticker import FixedLocator, NullFormatter, ScalarFormatter
    if T0s:
        ax[1].set_xlim(min(T0s) / 1.6, max(T0s) * 1.6)
        ax[1].xaxis.set_major_locator(FixedLocator(T0s))
        ax[1].xaxis.set_major_formatter(ScalarFormatter())
        ax[1].xaxis.set_minor_formatter(NullFormatter())
    ax[1].set_xlabel(r"$T_0$ (fs)")
    ax[1].set_ylabel(r"$|dk/d\xi|$")
    ax[1].legend(loc="upper right")
    return ps.save(plt, fig, outdir, "sweep_fig2_raman"), (
        "sweep_fig2_raman", "Figure S2",
        r"Eq. (1) with intrapulse Raman scattering and self-steepening at $N=%g$, $\delta_3=0.05$ (RK4IP). (a) Core spectral centroid $k_s(\xi)$ at "
        r"$T_0=%g$ fs for the Kerr+TOD model and its extensions. (b) Self-frequency-shift rate against $T_0$: initial rate of the full-field centroid "
        r"over $\xi\in[0,0.5]$ (open markers) and mean rate of the core centroid $k_s$ over $\xi\in[0,5]$ (small filled markers), with "
        r"Gordon's rate $\tfrac{8}{15}\tau_RN^4$ (dotted); dashed: the recoil-driven change of the drift velocity, $|dV_d/d\xi|$, of the Kerr+TOD run."
        % (N0, T0 if T0 else float("nan")))


VARIANT_LABEL = {"V1_shipped": "V1", "V2_prod": "V2", "V3_sg8_shipped": "V3a", "V3_sg8_prod": "V3b", "V4_d3draft_prod": "V4",
                 "V5_xi12_prod": "V5", "V1_shipped_fast": "V1 (fast)"}


def fig_ref25(plt, outdir, cfg, S):
    vs = (S.get("B8") or {}).get("variants") or []
    if not vs:
        return None
    fig, ax = _panel_pair(plt)
    x = np.arange(len(vs))
    for a, key, band, lab in ((ax[0], "DeltaS_range", (0.82, 2.42), r"$\Delta S_{\rm tot}$ (nats)"),
                              (ax[1], "eta_range", (0.52, 0.62), r"$\eta_{\rm GSL}$")):
        a.axhspan(band[0], band[1], color=ps.BAND, lw=0)
        for i, v in enumerate(vs):
            lo, hi = v[key]
            a.plot([i, i], [lo, hi], color="k", lw=1.6, solid_capstyle="butt")
            a.plot([i], [lo], marker="v", color="k", ms=4)
            a.plot([i], [hi], marker="^", color="k", ms=4)
        if key == "eta_range":
            a.axhline(0.5, color=ps.RED, ls="--", lw=ps.LW)
        else:
            a.axhline(0.0, color=ps.RED, ls="--", lw=ps.LW)
        a.set_xticks(x)
        a.set_xticklabels([VARIANT_LABEL.get(v["variant"], v["variant"].replace("_", " ")) for v in vs])
        a.set_xlim(-0.6, len(vs) - 0.4)
        a.set_ylabel(lab)
    return ps.save(plt, fig, outdir, "sweep_fig3_ref25_variants"), (
        "sweep_fig3_ref25_variants", "Figure S3",
        r"Ranges of (a) $\Delta S_{\rm tot}$ and (b) $\eta_{\rm GSL}$ over the configurations of each variant of the Ref. [25] robustness sweep "
        r"(V1: script as shipped, $N_t=2^{13}$, 4000 steps; V2: production resolution $N_t=2^{14}$, 6000 steps; V3a, V3b: super-Gaussian "
        r"$Ne^{-\tau^8/2}$ of the version of record at the two resolutions; V4: draft $\delta_3$ grid; V5: $\xi_{\max}=12$). "
        r"Shaded: the ranges stated in Ref. [25]; dashed: $\Delta S_{\rm tot}=0$ and $\eta_{\rm GSL}=1/2$.")


def fig_scatter(plt, outdir, cfg, S, res):
    fd = res.get("B7_findings")
    an = res.get("B7_scatter_N%s_k%s" % (_f(cfg["anchor"]["N"]), _f(cfg["anchor"]["k_op"])))
    if not fd and not an:
        return None
    fig, ax = _panel_pair(plt)
    if fd:
        k0 = fd["kappa0"]
        w = np.linspace(0.05, 1.6, 200)
        ax[0].semilogy(w, np.exp(-2 * np.pi * w / k0), color=ps.GREY, ls="--", label=r"$e^{-2\pi\omega/\kappa}$")
        P = res.get("B7_scatter_N%s_k%s" % (_f(fd["N"]), _f(fd["k_op"])))
        if P:
            rr = [x for x in P["omega_records"] if x["r_conf"] is not None]
            ax[0].semilogy([x["omega"] for x in rr], [x["r_conf"] for x in rr], color="k", label="solve")
        ax[0].semilogy([x["omega"] for x in fd["rows"]], [x["findings"] for x in fd["rows"]], ls="none", marker="o", mfc="none", color=ps.ORANGE,
                       label="FINDINGS")
        ax[0].set_xlabel(r"$\omega$")
        ax[0].set_ylabel(r"$|S_{Rv\leftarrow Ru}|^2/|S_{Ru\leftarrow Ru}|^2$")
        ax[0].legend(loc="lower left")
    if an:
        k0 = an["summary"]["kappa0"]
        rr = [x for x in an["omega_records"] if x["n_H"] is not None and x["n_H"] > 0]
        if rr:
            w = np.array([x["omega"] for x in rr])
            ax[1].semilogy(w, [x["n_H"] for x in rr], color="k", label=r"$|S_{Lu\leftarrow Rv}|^2$")
            ww = np.linspace(w.min(), w.max(), 200)
            ax[1].semilogy(ww, 1 / np.expm1(2 * np.pi * ww / k0), color=ps.GREY, ls="--", label=r"$(e^{2\pi\omega/\kappa}-1)^{-1}$")
            ax[1].legend(loc="lower left")
        ax[1].set_xlabel(r"$\omega$")
        ax[1].set_ylabel(r"$n_H$")
    return ps.save(plt, fig, outdir, "sweep_fig4_scattering"), (
        "sweep_fig4_scattering", "Figure S4",
        r"Scattering solve of the BdG operator on the soliton-plus-flow background at $\delta_3=0$, $\kappa_g=0.3$ (diagnostic, INV5). (a) Conversion "
        r"ratio on the supersonic side at $(N,k_{\rm op})=(2,3)$ with the FINDINGS values and the Boltzmann factor at the kinematic $\kappa$. "
        r"(b) Upstream positive-norm population from the incoming negative-norm mode at $(1,1.5)$, where $\omega>\mu$ opens the upstream channel, "
        r"with the Planck form at the kinematic $\kappa$.")


def fig_anchor(plt, outdir, cfg, S):
    a = S.get("B10")
    if not a:
        return None
    om = np.asarray(cfg["omega"], dtype=float)
    fig, ax = _panel_pair(plt)
    La = b_census.entanglement_layer(a["kappa_kin"], om)
    op = (S.get("B1") or {}).get("op")
    lo = a["mu"]
    ax[0].axvspan(lo, 1.6, color=ps.BAND, lw=0)
    ax[0].semilogy(om, La["EN"], color="k", label=r"$(N,k_{\rm op})=(1,1.5)$")
    if op:
        Lo = b_census.entanglement_layer(op["kappa_num"], om)
        ax[0].semilogy(om, Lo["EN"], color=ps.ORANGE, ls="--", label=r"$(2,1)$")
    ax[0].set_xlabel(r"$\omega$")
    ax[0].set_ylabel(r"$E_N$ (nats)")
    ax[0].legend(loc="lower left")
    for eta, ls in ((1.0, "-"), (0.9, "--"), (0.5, "-."), (0.1, ":")):
        y = La["EN"] if eta == 1.0 else La["EN_loss_eta%.1f" % eta]
        ax[1].semilogy(om, y, color="k", ls=ls, label=r"$\eta=%g$" % eta)
    ax[1].axvspan(lo, 1.6, color=ps.BAND, lw=0)
    ax[1].set_xlabel(r"$\omega$")
    ax[1].set_ylabel(r"$E_N$ (nats)")
    ax[1].legend(loc="lower left")
    return ps.save(plt, fig, outdir, "sweep_fig5_anchor"), (
        "sweep_fig5_anchor", "Figure S5",
        r"Entanglement layer at the anchor $(N,k_{\rm op},\delta_3,\kappa_g)=(1,1.5,0.05,0.3)$. (a) $E_N(\omega)$ at the kinematic $\kappa$, with the "
        r"operating point $(2,1)$ for comparison; shaded: the Hawking window $[\mu,1.6]$ of the anchor. (b) $E_N$ after symmetric pure loss "
        r"$\eta$, $\nu_-'=1-\eta(1-e^{-2r})$.")


def fig_loss(plt, outdir, cfg, S):
    red = S.get("_B3_red") or {}
    if not red:
        return None
    fig, ax = _panel_pair(plt)
    for d3, mk, c in ((0.02, "o", ps.BLUE), (0.05, "s", ps.RED)):
        pts = sorted([r for r in red.values() if abs(r["d3"] - d3) < 1e-12], key=lambda r: r["N"])
        if not pts:
            continue
        N = np.array([r["N"] for r in pts])
        ax[0].semilogy(N, [max(r["loss_0_40"], 1e-16) for r in pts], marker=mk, mfc="none", color=c, label=r"$\delta_3=%g$" % d3)
        ax[1].plot(N, [r["window_drift"][4] for r in pts], marker=mk, mfc="none", color=c, label=r"$\delta_3=%g$, $[30,40]$" % d3)
        ax[1].plot(N, [r["window_drift"][0] for r in pts], marker=mk, color=c, ls=":", label=r"$\delta_3=%g$, $[0,5]$" % d3)
        NN = np.linspace(N.min(), N.max(), 100)
        ax[1].plot(NN, -d3 * NN ** 2, color=c, ls="--", lw=0.8)
    ax[0].set_xlabel(r"$N$")
    ax[0].set_ylabel(r"$1-Q(40)/Q(0)$")
    ax[0].legend(loc="lower right")
    ax[1].set_xlabel(r"$N$")
    ax[1].set_ylabel(r"$d\tau_c/d\xi$")
    ax[1].legend(loc="lower left")
    return ps.save(plt, fig, outdir, "sweep_fig6_loss_drift"), (
        "sweep_fig6_loss_drift", "Figure S6",
        r"Launched soliton $N\,\mathrm{sech}(N\tau)$ under Eq. (1) (split step, absorbing edges). (a) Core-norm loss over $\xi\in[0,40]$ against $N$. "
        r"(b) Drift velocity of the peak over the first and last windows, with the leading order $-\delta_3N^2$ (dashed).")


def make_all(rundir, cfg, S, res):
    outdir = os.path.join(rundir, "figures")
    os.makedirs(outdir, exist_ok=True)
    plt = ps.setup()
    files, caps = [], []
    for fn in (fig_window_map, fig_raman, fig_ref25, fig_anchor, fig_loss):
        out = fn(plt, outdir, cfg, S)
        if out:
            files += out[0]
            caps.append(out[1])
    out = fig_scatter(plt, outdir, cfg, S, res)
    if out:
        files += out[0]
        caps.append(out[1])
    caps.sort(key=lambda c: c[0])
    if caps:
        files += ps.write_captions(outdir, caps, "SOTHE-P2-SWEEP diagnostic figures (RECORDED, HR-3)")
    return files
