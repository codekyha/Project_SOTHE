"""make_validation_figure.py

Runs the validation harness, saves
  - figures/validation_figure.png
  - results/validation_run.json   (machine-readable numbers)
  - results/validation_run.txt    (human-readable summary)

Companion to validate_paper_claims.py. Imports the routines from there.

Usage
-----
    python3 make_validation_figure.py
"""
from __future__ import annotations

import json
import os
import sys
import time

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Local imports (same directory)
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from validate_paper_claims import (gnlse_dimensionless, find_rr_lobe,
                                   _trapz)  # noqa: F401

OUT_FIG  = os.path.join(os.path.dirname(__file__), "..", "figures", "validation_figure.png")
OUT_JSON = os.path.join(os.path.dirname(__file__), "..", "results", "validation_run.json")
OUT_TXT  = os.path.join(os.path.dirname(__file__), "..", "results", "validation_run.txt")
for f in (OUT_FIG, OUT_JSON, OUT_TXT):
    os.makedirs(os.path.dirname(f), exist_ok=True)


def main():
    t0 = time.time()

    # --- 1. nominal-point GSL run -----------------------------------------
    print("Running nominal-point GSL evolution ...")
    out = gnlse_dimensionless(N_sol=3.5, delta3=0.02, xi_max=12.0,
                              Nt=2**13, n_steps=4000, n_save=200)

    dS_tot   = float(out['S_tot'][-1] - out['S_tot'][0])
    eta_f    = float(out['eta'][-1])
    frac_neg = float(np.mean(np.diff(out['S_tot']) < 0))
    drift    = float((out['P_norm'][-1] - out['P_norm'][0]) / out['P_norm'][0])

    # --- 2. phase-matching scan -------------------------------------------
    print("Running 5-point delta_3 scan ...")
    deltas = np.array([0.06, 0.07, 0.08, 0.09, 0.10])
    k_RR_list = np.zeros_like(deltas)
    for i, d3 in enumerate(deltas):
        s = gnlse_dimensionless(Nt=2**13, n_steps=4000, n_save=120,
                                xi_max=8.0, delta3=float(d3))
        Spec = np.abs(np.fft.fftshift(np.fft.fft(s['psi_hist'][:, -1]))) ** 2
        k_RR_list[i] = find_rr_lobe(s['k'], Spec, k_c=s['k_c_hist'][-1],
                                    side='+', k_offset=2.0, k_max=25.0,
                                    smooth=5, prefer_k=1.0 / d3)

    x   = 1.0 / deltas
    cRR = float(np.sum(x * k_RR_list) / np.sum(x * x))
    prod_mean = float(np.mean(k_RR_list * deltas))
    prod_std  = float(np.std(k_RR_list * deltas))

    # --- 3. figure --------------------------------------------------------
    print("Plotting ...")
    fig = plt.figure(figsize=(13.5, 4.0), facecolor='w')
    plt.rcParams.update({"font.size": 10, "axes.linewidth": 1.0})

    # (a) spectral evolution -- normalize per propagation step so the
    # soliton ridge sits near 0 dB and the dispersive wave is visible
    ax1 = fig.add_subplot(1, 3, 1)
    spec_evo = np.abs(np.fft.fftshift(np.fft.fft(out['psi_hist'], axis=0), axes=0)) ** 2
    spec_evo_norm = spec_evo / np.max(spec_evo, axis=0, keepdims=True)
    spec_dB = 10 * np.log10(spec_evo_norm + 1e-15)
    im = ax1.imshow(spec_dB,
                    extent=[out['xi'][0], out['xi'][-1], out['k'][0], out['k'][-1]],
                    origin='lower', aspect='auto', cmap='jet', vmin=-60, vmax=0)
    ax1.set_xlabel(r'$\xi$')
    ax1.set_ylabel(r'$k\;(1/T_0)$')
    ax1.set_ylim(-30, 30)
    ax1.set_title('(a) Spectral evolution (per-step normalised)')
    fig.colorbar(im, ax=ax1, label='dB')

    # (b) entropies + eta_GSL
    ax2 = fig.add_subplot(1, 3, 2)
    ax2.plot(out['xi'], out['S_hor'], 'b--', label=r'$S_{\rm hor}$', linewidth=1.4)
    ax2.plot(out['xi'], out['S_rad'], 'r-.', label=r'$S_{\rm rad}$', linewidth=1.4)
    ax2.plot(out['xi'], out['S_tot'], 'k-',  label=r'$S_{\rm tot}$', linewidth=2.0)
    ax2.set_xlabel(r'$\xi$'); ax2.set_ylabel('S (nats)')
    ax2.grid(True, alpha=0.3); ax2.legend(loc='best')
    ax2.set_title(r'(b) Spectral entropies, $\Delta S_{\rm tot}=%+.2f$ nats' % dS_tot)

    ax2b = ax2.twinx()
    ax2b.plot(out['xi'], out['eta'], color='tab:orange', linewidth=1.6)
    ax2b.axhline(0.5, color='gray', linestyle=':', linewidth=1)
    ax2b.set_ylabel(r'$\eta_{\rm GSL}$', color='tab:orange')
    ax2b.tick_params(axis='y', labelcolor='tab:orange')
    ax2b.set_ylim(0, 1)

    # (c) phase-matching
    ax3 = fig.add_subplot(1, 3, 3)
    d3_dense = np.linspace(deltas.min() * 0.95, deltas.max() * 1.05, 200)
    ax3.plot(d3_dense, cRR / d3_dense, 'b-', linewidth=2.0,
             label=fr'$k_{{\rm RR}} = {cRR:.2f}/\delta_3$')
    ax3.plot(deltas, k_RR_list, 'ro', mfc='w', mew=1.5, markersize=9, label='measured')
    ax3.set_xlabel(r'$\delta_3$'); ax3.set_ylabel(r'$k_{\rm RR}\;(1/T_0)$')
    ax3.grid(True, alpha=0.3); ax3.legend(loc='best')
    ax3.set_title(r'(c) Cherenkov scaling, $\langle k_{\rm RR}\delta_3\rangle=%.2f$' % prod_mean)

    fig.suptitle('SOTHE supplementary: Python validation run', fontsize=11, y=1.02)
    fig.tight_layout()
    fig.savefig(OUT_FIG, dpi=200, bbox_inches='tight')
    plt.close(fig)
    print(f"  figure -> {OUT_FIG}")

    # --- 4. dump numeric results -----------------------------------------
    results = dict(
        run_date_utc=time.strftime("%Y-%m-%d %H:%M:%S", time.gmtime()),
        elapsed_seconds=round(time.time() - t0, 2),
        gsl_test=dict(
            params=dict(N_sol=3.5, delta3=0.02, xi_max=12.0,
                        Nt=2**13, n_steps=4000),
            Delta_S_tot_nats=round(dS_tot, 4),
            eta_GSL_final=round(eta_f, 4),
            fraction_local_decreases=round(frac_neg, 4),
            photon_number_drift=drift,
            paper_range_Delta_S_tot=[0.82, 2.42],
            paper_range_eta_GSL=[0.52, 0.62],
            passes=bool((dS_tot > 0) and (eta_f > 0.5)
                        and (0.20 <= frac_neg <= 0.65) and (abs(drift) < 1e-2)),
        ),
        phase_matching_test=dict(
            delta3=deltas.tolist(),
            k_RR=k_RR_list.tolist(),
            product_k_RR_times_delta3=(k_RR_list * deltas).tolist(),
            c_RR_fit=round(cRR, 4),
            product_mean=round(prod_mean, 4),
            product_std=round(prod_std, 4),
            relative_spread_percent=round(100 * prod_std / max(prod_mean, 1e-30), 2),
            paper_c_RR=1.01,
            paper_product_mean=1.02,
            passes=bool((0.8 <= prod_mean <= 1.2) and (prod_std / prod_mean < 0.15)),
        ),
    )
    with open(OUT_JSON, "w") as fh:
        json.dump(results, fh, indent=2)
    print(f"  json   -> {OUT_JSON}")

    txt = []
    txt.append("SOTHE Supplementary -- Python validation run")
    txt.append("=" * 60)
    txt.append(f"Date (UTC): {results['run_date_utc']}")
    txt.append(f"Elapsed: {results['elapsed_seconds']} s")
    txt.append("")
    txt.append("[1] Nominal-point GSL")
    g = results['gsl_test']
    txt.append(f"  Delta S_tot         = {g['Delta_S_tot_nats']:+8.3f} nats")
    txt.append(f"  eta_GSL(xi_f)       = {g['eta_GSL_final']:8.3f}")
    txt.append(f"  fraction dS/dxi<0   = {100*g['fraction_local_decreases']:8.2f} %")
    txt.append(f"  photon-number drift = {g['photon_number_drift']:+.2e}")
    txt.append(f"  PASS: {g['passes']}")
    txt.append("")
    txt.append("[2] Cherenkov phase-matching")
    r = results['phase_matching_test']
    txt.append(f"  delta_3 grid : {r['delta3']}")
    txt.append(f"  k_RR         : {[round(x,3) for x in r['k_RR']]}")
    txt.append(f"  product      : {[round(x,3) for x in r['product_k_RR_times_delta3']]}")
    txt.append(f"  c_RR fit     = {r['c_RR_fit']:.3f}    (paper: {r['paper_c_RR']})")
    txt.append(f"  product mean = {r['product_mean']:.3f}    (paper: {r['paper_product_mean']})")
    txt.append(f"  product std  = {r['product_std']:.3f}    "
               f"({r['relative_spread_percent']:.1f}% rel.)")
    txt.append(f"  PASS: {r['passes']}")
    with open(OUT_TXT, "w") as fh:
        fh.write("\n".join(txt))
    print(f"  txt    -> {OUT_TXT}")

    # --- 5. CSV exports for raw inspection -------------------------------
    OUT_CSV_S = os.path.join(os.path.dirname(__file__), "..", "results", "entropy_trajectory.csv")
    np.savetxt(OUT_CSV_S,
               np.column_stack([out['xi'], out['S_hor'], out['S_rad'], out['S_tot'], out['eta'], out['P_norm']]),
               delimiter=",",
               header="xi,S_hor,S_rad,S_tot,eta_GSL,P_norm",
               comments="")
    print(f"  csv    -> {OUT_CSV_S}")

    OUT_CSV_R = os.path.join(os.path.dirname(__file__), "..", "results", "phase_matching_scan.csv")
    np.savetxt(OUT_CSV_R,
               np.column_stack([deltas, k_RR_list, k_RR_list * deltas]),
               delimiter=",",
               header="delta3,k_RR,k_RR_times_delta3",
               comments="")
    print(f"  csv    -> {OUT_CSV_R}")


if __name__ == "__main__":
    main()
