"""Figures 1, 2, 4, 5, 6 and S1 of the article "Bogoliubov-de Gennes spectrum, Cherenkov radiation and recoil of solitons with third-order dispersion underlying optical analogue horizons" (H. Oguz, 2026), from the archived runs of SOTHE-P2 1.0.0 and SOTHE-P2-SWEEP 1.0.0.
The file names and the figure numbers in this script follow the record text: fig2_bdg_krein is Fig. 1, fig3_radiative_loss Fig. 2,
fig1_kinematics Fig. 4, fig4_channels Fig. 5 and fig5_entanglement Fig. 6 of the article; Fig. 3 is drawn by SOTHE-P2-RAD 1.0.0
(rad.py figures).

Operating point (N, delta3, kappa_g, k_op) = (1, 0.05, 0.3, 1.5). Style: sothe_p2.pubstyle (IOP widths 8.5 / 15 cm, Computer
Modern, no titles, no annotation).
Data: SOTHE-P2-SWEEP 1.0.0, UHeM Altay job 530429 (run 20260928T131100Z): the B10 anchor CSV files, the B3 series and points, the
B4 census, the B5 spectrum at delta3 = 0.05 and the B7 scattering calculation; closed forms (surface gravity, E_N, the boost edge).
Fig. 2(a): SOTHE-P2 1.0.0, job 530047 (run 20260927T175224Z), stage S5 figure data.

Environment variables (or edit the constants below):
  SOTHE_SWEEP_PKG    folder of the SOTHE-P2-SWEEP package (default: ../SOTHE-P2-SWEEP)
  SOTHE_SWEEP_RUN    run folder 20260928T131100Z, unpacked from SOTHE-P2-SWEEP_results_20260928T131100Z.tar.gz
  SOTHE_P2_FIGDATA   S5_figures/figure_data.npz, unpacked from SOTHE-P2_results_20260927T175224Z.tar.gz
  SOTHE_ARTICLE_OUT  output folder (default: ./article_output)
See README.md, or run reproduce.py, which sets all of this up."""
import csv, json, math, os, shutil, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get('SOTHE_SWEEP_PKG') or os.path.join(HERE, '..', 'SOTHE-P2-SWEEP'))
from sothe_p2 import pubstyle as S
from sweep_p2 import b_kin

RUN = os.environ.get('SOTHE_SWEEP_RUN') or sys.exit('set SOTHE_SWEEP_RUN (see README.md)')
OUT = os.environ.get('SOTHE_ARTICLE_OUT') or os.path.join(os.getcwd(), 'article_output')
os.makedirs(OUT, exist_ok=True)
SUM = json.load(open(os.path.join(RUN, 'summary.json')))['results']
N, D3, KG, KOP = 1.0, 0.05, 0.3, 1.5
beta0 = KOP**2/2 - D3*KOP**3
KAP = N*beta0/math.sqrt(beta0**2 + KG**2)
plt = S.setup()
FD = {}

def rcsv(*p):
    return list(csv.DictReader(open(os.path.join(*p))))

# ------------------------------------------------------------------ Fig. 1: kinematic layer (a) drive scan, (b) flow
dr = rcsv(RUN, 'B10_anchor/paper2_data_anchor/kappa_drive_independence.csv')
d = np.array([float(r['drive']) for r in dr]); kd = np.array([float(r['kappa_kin']) for r in dr])
dd = np.linspace(0.7, 1.3, 6001)
saw = np.array([b_kin.eq7pp(N, D3, KG, x, KOP)[0] for x in dd])
lam, dt = b_kin.eq7pp(N, D3, KG, 1.0, KOP)[1:3]
lo, hi = -(7/12)*(lam*dt)**2, -(1/3)*(lam*dt)**2
fl = rcsv(RUN, 'B10_anchor/paper2_data_anchor/kinematic_kappa_flow.csv')
g = np.array([float(r['kappa_g']) for r in fl]); kf = np.array([float(r['kappa_kin']) for r in fl])
gg = np.linspace(0, 1, 401)
kcl = N*beta0/np.sqrt(beta0**2 + gg**2)
f = S.figure(plt, S.W2, 5.6)
a, b = f.subplots(1, 2)
xs = np.array([b_kin.eq7pp(N, D3, KG, x, KOP)[3] for x in d]) / dt
xx = np.linspace(0, 1, 401)
a.plot(xx, -1e5 * (1/3 + xx*(1 - xx)) * (lam*dt)**2, '-', color='k', lw=0.8, zorder=2)
a.plot(xs, 1e5*(kd/KAP - 1), 's', mfc=S.BLUE, mec=S.BLUE, ls='none', clip_on=False, zorder=3)
a.set_xlim(0, 1)
a.set_ylim(-4.4, -2.2)
a.set_xlabel(r'$s/\Delta\tau$')
a.set_ylabel(r'$10^{5}\,(\kappa_{\rm num}/\kappa-1)$')
S.label(a, '(a)'); S.frame(a)
b.plot(gg, kcl, '-', color='k', lw=0.8, zorder=2)
b.plot(g, kf, 'o', mfc=S.ORANGE, mec=S.ORANGE, ls='none', clip_on=False, zorder=3)
b.plot([0, 1], [N, N], 'k:', lw=0.8)
b.plot([0], [N], 's', ms=5.5, mec=S.RED, mfc=S.RED, ls='none', clip_on=False, zorder=4)
b.set_xlim(-0.03, 1.03)
yl = (0.65, 1.02)
b.set_ylim(*yl)
b.set_xlabel(r'$\kappa_g$'); b.set_ylabel(r'$\kappa$', rotation=0, labelpad=8)
S.frame(b, mirror=False); b.tick_params(which='both', top=True)
b2 = b.twinx(); b2.set_ylim(yl[0]/(2*math.pi), yl[1]/(2*math.pi)); b2.set_ylabel(r'$T_H=\kappa/2\pi$ (local)')
b2.tick_params(which='both', direction='in', right=True, labelsize=S.FS_SMALL); S.nice_ticks(b2)
for sp in b2.spines.values(): sp.set_linewidth(S.LW_AX)
S.label(b, '(b)')
S.save(plt, f, OUT, 'fig1_kinematics')
FD.update(fig1_drive=d, fig1_kappa=kd, fig1_kappa_closed=KAP, fig1_band=[lo, hi], fig1_kg=g, fig1_kappa_flow=kf)

# ------------------------------------------------------------------ Fig. 2: BdG. Panel (a): job 530047 (delta3 = 0 about psi0,
# S5 figure data of job 530047); panel (b): job 530429 block B5 at the operating value delta3 = 0.05 about phi
# (full spectrum with Krein signs, zero sector, psi0 quartet). Plotting only, with the marker convention of SOTHE-P2 1.0.0.
from sothe_p2 import figures as F2
Z = np.load(os.environ.get('SOTHE_P2_FIGDATA') or sys.exit('set SOTHE_P2_FIGDATA (see README.md)'), allow_pickle=False)
B5 = json.load(open(os.path.join(RUN, 'tasks', 'B5', 'B5_bdg_N1_d0.05.json')))['result']
wf = np.array(B5['spectrum']['re']) + 1j * np.array(B5['spectrum']['im'])
kf = np.array(B5['spectrum']['krein'], dtype=float)
zs = np.array([z['re'] + 1j * z['im'] for z in B5['summary']['zero_sector']])
q0 = np.array([z['re'] + 1j * z['im'] for z in B5['summary']['psi0_quartet']])
f2, info2 = F2.plot_fig3_bdg_krein(plt, Z['fig2_w_psi0_d3_0'], Z['fig2_krein_psi0_d3_0'], float(Z['fig2_mu']), wf, kf, 0.5, q0, zs)
ai = [a for a in f2.axes[1].child_axes][0]
ai.set_xlim(-0.2, 0.2); ai.set_ylim(-0.2, 0.2); ai.set_xticks([-0.1, 0.0, 0.1]); ai.set_yticks([-0.1, 0.0, 0.1])
S.save(plt, f2, OUT, 'fig2_bdg_krein')
pos, neg, neu = F2.krein_classes(wf, kf)
ann = (np.abs(wf.real) > 0.05) & (np.abs(wf.real) <= 1.5)
negE = ann & ((pos & (wf.real < 0)) | (neg & (wf.real > 0)))
print('fig2 (b): annulus', int(ann.sum()), 'pos', int((ann & pos).sum()), 'neg', int((ann & neg).sum()),
      'negative-energy', int(negE.sum()), np.round(wf[negE].real, 4).tolist(), 'maxIm annulus', float(np.max(np.abs(wf[ann].imag))), info2)
FD.update(fig2b_w=wf, fig2b_krein=kf, fig2b_zero_sector=zs, fig2b_psi0_quartet=q0)

# ------------------------------------------------------------------ Fig. 3: radiative loss
def series(n, d3):
    r = rcsv(RUN, 'B3_loss/series', 'loss_N%g_d%g.csv' % (n, d3))
    return {k: np.array([float(x[k]) for x in r]) for k in ('xi', 'Q', 'tau_c')}
Y2, Y1 = series(2, 0.05), series(1, 0.05)
pts = rcsv(RUN, 'B3_loss', 'loss_points.csv')
EQ10 = {(1.0, 0.02), (1.0, 0.03), (1.0, 0.04), (1.0, 0.05), (1.0, 0.06), (1.0, 0.08), (1.0, 0.1), (2.0, 0.025), (2.0, 0.05)}
P = [(float(p['N']), float(p['d3']), float(p['loss_0_40']), float(p['c_over_xi'])) for p in pts if (float(p['N']), round(float(p['d3']), 4)) in EQ10]
assert len(P) == 9
c2 = [p[3] for p in P if p[0] == 2.0 and abs(p[1]-0.05) < 1e-12][0]
f = S.figure(plt, S.W1, 13.0)
a, b, c = f.subplots(3, 1)
from matplotlib.ticker import FixedFormatter, FixedLocator, NullFormatter, NullLocator
for Y, ls, col, z in ((Y2, '-', 'k', 3), (Y1, '-.', S.ORANGE, 4)):
    s = Y['xi'] >= 0.5
    a.loglog(Y['xi'][s], Y['Q'][s]/Y['Q'][0], ls, color=col, zorder=z)
i10 = int(np.argmin(np.abs(Y2['xi'] - 10)))
xx = np.linspace(10, 40, 50)
a.loglog(xx, (Y2['Q'][i10]/Y2['Q'][0])*(xx/10)**(-c2), '--', color=S.RED, zorder=5)
a.set_xlim(0.5, 40)
ylo = math.floor(float(np.min(Y2['Q'][Y2['xi'] >= 0.5]/Y2['Q'][0]))*200)/200
a.set_ylim(ylo, 1.002)
yt = np.arange(math.ceil(ylo*100), 101)/100
xt = [0.5, 1, 2, 5, 10, 20, 40]
a.yaxis.set_major_locator(FixedLocator(yt)); a.yaxis.set_major_formatter(FixedFormatter(['%.2f' % v for v in yt]))
a.yaxis.set_minor_locator(NullLocator()); a.yaxis.set_minor_formatter(NullFormatter())
a.xaxis.set_major_locator(FixedLocator(xt)); a.xaxis.set_major_formatter(FixedFormatter(['%g' % v for v in xt]))
a.xaxis.set_minor_formatter(NullFormatter())
a.set_xlabel(r'$\xi$'); a.set_ylabel(r'$Q(\xi)/Q(0)$'); S.label(a, '(a)'); S.frame(a)
b.plot(Y2['xi'], Y2['tau_c'], 'k-', zorder=3)
b.plot(Y1['xi'], Y1['tau_c'], '-.', color=S.ORANGE, zorder=4)
b.set_xlim(0, 40); b.set_xlabel(r'$\xi$'); b.set_ylabel(r'$\tau_c$', rotation=0, labelpad=8); S.label(b, '(b)'); S.frame(b)
n1 = [p for p in P if p[0] == 1.0]; n2 = [p for p in P if p[0] == 2.0]
h2, = c.semilogy([p[0]*p[1] for p in n2], [p[2] for p in n2], 's', ms=3.5, mec='k', mfc='k', ls='none', clip_on=False, zorder=3)
h1, = c.semilogy([p[0]*p[1] for p in n1], [p[2] for p in n1], 'o', ms=5.5, mec='k', mfc='none', mew=0.7, ls='none', clip_on=False, zorder=4)
op = [p for p in n1 if abs(p[1]-0.05) < 1e-12][0]
h3, = c.semilogy([op[0]*op[1]], [op[2]], 'o', ms=10.0, mec=S.ORANGE, mfc=S.ORANGE, mew=0.7, ls='none', clip_on=False, zorder=2)  # larger disk behind the N=1 circle and the N=2 square
c.set_xlim(0, 0.11); c.set_xlabel(r'$N\delta_3$'); c.set_ylabel(r'$1-Q(40)/Q(0)$'); S.label(c, '(c)'); S.frame(c)
c.legend([h1, h2, h3], [r'$N=1$', r'$N=2$', r'operating point'], loc='lower right')
S.save(plt, f, OUT, 'fig3_radiative_loss')
FD.update(fig3_xi=Y2['xi'], fig3_Q_N2=Y2['Q'], fig3_tau_c_N2=Y2['tau_c'], fig3_Q_N1=Y1['Q'], fig3_tau_c_N1=Y1['tau_c'], fig3_c=c2,
          fig3_points=np.array([[p[0], p[1], p[2]] for p in P]))

# ------------------------------------------------------------------ Fig. 4: channels over (N, omega)
WT = {(r['N'], r['k_op']): r for r in SUM['B4']['window_table']}
Ns = np.array(sorted({r['N'] for r in SUM['B4']['window_table']}))
e1 = np.array([WT[(n, 1.0)]['late_edge'] for n in Ns]); e15 = np.array([WT[(n, 1.5)]['late_edge'] for n in Ns])
l15 = np.array([WT[(n, 1.5)]['lab_edge'] for n in Ns])
nn = np.linspace(0.5, 2.1, 400)
f = S.figure(plt, S.W1, 8.6)
a = f.add_subplot(1, 1, 1)
# window at k_op = 1.5: mu < omega < min(edge, 1.6), piecewise linear edge between grid points
edge_i = np.interp(nn, Ns, e15)
top = np.minimum(edge_i, 1.6); bot = np.maximum(nn**2/2, 0.05)
m = (nn >= Ns[0]) & (nn <= Ns[-1]) & (top > bot)
a.fill_between(nn[m], bot[m], top[m], color=S.FILL, lw=0, zorder=0)
a.plot(nn, nn**2/2, 'k-', lw=1.0, zorder=2)
a.plot([0.5, 2.1], [1.6, 1.6], 'k:', lw=0.8, zorder=2)
a.plot([0.5, 2.1], [0.05, 0.05], 'k:', lw=0.8, zorder=2)
h15, = a.plot(Ns, e15, 'o-', color=S.ORANGE, mec=S.ORANGE, mfc=S.ORANGE, lw=0.9, zorder=4)
hl, = a.plot(Ns, l15, 'o:', color=S.ORANGE, mec=S.ORANGE, mfc='none', lw=0.8, zorder=4)
h1, = a.plot(Ns, e1, 's--', color=S.BLUE, mec=S.BLUE, mfc=S.BLUE, lw=0.9, zorder=4)
# Galilean boost at k_op = 1.5 with the recorded late drifts (xi in [30, 40]) of block B3: closed-form edge (S2) with
# V = 2 c_match - V_d and the long-wavelength threshold min_k [k^2/2 + mu - d3 k^3 - V_d k] (plotting of closed forms)
from scipy.optimize import minimize_scalar
LP = {(float(r_['N']), round(float(r_['d3']), 4)): r_ for r_ in rcsv(RUN, 'B3_loss', 'loss_points.csv')}
b15 = 1.5 ** 2 / 2 - D3 * 1.5 ** 3; S15 = b15 / math.sqrt(b15 ** 2 + KG ** 2); c15 = abs(1.5 - 3 * D3 * 1.5 ** 2) * S15
def edge_boost(V, mu):
    ks = (math.sqrt(1 + 12 * D3 * V) - 1) / (6 * D3)
    return V * ks - 0.5 * ks ** 2 - mu - D3 * ks ** 3
Vd_late = np.array([float(LP[(float(n), 0.05)]['window_drift_5']) for n in Ns])
eb15 = np.array([edge_boost(2 * c15 - v, n * n / 2) for n, v in zip(Ns, Vd_late)])
thb = np.array([minimize_scalar(lambda k, n=n, v=v: 0.5 * k * k + n * n / 2 - D3 * k ** 3 - v * k, bounds=(-1 / (3 * D3), 1 / (3 * D3)),
                                method='bounded', options={'xatol': 1e-12}).fun for n, v in zip(Ns, Vd_late)])
hb, = a.plot(Ns, eb15, '^-.', color=S.RED, mec=S.RED, mfc=S.RED, ms=4.0, lw=0.9, zorder=4)
ht, = a.plot(Ns, thb, '+--', color='k', ms=5.0, lw=0.8, zorder=3)
a.set_xlim(0.5, 2.1); a.set_ylim(0, 2.8)
a.set_xlabel(r'$N$'); a.set_ylabel(r'$\omega$', rotation=0, labelpad=8)
f.legend([h15, hb, hl, h1, ht], [r'$k_{\rm op}=1.5$, drift-adjusted', r'$k_{\rm op}=1.5$, boost', r'$k_{\rm op}=1.5$, laboratory',
                                r'$k_{\rm op}=1$, drift-adjusted', r'threshold, boost'], loc='outside lower center', ncol=2,
         handlelength=2.4, columnspacing=1.0, labelspacing=0.3)
S.frame(a)
S.save(plt, f, OUT, 'fig4_channels')
print('fig4 boost edges', np.round(eb15, 3).tolist(), 'thresholds', np.round(thb, 4).tolist())
FD.update(fig4_N=Ns, fig4_edge_k1=e1, fig4_edge_k15=e15, fig4_lab_edge_k15=l15, fig4_boost_edge_k15=eb15, fig4_boost_threshold=thb,
          fig4_late_drift=Vd_late)

# ------------------------------------------------------------------ Fig. 5: entanglement at the operating point
w = np.linspace(0.05, 1.6, 60)
b2_ = 1/(np.exp(2*np.pi*w/KAP) - 1)
r = np.arcsinh(np.sqrt(b2_))
EN = 2*r
EN5 = np.maximum(0, -np.log(1 - 0.5*(1 - np.exp(-2*r))))
f = S.figure(plt, S.W1, 5.6)
a = f.add_subplot(1, 1, 1)
mu = N*N/2
wf = np.linspace(mu, 1.6, 200)
ENf = 2*np.arcsinh(np.sqrt(1/(np.exp(2*np.pi*wf/KAP) - 1)))
a.fill(np.concatenate([wf, wf[::-1]]), np.concatenate([ENf, np.zeros(wf.size)]), color=S.FILL, ec='none', zorder=0)
h1, = a.plot(w, EN, '-', color=S.ORANGE, lw=1.4, zorder=3)
h2, = a.plot(w, EN5, '--', color='k', lw=0.9, zorder=3)
a.set_xlim(0, 1.65); a.set_ylim(0, 1.08*float(np.max(EN)))
a.set_xlabel(r'$\omega$'); a.set_ylabel(r'$E_N$ (nats)')
a.legend([h1, h2], [r'$\eta=1$', r'$\eta=0.5$'], loc='upper right')
S.frame(a)
S.save(plt, f, OUT, 'fig5_entanglement')
FD.update(fig5_omega=w, fig5_EN=EN, fig5_EN_eta05=EN5, fig5_kappa=KAP)
np.savez_compressed(os.path.join(OUT, 'figure_data.npz'), **{k: np.asarray(v) for k, v in FD.items()})
print('kappa', KAP, 'band', lo, hi, 'c2', c2, 'op loss', op)
print(sorted(os.listdir(OUT)))

# ------------------------------------------------------------------ Fig. S1 (Supplement): scattering solve (job 530429 B7, plotting only)
SC = rcsv(RUN, 'B7_scatter', 'scatter_omega.csv')
def pt(n, k):
    return [r for r in SC if abs(float(r['N']) - n) < 1e-9 and abs(float(r['k_op']) - k) < 1e-9]
A = pt(2.0, 3.0); B = [r for r in pt(1.0, 1.5) if r['two_sided'] in ('1', 'True', 'true')]
wa = np.array([float(r['omega']) for r in A]); ra = np.array([float(r['r_conf']) for r in A]); ba = np.array([float(r['boltzmann']) for r in A])
wb = np.array([float(r['omega']) for r in B]); nb = np.array([float(r['n_H']) for r in B]); pb = np.array([float(r['planck']) for r in B])
f = S.figure(plt, S.W2, 5.6)
a, b = f.subplots(1, 2)
a.semilogy(wa, ra, 'k-', zorder=3)
a.semilogy(wa, ba, '--', color=S.GREY, zorder=2)
a.set_xlim(0, 1.65); a.set_ylim(1e-3, 1.5)
a.set_xlabel(r'$\omega$'); a.set_ylabel(r'$|S_{Rv\leftarrow Ru}|^2/|S_{Ru\leftarrow Ru}|^2$')
S.label(a, '(a)'); S.frame(a)
b.semilogy(wb, nb, 'k-', zorder=3)
b.semilogy(wb, pb, '--', color=S.GREY, zorder=2)
b.set_xlim(0.45, 1.65); b.set_ylim(1e-5, 1e-1)
b.set_xlabel(r'$\omega$'); b.set_ylabel(r'$n_H$', rotation=0, labelpad=10)
S.label(b, '(b)'); S.frame(b)
S.save(plt, f, OUT, 'figS1_scattering')
print('figS1: (2,3)', len(wa), 'points, r_conf', ra.min(), ra.max(), '; (1,1.5) two-sided', len(wb), 'points, nH/planck', (nb/pb).min(), (nb/pb).max())
