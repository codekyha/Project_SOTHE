"""Numbers quoted in the article "Bogoliubov-de Gennes spectrum, Cherenkov radiation and recoil of solitons with third-order dispersion underlying optical analogue horizons" (H. Oguz, 2026): operating point (N, delta3, kappa_g, k_op) = (1, 0.05, 0.3, 1.5).
Sources: closed forms; SOTHE-P2-SWEEP 1.0.0 job 530429 (run 20260928T131100Z); the job 530047 values (SOTHE-P2 1.0.0), which job
530429 reproduces bitwise; and evaluations of the B1 and B9 code paths of SOTHE-P2-SWEEP (the B1 and B9 tasks of job 530429 are
bitwise identical on the Altay MKL and on OpenBLAS).
Each number is written to article_numbers.json with its source.

Environment variables (or edit the constants below):
  SOTHE_SWEEP_PKG       folder of the SOTHE-P2-SWEEP package (default: ../SOTHE-P2-SWEEP)
  SOTHE_SWEEP_RUN       run folder 20260928T131100Z, unpacked from SOTHE-P2-SWEEP_results_20260928T131100Z.tar.gz
  SOTHE_ARTICLE_NUMBERS output file (default: ./article_numbers.json)
See README.md, or run reproduce.py, which sets all of this up."""
import csv, json, math, os, sys
import numpy as np
HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.environ.get('SOTHE_SWEEP_PKG') or os.path.join(HERE, '..', 'SOTHE-P2-SWEEP'))
from sweep_p2 import b_kin
from sothe_p2 import suite

RUN = os.environ.get('SOTHE_SWEEP_RUN') or sys.exit('set SOTHE_SWEEP_RUN (see README.md)')
S = json.load(open(os.path.join(RUN, 'summary.json')))['results']
T = {}
def put(k, v, src):
    T[k] = {'value': v, 'source': src}

N, d3, kg, kop = 1.0, 0.05, 0.3, 1.5
beta0 = kop**2/2 - d3*kop**3
Sv = beta0/math.sqrt(beta0**2 + kg**2)
vg0 = kop - 3*d3*kop**2
kap = N*Sv
put('beta0_kop', beta0, 'closed form k^2/2 - d3 k^3')
put('S_kop', Sv, 'closed form')
put('vg0_kop', vg0, 'closed form k - 3 d3 k^2')
put('kappa_closed', kap, 'closed form N S(k_op)')
put('T_H_closed', kap/(2*math.pi), 'closed form')
put('HP_closed', 2*math.pi/kap, 'closed form')
put('mu', N*N/2, 'closed form')
cm = abs(vg0)*Sv
put('c_match', cm, 'closed form |v_g0| S'); put('two_c_match', 2*cm, 'closed form')
lam = kap/cm; dt = 40/4095
put('lambda', lam, 'closed form N/|v_g0|'); put('dtau', dt, 'grid')
put('eq7pp_mid', -(7/12)*(lam*dt)**2, 'closed form, s = dtau/2')
put('eq7pp_node', -(1/3)*(lam*dt)**2, 'closed form, s = 0')
put('band_rel', 0.25*(lam*dt)**2, 'closed form'); put('band_abs', 0.25*(lam*dt)**2*kap, 'closed form')
put('cells_per_step', 0.075/(2*N*dt), 'closed form Delta d/(2 N dtau)')
put('tau_h_range', 0.3/(2*N), 'closed form max |tau_h| = 0.3/(2N)')
B10 = S['B10']
put('kappa_num', B10['kappa_kin'], 'job 530429 B10')
put('T_H_num', B10['T_H'], 'job 530429 B10')
put('deficit_obs', B10['deficit_rel'], 'job 530429 B10')
# drive scan
dr = [ [float(x) for x in r] for r in list(csv.reader(open(os.path.join(RUN,'B10_anchor/paper2_data_anchor/kappa_drive_independence.csv'))))[1:] ]
d_ = np.array([r[0] for r in dr]); kd = np.array([r[1] for r in dr])
put('drive_values', d_.tolist(), 'job 530429 B10 CSV'); put('kappa_drive', kd.tolist(), 'job 530429 B10 CSV')
put('drive_spread_abs', float(kd.max()-kd.min()), 'job 530429 B10 CSV')
put('drive_spread_rel', float((kd.max()-kd.min())/kap), 'job 530429 B10 CSV / closed kappa')
pred = np.array([b_kin.eq7pp(N, d3, kg, d, kop)[0] for d in d_])
kpred = kap*(1+pred)
put('eq7pp_drive_max_abs_dev', float(np.max(np.abs(kd-kpred))), 'B10 values against Eq. (7pp)')
sp_obs = kd.max()-kd.min(); sp_pred = kpred.max()-kpred.min()
put('eq7pp_spread_rel_dev', float(abs(sp_obs-sp_pred)/sp_obs), 'B10 values against Eq. (7pp)')
# mirror
mir = max(abs(kd[i]-kd[len(kd)-1-i]) for i in range(len(kd)))
put('mirror_max', float(mir), 'job 530429 B10 CSV (d and 2-d pairs)')
# flow
fl = [ [float(x) for x in r] for r in list(csv.reader(open(os.path.join(RUN,'B10_anchor/paper2_data_anchor/kinematic_kappa_flow.csv'))))[1:] ]
put('flow_kg', [r[0] for r in fl], 'job 530429 B10 CSV'); put('flow_kappa', [r[1] for r in fl], 'job 530429 B10 CSV')
put('flow_monotone', bool(all(fl[i+1][1] < fl[i][1] for i in range(len(fl)-1))), 'job 530429 B10 CSV')
# continuity at k_op = 1.5, kappa_g = 0, N = 1, 2, 3 (B1 code path)
cont = []
for Nc in (1.0, 2.0, 3.0):
    k, aux = b_kin.kappa_extract(Nc, d3, 0.0, 1.0, kop)
    p, l_, dt_, s_ = b_kin.eq7pp(Nc, d3, 0.0, 1.0, kop)
    cont.append({'N': Nc, 'kappa': k, 'extracted_deficit': Nc-k, 'eq7pp_deficit': -p*Nc, 'rel': abs((Nc-k)-(-p*Nc))/(-p*Nc)})
put('continuity', cont, 'B1 code path (kappa_extract) at k_op = 1.5; N = 1, 2 also in job 530429 B1')
# cross-check N = 1, 2 against job 530429 B1 kin_grid.csv
kg_rows = list(csv.DictReader(open(os.path.join(RUN,'B1_kinematics/kin_grid.csv'))))
chk = {}
for r in kg_rows:
    if abs(float(r['k_op'])-1.5)<1e-12 and abs(float(r['d3'])-0.05)<1e-12 and float(r['kg'])==0.0 and float(r['N']) in (1.0, 2.0):
        chk[float(r['N'])] = float(r['kappa_num'])
put('continuity_530429_check', {str(k): (v, [c['kappa'] for c in cont if c['N']==k][0]) for k, v in chk.items()}, 'job 530429 B1 kin_grid.csv vs container')
# parity table at the anchor (B9 code path)
par = b_kin.parity({'odd':[4095,8191,16383,32767,65535],'even':[4096,8192,16384,32768,65536],'cont_N':[1.0,2.0,3.0]}, {'N':N,'d3':d3,'kg':kg,'k_op':kop})
put('parity_table2', [(r['parity'], r['n_tau'], r['deficit_over_lamdt2']) for r in par['table2']], 'B9 code path at the anchor')
# thermality and entanglement
put('thermality_slope', B10['thermality_slope'], 'job 530429 B10'); put('thermality_R2', B10['thermality_R2'], 'job 530429 B10')
put('kappa_fit', B10['kappa_fit'], 'job 530429 B10'); put('flux_residual', B10['flux_residual'], 'job 530429 B10')
def EN(w, k):
    return 2*math.asinh(math.sqrt(1/(math.exp(2*math.pi*w/k)-1)))
for w, key in ((0.05, 'EN_005'), (0.5, 'EN_mu'), (1.6, 'EN_16')):
    put(key, EN(w, kap), 'closed form at closed kappa')
put('nu_max', math.exp(-EN(1.6, kap)), 'closed form')
put('EN_3kappa', EN(3*kap, kap), 'closed form, kappa-independent 2 e^{-3 pi}')
put('IR_asym_005', math.log(2*kap/(math.pi*0.05)), 'closed form ln(2 kappa/pi omega)')
put('IR_dev_005', abs(EN(0.05, kap)-math.log(2*kap/(math.pi*0.05))), 'closed form')
for w, key in ((0.05, 'nbar_peak'), (0.5, 'nbar_mu'), (1.6, 'nbar_band')):
    put(key, (math.exp(EN(w, kap))-1)/2, 'closed form')
def ENeta(w, k, eta):
    r = math.asinh(math.sqrt(1/(math.exp(2*math.pi*w/k)-1)))
    return max(0.0, -math.log(1-eta*(1-math.exp(-2*r))))
put('EN_mu_eta', [ENeta(0.5, kap, e) for e in (1, 0.9, 0.5, 0.1)], 'closed form, loss channel')
en = list(csv.DictReader(open(os.path.join(RUN,'B10_anchor/paper2_data_anchor/partner_log_negativity.csv'))))
put('EN_cross_max', max(abs(float(r['E_N'])-float(r['E_N_cross'])) for r in en), 'job 530429 B10 CSV')
put('band_frac_above_mu', (1.6-0.5)/(1.6-0.05), 'closed form')
put('omega_thermal_n_at_mu', 1/(math.exp(2*math.pi*0.5/kap)-1), 'closed form')
# census
cz = {c['frame']: c for c in B10['census']}
put('census', B10['census'], 'job 530429 B10')
# loss
L = B10['loss']
put('loss_0_40', L['loss_0_40'], 'job 530047 = 530429 B3'); put('loss_HP_after10', L['loss_one_anchor_HP_after_10'], 'job 530429 B10')
put('drift_windows', L['window_drift'], 'job 530047 = 530429 B3')
put('HP_in_40', 40/(2*math.pi/kap), 'closed form')
B3 = S['B3']
put('op_old_loss', {'c': B3['op']['c_over_xi'], 'loss_0_40': B3['op']['loss_0_40'], 'loss_first_hp': B3['op']['loss_first_hp'], 'drift': [B3['op']['window_drift'][0], B3['op']['window_drift'][-1]]}, 'job 530047 = 530429 B3 at (2, 0.05)')
# platform dictionary at N = 1
T0 = 50e-15; b2 = 15e-27; g = 0.1e-3
P0 = N**2*b2/(g*T0**2)
put('P0_W', P0, 'closed form N^2 |beta2|/(gamma T0^2)'); put('E_pulse_J', 2*P0*T0, 'closed form 2 P0 T0')
put('f_kop_THz', kop/(2*math.pi*T0)/1e12, 'closed form'); put('L_D_m', T0**2/b2, 'closed form'); put('beta3_ps3_km', 6*15*0.05*d3, 'closed form 6 |beta2| T0 delta3')
json.dump(T, open(os.environ.get('SOTHE_ARTICLE_NUMBERS') or os.path.join(os.getcwd(), 'article_numbers.json'), 'w'), indent=1, default=float)
for k, v in T.items():
    val = v['value']
    s = json.dumps(val, default=float)
    print('%-26s %s' % (k, s[:150]))
