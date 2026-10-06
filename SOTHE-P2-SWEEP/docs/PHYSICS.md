# SOTHE-P2-SWEEP: equations, protocols and checks of each block

Every number of the sweep is RECORDED Python output (HR-3: only MATLAB R2025b output is canonical). Block B7 is diagnostic: scattering-derived $\alpha$ and $\beta$ are not Paper-2 claims (INV5).

## Conventions

The manuscript's model, Eq. (1), in the retarded frame with $\tau=-T/T_0$ (so that $\beta_3>0$ puts the resonant radiation at $k>0$):

$$\partial_\xi\psi=i\left(\tfrac12\psi_{\tau\tau}+|\psi|^2\psi+i\delta_3\psi_{\tau\tau\tau}\right),\qquad \delta_3=\frac{\beta_3}{6|\beta_2|T_0},$$

with $\partial_\tau\leftrightarrow ik$ under the FFT. Launched soliton $\psi_0=N\,\mathrm{sech}(N\tau)$, $\mu=N^2/2$.

Kinematic layer (suite v0.1.0): $\beta_0(k)=\tfrac12k^2-\delta_3k^3$, $v_{g0}=k-3\delta_3k^2$, $S=|\beta_0|/\sqrt{\beta_0^2+\kappa_g^2}$, $c=\max(|v_{g0}(k_{\rm op})|S,10^{-3})$, $\kappa_0=NS(k_{\rm op})$, flow $v(\tau)=c\,[1+\tanh(\kappa_0(\tau-\tau_h)/c)]$ with $\tau_h=-(d-1)/2N$, extracted on 4096 nodes of $[-20,20]$.

## Grids (full sweep)

| Axis | Values |
|---|---|
| $N$ | 0.6, 0.8, …, 2.0 (8) |
| $k_{\rm op}$ | 1, 1.25, 1.5, 1.75, 2, 2.5, 3 (7) |
| $\delta_3$ (kinematics, census) | 0, 0.02, 0.05 |
| $\kappa_g$ | 0, 0.3, 0.6 |
| band $\omega$ | 60 points on $[0.05,1.6]$ (census: 311 points) |
| drive | 9 points on $[0.7,1.3]$; $\kappa_g$ flow: 21 points on $[0,1]$ |

`python sweep.py submit --fast` uses a reduced grid (smoke test, not protocol values).

## B1. Kinematic layer

For every $(N,k_{\rm op},\delta_3,\kappa_g)$: $\kappa_0$, the extracted $\kappa$ at drive 1, the drive scan and its mirror $d\to2-d$, the $\kappa_g$ flow, and Eq. (7″):

$$\frac{\kappa_{\rm num}}{\kappa_0}-1=-\lambda^2\left[\frac{\Delta\tau^2}{3}+s(\Delta\tau-s)\right],\qquad \lambda=\frac{\kappa_0}{c},$$

$s$ being the offset of the crossing above the node below. The drive spread is bounded by the band $\lambda^2\Delta\tau^2/4$.

## B9. Parity tables

Supplement Table 2: $(\kappa_{\rm num}/\kappa_0-1)/(\lambda\Delta\tau)^2$ at drive 1 for odd and even $n_\tau$ (4095 … 65535, 4096 … 65536). Table 3: $N-\kappa_{\rm num}$ at $\kappa_g=0$ for $N=1,2,3$ against Eq. (7″).

## B3. Radiative loss and drift

`kit.splitstep_loss` of SOTHE-P2 1.0.0 (symmetric split step; absorber $e^{-\sigma\,d\xi}$ beyond $|\tau|>0.35L$; $d\xi=0.002$, $n=8192$, $L=400$, $\xi\le40$; records every 0.5; sub-grid peak). Points: $N\times\{0.02,0.05\}$ and the nine points of job 530047. Convergence runs at $d\xi/2$ and $2n$. The reduction is `g4c.run` of SOTHE-P2 (window rates and drifts on $[0,5],[5,10],[10,20],[20,30],[30,40]$; loss over $[0,40]$; loss over one Hawking period after $\xi=10$ with $\kappa=NS(1)$ at $\kappa_g=0.3$; recoil check).

## B4. Comoving channel census and Hawking window

`kit.channel_census` (real roots of $\omega=Vk+\Omega_\pm(k)$, low-$k$ branch $|k|<1/3\delta_3$) in the laboratory frame and in the frame comoving with each measured drift $V_d$ of B3: $c'=|v_{g0}(k_{\rm op})-V_d|\,S$, $V=2c'$, closed-form pair edge $\omega^{(-)}_{\max}$ ($V^2/2-\mu$ at $\delta_3=0$). Hawking window (REVIEW_G5 S1-1):

$$\mathcal W=\left[\max(\mu,0.05),\ \min(\omega^{(-)}_{\max},1.6)\right],\qquad |\beta_\mu|^2=\frac{1}{e^{\pi N/S(k_{\rm op})}-1},\qquad E_N(\mu)=2\,\mathrm{arcsinh}|\beta_\mu|.$$

Escape inside the band requires $N<\sqrt{2\cdot1.6}=1.789$.

## B2. Entanglement layer

At the extracted $\kappa$ of each grid point: $|\beta|^2=1/(e^{2\pi\omega/\kappa}-1)$, $E_N=2\,\mathrm{arcsinh}|\beta|$, $\nu_-=e^{-E_N}$, $\bar n_{\max}=\tfrac12(e^{E_N}-1)$, on the band and at the window edges. Symmetric pure loss $\eta$: $\nu_-'=1-\eta(1-e^{-2r})$, $E_N'=-\ln\nu_-'$.

## B5. BdG about the stationary TOD soliton

`kit.tod_soliton` (gauge-fixed Newton, periodic 4th-order stencils, $n_\tau=500$ on $[-16,16]$) and `kit.bdg_spectrum_bg`:

$$M=\begin{pmatrix}H+T&-\phi^2\\ \bar\phi^2&-H+T\end{pmatrix},\qquad H=-\tfrac12D_2-2|\phi|^2+\mu,\qquad T=-i\delta_3D_3 .$$

Reported: Newton residual, tail $\max_{|\tau|>0.8\tau_{\max}}|\phi|/N$, annulus $\max|{\rm Im}\,\omega|$ over $0.1\mu<|{\rm Re}\,\omega|\le3\mu$, Krein census, zero sector (four smallest $|\omega|$), $\sqrt{\epsilon\rho}$, relative singular values of $M$ and $M^2$ with $\dim\ker M$ and $\dim\ker M^2$ (threshold $10^{-12}$), and the $\psi_0$-linearization artefact (the v0.1.0 rate table, ${\rm Im}\,\omega\propto\sqrt{\delta_3}$). The manuscript's BdG scope is $N\delta_3\le0.05$; larger $N\delta_3$ is reported, not claimed.

## B6. Raman and self-steepening (new numerics)

$$\partial_\xi\psi=i\left(\tfrac12\psi_{\tau\tau}+i\delta_3\psi_{\tau\tau\tau}\right)+i\left(1-is\,\partial_\tau\right)\left[\psi\left((1-f_R)|\psi|^2+f_R\,(h\star|\psi|^2)\right)\right],\qquad (h\star I)(\tau)=\int_0^\infty h(\sigma)I(\tau+\sigma)\,d\sigma,$$

integrated by RK4IP (Hult 2007) on the grid of B3, with the field kept inside the window by exact integer-point rolls when the peak moves more than $0.1L$ from the centre. Models: `none` ($f_R=s=0$); `bw` (Blow and Wood: $f_R=0.18$, $\tau_1=12.2$ fs, $\tau_2=32$ fs, analytic transfer function $\tilde h(\omega)=\frac{(\tau_1^2+\tau_2^2)/(\tau_1^2\tau_2^2)}{(1/\tau_2-i\omega)^2+1/\tau_1^2}$ evaluated at $\omega=k/T_0$, first moment $M_1=2\tau_1^2\tau_2/(\tau_1^2+\tau_2^2)=8.122$ fs, $T_R=f_RM_1=1.462$ fs); `lin3` (linear $i\tau_R\psi\,\partial_\tau|\psi|^2$ with $T_R=3$ fs, the manuscript's estimate); `ss` ($s=\lambda_0/(2\pi cT_0)$ at 1550 nm); `bw_ss`. $T_0\in\{50,100,200,500,1000\}$ fs at $(N,\delta_3)=(2,0.05)$ and $(1,0.05)$.

Checks: at $N=1$, $\delta_3=0$ the full-field centroid rate over $\xi\in[1,5]$ equals Gordon's $dk/d\xi=-\tfrac{8}{15}\tau_RN^4$ (linear $\tau_R=0.005,0.01,0.02$; Blow–Wood at $T_0=1000$ fs); $f_R=s=0$ reproduces `kit.splitstep_loss`. Reported: initial rate ($\xi\in[0,0.5]$), mean rate on $[0,5]$, shift over the first Hawking period and over $\xi=40$, window drifts, and the comoving census along each run.

## B7. Mode-conversion scattering (new numerics; D-06)

Frequency-domain BdG operator on the soliton-plus-flow background at $\delta_3=0$ (the FINDINGS configuration):

$$\begin{pmatrix}H+F&-P\\ P&-H+F\end{pmatrix}\begin{pmatrix}u\\ v\end{pmatrix}=\omega\begin{pmatrix}u\\ v\end{pmatrix},\quad H=-\tfrac12\partial_\tau^2+\mu-2P,\quad P=N^2\mathrm{sech}^2(N\tau),\quad F=-\tfrac i2\left(V\partial_\tau+\partial_\tau V\right),$$

$V(\tau)=c\,[1+\tanh(\kappa_0\tau/c)]$ with $\kappa_0=NS(k_{\rm op};0,0.3)$, $c=k_{\rm op}S$. Uniform grid on $[-L,L]$ ($L=25$, $h=0.01$), second-order stencils; outside, the coefficients are constant ($V=0$, $P=0$ on the left; $V=2c$, $P=0$ on the right) and the solution is a sum of the exact discrete exterior modes $z^j$ (roots of the $u$- and $v$-quadratics). Admissible exterior modes are outgoing (discrete group velocity $v_g=\pm\sin(kh)/h+V\cos(kh)$, upper sign for $u$, lower for $v$, pointing away from the grid) or decaying; the square sparse system is factorized once per $\omega$ (SuperLU) and solved for each incoming channel with unit flux. Unitarity is the Krein-signed flux balance

$$\sum_{\rm out}s_j|S_{j\leftarrow i}|^2=s_i,\qquad s=+1\ (u),\ -1\ (v).$$

Reported: channel sets, $|S|^2$ by label (side L/R, component u/v), the conversion ratio $r=|S_{Rv\leftarrow Ru}|^2/|S_{Ru\leftarrow Ru}|^2$, the upstream population $n_H=|S_{Lu\leftarrow Rv}|^2$ where the upstream channel is open ($\omega>\mu$), the Boltzmann and Planck forms at $\kappa_0$, a fit of $\ln n_H$ against $\omega$. Checks: Kaup (N = 1, V = 0: reflectionless, $|R|^2\to0$ as $h^4$), the no-soliton null ($P=0$: no $u\leftrightarrow v$ conversion, exactly zero), convergence under $h\to h/2$ and $L\to35$, the FINDINGS table at $(N,k_{\rm op})=(2,3)$.

## B8. Ref. [25] robustness sweep variants

The shipped script (`run_robustness_sweep.m`, ported statement by statement in `sothe_p2/ref25/src.py`): $N\in\{2.5,\dots,4.5\}$, $\delta_3=$ linspace(0.01, 0.1, 5), pulses sech / Gaussian / super-Gaussian $Ne^{-\tau^4/4}$ / chirped sech ($C=0.5$), centroid-tracking mask at width/order 2.5/12, 3/10, 4/8; $\Delta S_{\rm tot}$ and $\eta_{\rm GSL}$ at $\xi_{\max}$.

| Variant | Change against the script |
|---|---|
| V1_shipped | none ($N_t=2^{13}$, 4000 steps, $\xi_{\max}=10$); control against job 530047 |
| V2_prod | $N_t=2^{14}$, 6000 steps (the production resolution of Ref. [25]) |
| V3_sg8_shipped / V3_sg8_prod | super-Gaussian of the version of record, $Ne^{-\tau^8/2}$ ("order $m=4$"), only |
| V4_d3draft_prod | $\delta_3\in\{0.02,0.04,0.06,0.08,0.10\}$ (the Paper-1 draft text) |
| V5_xi12_prod | $\xi_{\max}=12$ |

The "static" and "comoving-envelope" masks of the version of record are not in the package, so no variant can test them.

## B10. Anchor

At $(N,k_{\rm op},\delta_3,\kappa_g)=(1,1.5,0.05,0.3)$ (REVIEW_G5 S1-1 option A): the four paper2-format CSV files of the thermal and entanglement layers (`thermality_ratio.csv`, `kappa_drive_independence.csv`, `kinematic_kappa_flow.csv`, `partner_log_negativity.csv`, written by `sothe_p2.phase2.write_csv`), and `anchor_summary.json` with the kinematic values, the census per drift epoch, the loss at $(1,0.05)$, the BdG block at $(1,0.05)$ and the $\delta_3=0$ scattering at $(1,1.5)$.
