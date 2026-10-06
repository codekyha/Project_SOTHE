# SOTHE-P2-RAD: methods

Mathematics is LaTeX. Equation numbers of the article are not used here; the article cites this package and its results
archive for every number that it quotes from them.

## 1. Equation and conventions

$$\partial_\xi\psi=i\left(\tfrac12\partial_\tau^2\psi+|\psi|^2\psi+i\delta_3\partial_\tau^3\psi\right),$$

the dimensionless generalized nonlinear Schrödinger equation (GNLSE) of the article, with $\xi=z/L_D$, $\tau=-T/T_0$ and
$\delta_3=\beta_3/(6|\beta_2|T_0)$. A plane wave $e^{ik\tau-i\beta_0(k)\xi}$ has $\beta_0(k)=\tfrac12k^2-\delta_3k^3$ and group
velocity $v_{g0}(k)=k-3\delta_3k^2$; $k>0$ is the blue side. The fundamental soliton of the integrable equation is
$\psi_0=N\operatorname{sech}(N\tau)$ with propagation constant $\mu=N^2/2$ and norm $Q=2N$.

## 2. Numerical scheme

Symmetric split step on $n=8192$ points over $\tau\in[-200,200)$, $\Delta\xi=0.002$: half a linear step
$\exp[i(-\tfrac12k^2+\delta_3k^3+Vk)\Delta\xi/2]$ in Fourier space, the Kerr step $\psi\to\psi e^{i|\psi|^2\Delta\xi}$, half a
linear step, then the absorber $e^{-\sigma\Delta\xi}$ with $\sigma=8[(|\tau|-a)/(L/2-a)]^2$ for $|\tau|>a=0.35L$. With $V=0$ and a
constant $\delta_3$ this is the scheme, grid and absorber of the radiative-loss runs of SOTHE-P2 1.0.0
(`kit.splitstep_loss`), whose recorded loss at $(N,\delta_3)=(1,0.05)$ the unit tests reproduce to $10^{-6}$ in relative terms.
$V\neq0$ integrates the equation in the frame $\tau'=\tau-V\xi$ (a translation; the unit tests check that the spectra of the
two frames agree to round-off). A prepared soliton is the integrable soliton launched at $\delta_3=0$ while $\delta_3$ rises as
$\delta_3(\xi)=\delta_3\,[1-\cos(\pi\xi/\xi_r)]/2$ over $\xi_r=100$ and is then held; $\delta_3$ is evaluated at the midpoint of
each step.

## 3. Local soliton and its scaled third-order coefficient

About its mean wavenumber $k_s$ (momentum centroid), the field $\psi=e^{ik_s\tau}f$ obeys
$$\partial_\xi f=-i\beta_0(k_s)f-v_{g0}(k_s)\partial_\tau f+i\left(\tfrac12D_2\partial_\tau^2f+|f|^2f\right)-\delta_3\partial_\tau^3f,\qquad D_2=1-6\delta_3k_s ,$$
exactly (the dispersion is cubic). Rescaling $\tau$ by $\sqrt{D_2}$ and by the amplitude $N$ maps the local soliton
$f=N\operatorname{sech}(N\tau/\sqrt{D_2})$ onto the zero-carrier soliton of amplitude 1 of the original equation with
$$\tilde\delta=\frac{N\delta_3}{D_2^{3/2}} ,\qquad Q=2N\sqrt{D_2},\qquad \xi\to N^2\xi .$$
The centre of mass moves at the spectrally averaged group velocity (exactly, for this Hamiltonian equation), which for the
local sech is $V=v_{g0}(k_s)-\delta_3N^2/D_2$; at $k_s=0$ it is the drift $-\delta_3N^2$ of the article.

## 4. Resonance

A soliton $a(\tau-V\xi)\,e^{ik_s\tau+i\theta\xi}$ drives the linear waves at the frequency $\theta-(k-k_s)V$; the radiation is
phase-matched where
$$\beta_0(k)+\theta-(k-k_s)V=0 .$$
With $\theta=-\beta_0(k_s)+N^2/2$ and the velocity of section 3 this becomes, in the scaled wavenumber
$q=(k-k_s)\sqrt{D_2}/N$,
$$\tilde\delta q^3-\tfrac12q^2-\tilde\delta q-\tfrac12=0\qquad\text{(drift-corrected)},$$
whose largest root $q_r(\tilde\delta)$ is the resonant (Cherenkov) wavenumber; without the drift term,
$\tilde\delta q^3-\tfrac12q^2-\tfrac12=0$. The measured version uses the phase rotation rate $\varphi_c$ at the soliton centre and
the measured velocity: $\delta_3k^3-\tfrac12k^2+Vk-\varphi_c=0$. The wavenumber of the radiation found behind the soliton
(section 6) is compared with both.

## 5. First-order (Born) rate

The third-order term acting on the sech soliton is the source $S=-\delta_3\psi_0'''e^{i\theta\xi}$, with Fourier amplitude
$|\hat S(k)|=\delta_3|k|^3\,\pi\operatorname{sech}(\pi k/2N)$. Driven at the phase-matched wavenumber, the radiated norm grows
linearly in $\xi$ (golden rule), at the rate $|\hat S(k_r)|^2/|v_{\rm rel}|$ with $v_{\rm rel}$ the group velocity of the radiation
relative to the soliton. Relative to the norm, in scaled units,
$$R_{\rm Born}(\tilde\delta)=\frac{\pi^2\tilde\delta^2q_r^6\operatorname{sech}^2(\pi q_r/2)}{2\,|q_r-3\tilde\delta q_r^2+\tilde\delta|},\qquad R=N^2R_{\rm Born}(\tilde\delta)\ \text{per unit}\ \xi .$$
The exponent $\pi q_r$ is fixed by the pole of the sech nearest the real axis, at $\tau=i\pi/2N$. The prefactor of the first-order
formula is not reliable at the phase-matched wavenumber, where the third-order term is as large as the second-order one; the
measured rates (section 7) fix it.

## 6. Estimators of the emission rate

* **Flux.** Behind the soliton the resonant radiation of a steady emission is a plane wave of uniform intensity $I$. The tail
  window $[\tau_c-120,\tau_c-20]$ (inside the region without absorption) is tapered (Tukey, 20%), filtered to
  $|k-k_c|<2$ around the predicted resonant wavenumber $k_c$, and $I$ is the mean filtered intensity over the flat part, also
  in eight bins along the window. The emitted norm flux is $F=I\,|v_{g0}(k_b)-V|$, with $k_b$ the power-weighted wavenumber in the
  band; radiation at distance $d$ behind the soliton was emitted at the retarded time $\xi_e=\xi-d/|v_{g0}(k_b)-V|$, and
  $R(\xi_e)=F/Q(\xi_e)$. The nearest bin is used. A uniform intensity over the eight bins identifies a steady emission.
* **Core.** $R=-d\ln Q/d\xi$ from a local linear fit over $\xi\pm2.5$, $Q$ the core norm $\int_{|\tau-\tau_m|<6/N}|\psi|^2d\tau$
  about the grid maximum $\tau_m$ (the definition of SOTHE-P2). The internal oscillation of a launched soliton limits it to
  rates above about $10^{-6}$.
* **Box.** $-d\ln Q_{\rm box}/d\xi$ of the norm in the region without absorption; in a steady state it equals the emission rate.
  The norm of the scheme drifts by about one unit in the last place per step (round-off; about $10^{-13}$ per unit $\xi$ at
  $\Delta\xi=0.002$), which limits this estimator to rates above about $10^{-12}$.
* **Momentum.** $k_s$ is the momentum centroid of the core, $\int_{\rm core}\mathrm{Im}(\psi^*\partial_\tau\psi)\,d\tau/Q$ (spectral
  derivative); the Hann-windowed spectral centroid of SOTHE-P2 is recorded as well.

## 7. Scaled rate function

If the emission is controlled by the local soliton, $R/N^2$ depends on the state only through $\tilde\delta$. The hold phase of
a prepared soliton is $\xi\ge\xi_r+40$ (radiation emitted after $\xi_r+5$); a prepared run is *steady* when its flux rate varies
by less than 5 percent (relative standard deviation) over the hold phase, which holds for $\delta_3\le0.07$. At larger $\delta_3$
the soliton loses enough norm during the ramp that $\tilde\delta$ and the rate still change over the hold phase.

*Calibration.* The form
$$R_s(\tilde\delta)=C\,\tilde\delta^{\,p}\,e^{\,g\tilde\delta-\pi q_r(\tilde\delta)}$$
is fitted (weighted least squares in $\ln R_s+\pi q_r$) to every hold-phase rate of the prepared solitons, each with its own
$\tilde\delta$ at the retarded time (from the measured $Q$ and $k_s$), with equal weight per run, so that slowly evolving states
enter pointwise. It is an interpolation formula over its range of $\tilde\delta$: $C$, $p$ and $g$ are strongly correlated and have
no separate meaning. For comparison, two- and three-parameter forms (with the drift-corrected and the no-drift root) are fitted to
the mean rates of the steady runs, and the ratio of the steady rates to the first-order rate and the slope of $\ln R$ against
$1/(N\delta_3)$ are reported.

*Test.* The rates of the launched solitons ($N=1$, $1.6$, $2$), which are not used in the fit, are compared with $R_s$ after the
launch transient (radiation emitted at $\xi_e\ge20$). The slope of the core norm of the prepared runs, where it resolves the rate,
is a cross-check of the flux estimator.

## 8. Recoil model

Norm and momentum balance for a soliton that radiates at $k_r$:
$$\frac{dQ}{d\xi}=-RQ,\qquad \frac{d(k_sQ)}{d\xi}=k_r\frac{dQ}{d\xi}\ \Longrightarrow\ \frac{dk_s}{d\xi}=-(k_r-k_s)R ,$$
closed by the local soliton of section 3: $N=Q/(2\sqrt{D_2})$, $\tilde\delta=N\delta_3/D_2^{3/2}$, $k_r-k_s=q_r(\tilde\delta)N/\sqrt{D_2}$ and
$R=N^2R_s(\tilde\delta)$, with the rate function of section 7 and no other parameter. The model is started from the state of each
launched run at $\xi=10$ and integrated to the end of the run ($\xi=40$ or $160$; DOP853, relative tolerance $10^{-11}$); the norm
lost since $\xi=10$, the total loss, $k_s$ and the velocity $V$ are compared with the run. Part of the launch transient of a
launched soliton is still decaying at $\xi=10$, and the model does not contain it. Two variants are evaluated on the runs to
$\xi=160$, each started at $\xi=40$ with a multiplier fixed by the flux rate emitted at $\xi=40\pm2$: the first-order rate law;
and the amplitude law $dN/d\xi=-NR$ with $R\propto e^{-\pi q_r(\tilde\delta)}$ (no-drift root). The velocity of section 3 is
of first order in $\delta_3$; the measured velocity is larger in magnitude by a few percent.

Asymptotics. With $\dot R=-KR^2$, $R\to1/(K\xi)$ and $-d\ln Q/d\xi\to c/\xi$ with $c=1/K$; $K=d(1/R)/d\xi$ is evaluated along the
model solution and compared with the window exponents $c=-\Delta\ln Q/\Delta\ln\xi$ of the runs over $[10,20]$, $[20,40]$,
$[40,80]$ and $[80,160]$. The momentum balance itself is checked along the runs: the median ratio of the measured $dk_s/d\xi$ to
$-(k_b-k_s)R$ ($k_b$ the band wavenumber, $R$ the core rate).

## 9. Catalogue

`python rad.py list` prints it (`sothe_rad/experiments.py`): the ten launched points of the loss table of the article
(laboratory frame, $\xi\le40$); prepared solitons at $\delta_3=0.035$ to $0.12$ ($N=1$); launched solitons at $\delta_3=0.06$ to
$0.12$ to $\xi=160$; two $N=2$ images of these under the exact scaling $\psi=N\Psi(N\tau,N^2\xi)$; and numerical checks
($\Delta\xi/2$, $2n$, a box of $L=600$, ramps of $50$ and $200$, a spectral filter at $0.8\,k_{\max}$, the laboratory frame).

## 10. Limits

The equation has no Raman term and no self-steepening (the article's Raman-free Kerr medium); the dispersion is truncated
at third order, which matters at the resonant wavenumbers $k\simeq5$ to $14$ (frequency offsets of $16$ to $45$ THz at
$T_0=50$ fs). The recoil model is adiabatic: it assumes that the soliton stays a local sech between emissions and that the
radiation leaves at a single wavenumber. Rates below about $10^{-15}$ per unit $\xi$ are not measured.
