%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% SOTHE: Solitonic Optical Thermodynamics of Event Horizons
% Full SSFM simulation + entropy diagnostics
%
% Author: H. Oguz
% Purpose: Numerical validation of generalized second law for optical horizons
%
% Governing equation:
% dA/dz = -i beta2/2 d^2A/dT^2 + beta3/6 d^3A/dT^3 + i gamma |A|^2 A
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

clear; clc;

%% ---------------- Physical constants ----------------
c     = 2.99792458e8;        % m/s
hbar  = 1.054571817e-34;     % J*s
kB    = 1.380649e-23;        % J/K

%% ---------------- Fiber parameters ------------------
lambda0 = 800e-9;                    % m
omega0  = 2*pi*c/lambda0;            % rad/s

beta2 = -15e-27;                     % s^2/m
beta3 =  0.1e-39;                    % s^3/m
gamma =  1.2e-3;                     % 1/(W m)

%% ---------------- Soliton parameters ----------------
T0 = 50e-15;                         % s
P0 = abs(beta2)/(gamma*T0^2);        % Fundamental soliton power
Nsoliton = 1;                        % Soliton order (1 or 3)

P0 = Nsoliton^2 * P0;

%% ---------------- Temporal grid ---------------------
NT   = 2^14;                         % Grid points
Tmax = 10*T0;                        % Window half-width
dt   = 2*Tmax/NT;
t    = (-NT/2:NT/2-1)*dt;

%% ---------------- Frequency grid --------------------
dw = 2*pi/(NT*dt);
omega = fftshift((-NT/2:NT/2-1))*dw;

%% ---------------- Initial condition -----------------
A = sqrt(P0)*sech(t/T0);

%% ---------------- Propagation parameters ------------
LD  = T0^2/abs(beta2);
LNL = 1/(gamma*P0);
dz  = min([LD, LNL])/50;             % Stability-controlled step
Nz  = 3000;                          % Number of steps
z   = (0:Nz-1)*dz;

%% ---------------- Linear operator -------------------
Lw = 1i*( beta2/2*omega.^2 - beta3/6*omega.^3 );

%% ---------------- Spectral filters -----------------
DeltaOmegaS = 5e13;                  % Soliton bandwidth
WS = exp(-((omega-omega0)/DeltaOmegaS).^8);
WR = 1 - WS;

%% ---------------- Storage arrays -------------------
Srad = zeros(1,Nz);
Shor = zeros(1,Nz);
Stot = zeros(1,Nz);

EnergySol = zeros(1,Nz);
EnergyRad = zeros(1,Nz);

%% ---------------- Main SSFM loop --------------------
for n = 1:Nz

    % Linear half-step
    A = ifft( exp(Lw*dz/2) .* fft(A) );

    % Nonlinear full-step
    A = A .* exp(1i*gamma*abs(A).^2*dz);

    % Linear half-step
    A = ifft( exp(Lw*dz/2) .* fft(A) );

    % Spectrum
    Aw = fftshift(fft(A));
    Iw = abs(Aw).^2;

    % ---------------- Radiation entropy ----------------
    Prad = WR .* Iw;
    Prad = Prad / trapz(omega,Prad);

    Srad(n) = -kB * trapz(omega, Prad .* log(Prad + 1e-20));

    % ---------------- Horizon entropy ------------------
    Nsol = trapz(omega, WS .* Iw ./ (hbar*abs(omega)));
    eta  = 1;   % dimensionless geometric factor
    Shor(n) = eta * Nsol;

    % ---------------- Energies (diagnostic) ------------
    EnergySol(n) = trapz(omega, WS .* Iw);
    EnergyRad(n) = trapz(omega, WR .* Iw);

    % ---------------- Total entropy --------------------
    Stot(n) = Srad(n) + Shor(n);

    % Progress update
    if mod(n,200)==0
        fprintf('Step %d / %d completed\n',n,Nz);
    end
end

%% ---------------- Entropy production rate -----------
dStot_dz = diff(Stot)/dz;

%% ===================== VISUALIZATION =====================

%% Figure 1: Temporal evolution
figure;
imagesc(z*1e3, t*1e15, abs(A).^2);
axis xy;
xlabel('z (mm)');
ylabel('T (fs)');
title('Temporal Intensity Evolution');
colorbar;

%% Figure 2: Spectral evolution (final)
figure;
plot((omega-omega0)/2/pi/1e12, Iw/max(Iw),'k','LineWidth',1.5);
xlabel('Frequency detuning (THz)');
ylabel('Normalized spectral intensity');
title('Final Optical Spectrum');
grid on;

%% Figure 3: Entropy evolution
figure;
plot(z*1e3, Srad,'r','LineWidth',1.5); hold on;
plot(z*1e3, Shor,'b','LineWidth',1.5);
plot(z*1e3, Stot,'k','LineWidth',2);
xlabel('z (mm)');
ylabel('Entropy (arb. units)');
legend('Radiation','Horizon','Total','Location','best');
title('Entropy Flow and Generalized Second Law');
grid on;

%% Figure 4: GSL diagnostic
figure;
plot(z(1:end-1)*1e3, dStot_dz,'k','LineWidth',1.5);
xlabel('z (mm)');
ylabel('dS_{tot}/dz');
title('Generalized Second Law Diagnostic');
yline(0,'--');
grid on;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% END OF FILE
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
