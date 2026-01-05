%% Optical Entropy and Soliton Fission Simulation
%  Solves GNLSE with Beta2, Beta3 and Kerr Nonlinearity
%  Computes Spectral Shannon Entropy (Radiation) vs Horizon Entropy (Soliton)
%  Author: H. Oguz & A. I. Collaborator

clear; clc; close all;

%% 1. Physics Parameters
c_light = 299792.458;       % Speed of light [nm/ps]
lambda0 = 800;              % Pump wavelength [nm]
omega0  = 2*pi*c_light/lambda0; % Angular frequency [rad/ps]

% Fiber Parameters (PCF)
beta2   = -15.0;            % GVD [ps^2/km] (Anomalous)
beta3   = +0.1;             % TOD [ps^3/km] (Driver of radiation)
gamma   = 0.1;              % Nonlinearity [1/W/km]

% Soliton Parameters
N_sol   = 3.0;              % Soliton order (N=3 ensures fission)
T0      = 0.050;            % Pulse width [ps] (50 fs)

%% 2. Grid and Computational Domain
n_points = 2^14;            % Grid size (16384 pts)
T_window = 10.0;            % Time window [ps]

dt = T_window / n_points;
t  = (-n_points/2 : n_points/2 - 1) * dt; % Time vector

% Frequency Grid (shifted for FFT)
dw = 2*pi / T_window;
w  = (0:n_points-1)' * dw;
w  = w - n_points*dw/2;     % Center around 0
w  = fftshift(w);           % Shift for FFT operations

%% 3. Initial Conditions
LD = (T0^2) / abs(beta2) * 1e-3;    % Dispersion length [km]
P0 = (N_sol^2 * abs(beta2)) / (gamma * T0^2); % Peak Power [W]

% Initial Sech Pulse
U = sqrt(P0) * sech(t / T0)';

% Dispersion Operator (Linear Step)
% D = i * (beta2/2 * w^2 + beta3/6 * w^3) - alpha/2
D_op = 1i * (beta2/2 * w.^2 + beta3/6 * w.^3);

%% 4. Spectral Filter Definition (The Thermodynamic Partition)
% We define a Super-Gaussian filter to separate the Soliton (Low Freq)
% from the Resonant Radiation (High Freq/Dispersive).
filter_width = 4.0 * (1/T0); % Width in frequency domain
mask_sol     = exp( -((w)/(filter_width)).^8 ); 
mask_rad     = 1 - mask_sol;

% Plot filter to verify
figure(1); plot(fftshift(w), fftshift(mask_sol)); 
title('Spectral Partition Filter'); xlabel('\omega'); ylabel('Transmission');
drawnow;

%% 5. SSFM Main Loop
L_prop  = 10 * LD;          % Propagation distance [km]
n_steps = 4000;             % Number of z-steps
dz      = L_prop / n_steps;

% Pre-compute Operators
half_step_disp = exp(D_op * dz / 2);

% Storage for Plotting
n_saves = 200;              % How many snapshots to save
save_idx = round(linspace(1, n_steps, n_saves));
U_save   = zeros(n_points, n_saves);
z_save   = linspace(0, L_prop, n_saves);

% Entropy Storage
S_rad_hist = zeros(1, n_saves);
S_hor_hist = zeros(1, n_saves);

counter = 1;
fprintf('Starting Simulation (L = %.4f km, N = %.1f)...\n', L_prop, N_sol);

for k = 1:n_steps
    
    % --- Step 1: Dispersion (Half) ---
    U_f = fft(U);
    U_f = U_f .* half_step_disp;
    U   = ifft(U_f);
    
    % --- Step 2: Nonlinearity (Full) ---
    U   = U .* exp(1i * gamma * abs(U).^2 * dz);
    
    % --- Step 3: Dispersion (Half) ---
    U_f = fft(U);
    U_f = U_f .* half_step_disp;
    U   = ifft(U_f);
    
    % --- Data Extraction ---
    if ismember(k, save_idx)
        U_save(:, counter) = U;
        
        % --- ENTROPY CALCULATION ---
        % 1. Power Spectrum
        Spec = abs(fftshift(fft(U))).^2;
        W_shifted = fftshift(w); % Axis for integration
        
        % 2. Horizon Entropy (Proportional to trapped energy/area)
        % We integrate the spectrum INSIDE the soliton mask
        Mask_shifted = fftshift(mask_sol);
        E_sol = trapz(W_shifted, Spec .* Mask_shifted);
        S_hor_hist(counter) = E_sol; % Proportional to N_sol
        
        % 3. Radiation Entropy (Shannon Spectral Entropy)
        % We look at the spectrum OUTSIDE the soliton mask
        Spec_rad = Spec .* (1 - Mask_shifted);
        
        % Normalize to treat as probability distribution
        total_rad_power = trapz(W_shifted, Spec_rad);
        
        if total_rad_power > 1e-20
            p_rad = Spec_rad / total_rad_power;
            
            % Compute -Sum(p ln p) avoiding log(0)
            valid_idx = p_rad > 1e-15;
            S_rad_val = -sum(p_rad(valid_idx) .* log(p_rad(valid_idx)));
            S_rad_hist(counter) = S_rad_val;
        else
            S_rad_hist(counter) = 0;
        end
        
        counter = counter + 1;
    end
end

fprintf('Simulation Complete.\n');

%% 6. Visualization and Analysis

% --- Prepare Data for Plotting ---
% Normalize Entropies for comparison (Arbitrary Units)
% We want to show the TREND, not absolute values (which depend on units)
S_rad_norm = S_rad_hist / max(S_rad_hist);
S_hor_norm = S_hor_hist / max(S_hor_hist);

% Total Entropy Proxy: 
% We weight them to show the conservation/growth principle.
% In a rigorous thermodynamic treatment, coefficients depend on T_Hawking.
% Here we demonstrate the phenomenological trade-off.
S_total = S_rad_norm + S_hor_norm; 

% Frequency axis for plotting (convert to THz)
freq_THz = fftshift(w) / (2*pi); 

% Convert field to dB
Spectrogram_dB = 10*log10(abs(fftshift(fft(U_save, [], 1), 1)).^2 + 1e-10);

% --- Figure generation ---
figure('Position', [100, 100, 1000, 400], 'Color', 'w');

% Subplot 1: Spectral Evolution
subplot(1, 2, 1);
imagesc(z_save, freq_THz, Spectrogram_dB);
axis xy;
colormap(parula);
c = colorbar; c.Label.String = 'Power Spectral Density (dB)';
ylim([-100, 100]); % Zoom in on relevant range
xlabel('Propagation Distance z (km)');
ylabel('Frequency Detuning \nu (THz)');
title('(a) Spectral Evolution (Soliton Fission)');

% Subplot 2: Entropy Dynamics
subplot(1, 2, 2);
hold on;
plot(z_save, S_hor_norm, 'r--', 'LineWidth', 2);
plot(z_save, S_rad_norm, 'b-', 'LineWidth', 2);
plot(z_save, S_total, 'k-', 'LineWidth', 1.5);
hold off;

grid on;
legend('Horizon Area (S_{BH})', 'Radiation Entropy (S_{rad})', 'Total Entropy', ...
       'Location', 'best');
xlabel('Propagation Distance z (km)');
ylabel('Normalized Entropy (a.u.)');
title('(b) Generalized Second Law');
xlim([0, L_prop]);

% Save high-res image for paper
print('Entropy_Analysis_Figure','-dpng','-r300');
fprintf('Figure saved as Entropy_Analysis_Figure.png\n');