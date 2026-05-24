function out = gnlse_dimensionless(varargin)
% GNLSE_DIMENSIONLESS  Symmetric split-step Fourier integrator for the
% dimensionless higher-order NLSE used in the SOTHE paper.
%
%   out = gnlse_dimensionless('Name', Value, ...)
%
%   Integrates
%
%       \partial_\xi \psi = i ( 1/2 \partial_\tau^2 \psi + |\psi|^2 \psi
%                                              + i \delta_3 \partial_\tau^3 \psi )
%
%   on a periodic tau-grid, returning the on-grid history (psi(tau, xi))
%   plus the bipartite spectral entropies S_hor(xi), S_rad(xi) and the
%   coarse-grained GSL efficiency eta_GSL(xi).
%
%   Spectral operator under the MATLAB fft convention (kernel e^{-ikt}):
%
%       D(k) = i ( -1/2 k^2  +  delta3 k^3 )
%
%   ** TOD SIGN NOTE **
%   The original 2025 submission code used D(k) = i (1/2 k^2 - delta3 k^3),
%   which corresponds to -i delta3 partial_tau^3 psi -- the wrong sign for
%   the integrability-consistent anomalous-dispersion convention. The
%   revision flips the sign, and so does this routine. See NOTES.md.
%
%   Name-Value pairs (defaults in brackets)
%   ---------------------------------------
%     'N_sol'    soliton order N                                  [3.5]
%     'delta3'   dimensionless third-order dispersion             [0.02]
%     'xi_max'   propagation distance in L_D                      [12.0]
%     'Nt'       temporal grid size, power of 2                   [2^14]
%     'TauW'     temporal window in T_0                           [20.0]
%     'n_steps'  number of SSFM steps                             [6000]
%     'n_save'   number of snapshots stored                       [200]
%     'mask_w'   soliton-mask half-width (in 1/T_0)               [3.0]
%     'mask_n'   super-Gaussian order                             [10]
%     'verbose'  print progress                                   [false]
%
%   Returns a struct with fields:
%     tau, xi, k, psi_hist, S_hor, S_rad, S_tot, eta, k_c_hist,
%     P_norm, params.

    p = inputParser;
    p.addParameter('N_sol',   3.5,  @isnumeric);
    p.addParameter('delta3',  0.02, @isnumeric);
    p.addParameter('xi_max',  12.0, @isnumeric);
    p.addParameter('Nt',      2^14, @(x) (x==2^round(log2(x))));
    p.addParameter('TauW',    20.0, @isnumeric);
    p.addParameter('n_steps', 6000, @isnumeric);
    p.addParameter('n_save',  200,  @isnumeric);
    p.addParameter('mask_w',  3.0,  @isnumeric);
    p.addParameter('mask_n',  10,   @isnumeric);
    p.addParameter('verbose', false,@islogical);
    p.parse(varargin{:});
    P  = p.Results;

    % ---- grid ----
    Nt   = P.Nt;
    dtau = P.TauW / Nt;
    tau  = (-Nt/2 : Nt/2 - 1).' * dtau;

    dk = 2*pi / P.TauW;
    k  = (-Nt/2 : Nt/2 - 1).' * dk;
    k  = ifftshift(k);                  % fft-aligned for the kernel
    k_shift = fftshift(k);              % display / entropy axis

    % ---- dispersion operator (CORRECTED SIGN) ----
    D = 1i * (-0.5 * k.^2 + P.delta3 * k.^3);

    % ---- initial condition: exact fundamental soliton scaled by N ----
    psi = P.N_sol * sech(tau);
    psi_f = fft(psi);

    % ---- propagation bookkeeping ----
    dxi      = P.xi_max / P.n_steps;
    save_id  = round(linspace(1, P.n_steps, P.n_save));
    xi       = linspace(0, P.xi_max, P.n_save).';

    psi_hist = zeros(Nt, P.n_save);
    S_hor    = zeros(P.n_save, 1);
    S_rad    = zeros(P.n_save, 1);
    k_c_hist = zeros(P.n_save, 1);
    P_norm   = zeros(P.n_save, 1);

    if P.verbose, fprintf('gnlse_dimensionless: %d steps, dxi = %.3e\n', P.n_steps, dxi); end

    c = 1;
    half = exp(D * dxi/2);
    for n = 1:P.n_steps
        % Symmetric SSFM: D/2 -- N -- D/2
        psi_f = psi_f .* half;
        psi_t = ifft(psi_f);
        psi_t = psi_t .* exp(1i * abs(psi_t).^2 * dxi);
        psi_f = fft(psi_t);
        psi_f = psi_f .* half;

        if n == save_id(c)
            psi_hist(:, c) = ifft(psi_f);

            Pspec = abs(fftshift(psi_f)).^2;
            [mh, mr, k_c] = soliton_mask(k_shift, Pspec, P.mask_w, P.mask_n);

            S_hor(c)    = spectral_entropy(k_shift, Pspec, mh);
            S_rad(c)    = spectral_entropy(k_shift, Pspec, mr);
            k_c_hist(c) = k_c;
            P_norm(c)   = trapz(k_shift, Pspec);   % photon-number proxy

            c = c + 1;
            if P.verbose && mod(c, 20) == 0
                fprintf('  step %d/%d, xi = %.2f, S_tot = %.3f nats\n', ...
                        n, P.n_steps, xi(c-1), S_hor(c-1) + S_rad(c-1));
            end
        end
    end

    S_tot = S_hor + S_rad;
    eta   = eta_gsl(xi, S_tot);

    out = struct( ...
        'tau',      tau, ...
        'xi',       xi, ...
        'k',        k_shift, ...
        'psi_hist', psi_hist, ...
        'S_hor',    S_hor, ...
        'S_rad',    S_rad, ...
        'S_tot',    S_tot, ...
        'eta',      eta, ...
        'k_c_hist', k_c_hist, ...
        'P_norm',   P_norm, ...
        'params',   P);
end
