function S = spectral_entropy(k, P, mask)
% SPECTRAL_ENTROPY  Shannon entropy of a normalized spectral distribution.
%
%   S = spectral_entropy(k, P, mask) computes
%
%       S = - \int_{Omega} p(k) ln p(k) dk,
%
%   where p(k) = P(k) * mask(k) / \int P(k) * mask(k) dk is the normalized
%   power spectral density restricted to the support specified by mask.
%
%   Inputs
%   ------
%   k     : column vector, fftshifted angular-frequency grid (rad/T0 in
%           the dimensionless GNLSE, or rad/ps in physical units).
%   P     : |hat{psi}(k)|^2, same shape as k.
%   mask  : weighting in [0,1] selecting the subsystem (1 inside, 0 outside,
%           with a smooth super-Gaussian edge in practice).
%
%   Output
%   ------
%   S     : Shannon entropy in nats. Returns 0 if the masked total power is
%           below 1e-18 (numerically empty subsystem).
%
%   Notes
%   -----
%   This is Eq. (2)/(5) of:
%       H. Oguz, "Generalized Thermodynamics of Solitonic Event Horizons
%       in Dispersive Field Theories", CQG-114816 (revised, 2026).
%
%   Both S_hor and S_rad in the paper are computed by this routine; the
%   only difference between them is the mask (soliton-following vs
%   radiation-domain). The exponential exp(S) is the spectral
%   participation number N_S used as an interpretive aid (Sec. II.B).

    if nargin < 3
        error('spectral_entropy:argChk', 'Three arguments required: k, P, mask.');
    end

    k    = k(:);
    P    = P(:);
    mask = mask(:);

    Pw  = P .* mask;
    nrm = trapz(k, Pw);

    if nrm < 1e-18
        S = 0;
        return
    end

    p   = Pw ./ nrm;
    idx = p > 1e-20;        % drop machine-noise floor to keep log finite
    S   = -trapz(k(idx), p(idx) .* log(p(idx)));
end
