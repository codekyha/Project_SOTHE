function res = rr_phase_match_fit(delta3_vec, k_RR_vec)
% RR_PHASE_MATCH_FIT  Akhmediev-Karlsson Cherenkov scaling fit.
%
%   res = rr_phase_match_fit(delta3_vec, k_RR_vec)
%
%   Fits the measured resonant-radiation wavenumbers k_RR to the
%   integrable-perturbation phase-matching law
%
%       k_RR(delta3) = c_RR / delta3,
%
%   equivalently k_RR * delta3 = c_RR. The fit is through the origin in
%   the (1/delta3, k_RR) plane and so contains a single parameter c_RR.
%   This is Eq. (15) of the revised paper and replaces the (withdrawn)
%   Boltzmann single-slope fit that the original 2025 submission used to
%   extract an effective T_H.
%
%   The paper reports c_RR = 1.01, mean k_RR*delta3 = 1.02, R^2 = 0.93
%   on a 5-point scan over delta3 in [0.06, 0.10] at N = 3.5.
%
%   Inputs
%   ------
%   delta3_vec : row or column vector of dimensionless TOD strengths.
%   k_RR_vec   : measured RR peak positions, same length.
%
%   Output fields
%   -------------
%   c_RR     : fitted constant (intercept-free LS).
%   product  : column vector of k_RR_i * delta3_i values.
%   R2       : coefficient of determination for the k_RR vs c_RR/delta3 fit.
%   delta3   : input delta3_vec (column).
%   k_RR     : input k_RR_vec (column).

    delta3 = delta3_vec(:);
    kRR    = k_RR_vec(:);

    % Through-origin least squares: c_RR = sum(x*y)/sum(x*x), x = 1/delta3
    x     = 1 ./ delta3;
    c_RR  = sum(x .* kRR) / sum(x .* x);

    yhat  = c_RR ./ delta3;
    SSres = sum((kRR - yhat).^2);
    SStot = sum((kRR - mean(kRR)).^2);
    R2    = 1 - SSres / max(SStot, eps);

    res.c_RR    = c_RR;
    res.product = kRR .* delta3;
    res.R2      = R2;
    res.delta3  = delta3;
    res.k_RR    = kRR;
end
