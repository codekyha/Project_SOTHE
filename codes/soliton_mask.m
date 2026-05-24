function [mask_h, mask_r, k_c] = soliton_mask(k, P, width, order)
% SOLITON_MASK  Adaptive super-Gaussian spectral mask following the soliton.
%
%   [mask_h, mask_r, k_c] = soliton_mask(k, P, width, order)
%
%   Returns a smooth super-Gaussian mask centred on the soliton's
%   instantaneous spectral centroid k_c, plus its complement for the
%   radiation domain.
%
%   Inputs
%   ------
%   k      : column vector, fftshifted spectral grid.
%   P      : |hat{psi}(k)|^2, same shape as k.
%   width  : half-width of the mask (in same units as k). Default 3.0.
%   order  : super-Gaussian order (even integer >= 2). Default 10
%            (a near-rectangular roll-off, smooth enough that aliasing in
%            entropy estimates stays below 1e-3 nats).
%
%   Outputs
%   -------
%   mask_h : soliton (horizon) mask, in (0,1].
%   mask_r : radiation mask, 1 - mask_h.
%   k_c    : spectral centroid (intensity-weighted first moment) at which
%            mask_h is centred. Equal to 0 when P is symmetric.
%
%   Why a dynamic centroid?
%   -----------------------
%   The original (static-mask) approach fails once Raman or self-steepening
%   shifts the soliton off the input frequency. Even in the truncated GNLSE
%   used here, the centroid drifts by O(1) over xi in [0, 12]. Tracking
%   the centroid is the "dynamic" scheme in summary.txt (Run 1).

    if nargin < 3, width = 3.0; end
    if nargin < 4, order = 10;  end
    if mod(order, 2) ~= 0
        error('soliton_mask:order', 'order must be an even positive integer.');
    end

    k = k(:);
    P = P(:);

    Ptot = trapz(k, P);
    if Ptot < 1e-18
        k_c = 0;
    else
        k_c = trapz(k, k .* P) / Ptot;
    end

    mask_h = exp(-((k - k_c) / width).^order);
    mask_r = 1 - mask_h;
end
