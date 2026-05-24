function [k_RR, amp_RR] = find_rr_lobe(k, Spec, k_c, side, k_offset)
% FIND_RR_LOBE  Locate the resonant-radiation peak position robustly.
%
%   [k_RR, amp_RR] = find_rr_lobe(k, Spec, k_c, side, k_offset)
%
%   Returns the wavenumber and amplitude of the most prominent local
%   maximum of |hat{psi}(k)|^2 in the radiation domain (i.e. away from
%   the soliton centroid k_c by at least k_offset). The "side" argument
%   selects the +k or -k lobe.
%
%   This routine is the operational replacement for the original
%   2025-submission's "fit a single-slope Boltzmann to the tail" step,
%   which is withdrawn in the revision (see NOTES.md, Section 2).
%
%   Inputs
%   ------
%   k        : fftshifted spectral grid (column).
%   Spec     : |hat{psi}(k)|^2 on the same grid (column).
%   k_c      : soliton spectral centroid; pass 0 if unknown.
%   side     : '+' or '-' (which Cherenkov lobe).            [default '+']
%   k_offset : minimum |k - k_c| considered radiation.       [default 2.0]

    if nargin < 4, side     = '+'; end
    if nargin < 5, k_offset = 2.0; end

    k    = k(:);
    Spec = Spec(:);

    switch side
        case '+'
            sel = k > (k_c + k_offset);
        case '-'
            sel = k < (k_c - k_offset);
        otherwise
            error('find_rr_lobe:side', 'side must be ''+'' or ''-''.');
    end

    if ~any(sel)
        k_RR = NaN; amp_RR = 0; return
    end

    Sw     = Spec(sel);
    kw     = k(sel);
    [amp_RR, j] = max(Sw);
    k_RR   = kw(j);

    % Local parabolic refinement (sub-grid)
    if j > 1 && j < numel(Sw)
        y1 = log(Sw(j-1) + eps);
        y2 = log(Sw(j  ) + eps);
        y3 = log(Sw(j+1) + eps);
        denom = (y1 - 2*y2 + y3);
        if abs(denom) > 0
            shift = 0.5 * (y1 - y3) / denom;
            dk    = kw(j) - kw(max(j-1,1));
            k_RR  = k_RR + shift * dk;
        end
    end
end
