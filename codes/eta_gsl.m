function eta = eta_gsl(xi, S_tot)
% ETA_GSL  Coarse-grained generalised-second-law efficiency.
%
%   eta = eta_gsl(xi, S_tot)
%
%   Returns the cumulative fraction of positive entropy increments along
%   the propagation, weighted by |dS|. The diagnostic was introduced in
%   the revised paper to replace the (withdrawn) pointwise inequality
%   dS_tot/dxi >= 0 with a meaningful coarse-grained measure:
%
%       eta_GSL(xi) =      Sum_{xi' <= xi} max(0,  dS(xi'))
%                     -------------------------------------------
%                          Sum_{xi' <= xi}  |dS(xi')|
%
%   eta = 1/2 is the symmetric-random-walk baseline. The paper reports
%   eta in [0.52, 0.62] across the 300-config sweep; eta > 1/2 by a
%   statistically clear margin is the operational "GSL holds" criterion.
%
%   Inputs
%   ------
%   xi    : monotone-increasing column vector of propagation coordinates.
%   S_tot : same shape as xi; S_hor + S_rad.
%
%   Output
%   ------
%   eta   : same shape as xi; eta(1) = NaN by convention (no step yet),
%           eta(end) is the global efficiency.

    xi    = xi(:);
    S_tot = S_tot(:);
    dS    = diff(S_tot);

    pos     = max(dS, 0);
    abs_dS  = abs(dS);

    cum_pos = cumsum(pos);
    cum_abs = cumsum(abs_dS);
    cum_abs(cum_abs == 0) = eps;       % avoid 0/0 at the very start

    eta            = [NaN; cum_pos ./ cum_abs];
end
