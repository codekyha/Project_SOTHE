function out = compute_dxic(xi, Stot, opts)
%COMPUTE_DXIC  Operational coarse-graining scale and window robustness.
%
%   out = COMPUTE_DXIC(xi, Stot)
%   out = COMPUTE_DXIC(xi, Stot, opts)
%
%   Implements the operational definition of Delta xi_c stated in Sec. 3.2
%   of the R2 revision of CQG-114816, and produces every number required by
%   the manuscript placeholders:
%
%     \DXICVAL  -> out.dxic            (quote the autocorrelation value;
%                                       cross-checked against out.dxic_fft)
%     \DXICROB  -> out.eta_var_percent (max % variation of eta_GSL over
%                                       windows in [dxic, 4*dxic])
%
%   plus the verification flags backing the two robustness sentences:
%
%     out.all_windows_nonneg  -> true iff <dS/dxi>_W >= -tol pointwise for
%                                EVERY window W in [dxic, 4*dxic]
%     out.frac_neg_worst      -> worst-case fraction of window-averaged
%                                samples below -tol (should be 0)
%
%   INPUT
%     xi    : propagation coordinate, vector (may be non-uniform; the
%             script resamples to a uniform grid internally). Restrict to
%             the post-fission window BEFORE calling, or set opts.xi_min.
%     Stot  : total coarse-grained entropy S_hor + S_rad on xi, same size.
%     opts  : (optional struct)
%             .xi_min   post-fission cut applied internally (default 2)
%             .n_uniform number of uniform resampling points (default:
%                        4 x numel(xi), min 2048)
%             .tol      tolerance for "non-negative" window averages,
%                        relative to max|dS| (default 1e-3)
%             .n_windows number of window sizes scanned in [dxic, 4 dxic]
%                        (default 13)
%
%   METHOD
%     1) Resample Stot(xi) to a uniform grid; dS = gradient(Stot, dxi).
%     2) Autocorrelation (manual -- no toolbox needed) of the mean-removed
%        dS; first zero crossing at lag L0  =>  period ~ 4*L0 for a
%        quasi-sinusoidal oscillation  =>  Delta xi_c = 4*L0*dxi.
%     3) Cross-check: FFT of mean-removed dS; dominant low-frequency peak
%        f*  =>  Delta xi_c (FFT) = 1/f*. If the two estimates differ by
%        more than 25%, quote the FFT value and inspect dS visually.
%     4) Robustness scan: for W in linspace(dxic, 4*dxic, n_windows),
%        window-average dS with movmean; record (a) whether the averaged
%        rate is >= -tol everywhere, (b) eta_GSL evaluated on the averaged
%        rate. Report max relative variation of eta over the scan.
%
%   Base MATLAB only (movmean requires R2016a+). Author: revision support
%   tooling for CQG-114816 R2; include in the public repository.

    if nargin < 3, opts = struct(); end
    if ~isfield(opts, 'xi_min'),     opts.xi_min     = 2;    end
    if ~isfield(opts, 'n_uniform'),  opts.n_uniform  = max(4*numel(xi), 2048); end
    if ~isfield(opts, 'tol'),        opts.tol        = 1e-3; end
    if ~isfield(opts, 'n_windows'),  opts.n_windows  = 13;   end

    xi   = xi(:);  Stot = Stot(:);
    assert(numel(xi) == numel(Stot), 'xi and Stot must have equal length.');

    % ---- 0) post-fission cut + uniform resampling --------------------------
    m    = xi >= opts.xi_min;
    xi   = xi(m);  Stot = Stot(m);
    assert(numel(xi) > 32, 'Too few samples after the post-fission cut.');
    xiu  = linspace(xi(1), xi(end), opts.n_uniform).';
    Su   = interp1(xi, Stot, xiu, 'pchip');
    dxi  = xiu(2) - xiu(1);

    % ---- 1) instantaneous rate ---------------------------------------------
    dS   = gradient(Su, dxi);
    f    = dS - mean(dS);                       % oscillatory component

    % ---- 2) autocorrelation period (manual, unbiased-normalized) ----------
    nmax = floor(numel(f)/2);
    r    = zeros(nmax,1);
    r0   = sum(f.^2);
    for L = 1:nmax
        r(L) = sum(f(1:end-L).*f(1+L:end)) / r0;
    end
    iz = find(r <= 0, 1, 'first');              % first zero crossing
    assert(~isempty(iz), 'No zero crossing of the autocorrelation found.');
    dxic_acf = 4 * iz * dxi;

    % ---- 3) FFT cross-check -------------------------------------------------
    F    = abs(fft(f .* hannwin(numel(f)))).^2;
    fr   = (0:numel(f)-1).'/(numel(f)*dxi);     % cycles per unit xi
    half = 2:floor(numel(f)/2);                 % skip DC
    [~, ip] = max(F(half));
    fstar    = fr(half(ip));
    dxic_fft = 1/fstar;

    rel = abs(dxic_acf - dxic_fft) / dxic_fft;
    if rel > 0.25
        warning(['Autocorrelation (%.4g) and FFT (%.4g) estimates differ ', ...
                 'by %.0f%%: quote the FFT value and inspect dS(xi).'], ...
                 dxic_acf, dxic_fft, 100*rel);
        dxic = dxic_fft;
    else
        dxic = dxic_acf;
    end

    % ---- 4) robustness scan over window sizes ------------------------------
    Ws   = linspace(dxic, 4*dxic, opts.n_windows);
    tol  = opts.tol * max(abs(dS));
    eta  = zeros(size(Ws));
    okW  = false(size(Ws));
    fneg = zeros(size(Ws));
    for j = 1:numel(Ws)
        nw  = max(3, round(Ws(j)/dxi));
        mS  = movmean(dS, nw);
        v   = mS(nw:end-nw+1);                  % drop edge-contaminated part
        okW(j)  = all(v >= -tol);
        fneg(j) = mean(v < -tol);
        eta(j)  = trapz(max(v,0)) / trapz(abs(v));
    end
    eta_var = 100 * (max(eta) - min(eta)) / mean(eta);

    % ---- 5) report -----------------------------------------------------------
    out = struct('dxic', dxic, 'dxic_acf', dxic_acf, 'dxic_fft', dxic_fft, ...
                 'eta_windowed', eta, 'windows', Ws, ...
                 'eta_var_percent', eta_var, ...
                 'all_windows_nonneg', all(okW), ...
                 'frac_neg_worst', max(fneg));

    fprintf('\n=========== compute_dxic: results ===========\n');
    fprintf('Delta xi_c (autocorr) : %.4g\n', dxic_acf);
    fprintf('Delta xi_c (FFT)      : %.4g\n', dxic_fft);
    fprintf('Quoted Delta xi_c     : %.4g   -> \\DXICVAL\n', dxic);
    fprintf('eta_GSL over windows  : [%.3f, %.3f], variation %.1f%%  -> \\DXICROB\n', ...
            min(eta), max(eta), eta_var);
    fprintf('<dS>_W >= 0 for all W : %d   (worst neg. fraction %.2e)\n', ...
            out.all_windows_nonneg, out.frac_neg_worst);
    fprintf('Robustness sentence holds as phrased: %d\n', ...
            out.all_windows_nonneg && eta_var < 100);
    fprintf('==============================================\n\n');
end

function w = hannwin(n)
% Hann window without the Signal Processing Toolbox.
    w = 0.5 * (1 - cos(2*pi*(0:n-1).'/(n-1)));
end

