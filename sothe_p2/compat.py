"""MATLAB semantics that NumPy does not share, used wherever the MATLAB source relies on them.

  mlinspace     MATLAB linspace(a, b, n): a + (k*(b - a))/(n - 1), endpoints exact.  numpy.linspace computes
                a + k*((b - a)/(n - 1)) and differs in the last bit at about a third of the points; the
                MATLAB R2025b CSV products of the suite match mlinspace at every point (docs/PORT_MAP.md)
  mround        MATLAB round: halves go away from zero (Python and NumPy round halves to even)
  mgradient     MATLAB gradient(F, x): one-sided first differences at the ends, (F(i+1)-F(i-1))/(x(i+1)-x(i-1))
  mtrapz        MATLAB trapz(x, y) for vectors: diff(x).' * (y(1:end-1) + y(2:end)) / 2
  mcumtrapz     MATLAB cumtrapz(x, y) for vectors: [0; cumsum(diff(x)/2 .* (y(1:end-1) + y(2:end)))]
  interp1       MATLAB interp1(x, y, xq) 'linear': NaN outside [x(1), x(end)]
  interp1_extrap  interp1(..., 'linear', 'extrap')
  conv_same     conv(u, v, 'same')
  sech          1 ./ cosh(x), as MATLAB computes it
  trapezoid     numpy.trapezoid (NumPy >= 2.0) or numpy.trapz (NumPy 1.x, e.g. 1.26.4 on Altay)
  mfmt          sprintf with MATLAB's spelling of non-finite values (NaN, Inf, -Inf)
Bit-for-bit agreement with MATLAB is not a goal where it depends on the math library (tanh, exp, pow) or on
LAPACK (eig, least squares): those differ at round-off level between any two runtimes."""
import math
import re

import numpy as np

if not hasattr(np, "trapezoid"):              # NumPy < 2.0: the Ref. [25] validation port calls np.trapezoid
    np.trapezoid = np.trapz                    # noqa: NPY201  (only reached on NumPy 1.x)


def trapezoid(y, x=None, axis=-1):
    return np.trapezoid(y, x, axis=axis)


def mlinspace(a, b, n):
    """MATLAB linspace(a, b, n) (R2025b formula): row of n points, y(1) = a, y(n) = b exactly."""
    n = int(math.floor(float(n)))
    a = float(a)
    b = float(b)
    if n < 1:
        return np.zeros(0)
    if n == 1:
        return np.array([b])                   # MATLAB: linspace(a, b, 1) = b
    k = np.arange(n, dtype=float)
    y = a + (k * (b - a)) / (n - 1)
    if a == b:
        y[:] = a
    else:
        y[0] = a
        y[-1] = b
    return y


def mround(x):
    """MATLAB round for a scalar: exact, halves away from zero."""
    x = float(x)
    if not math.isfinite(x):
        return x
    a = abs(x)
    f = math.floor(a)
    r = f + 1.0 if a - f >= 0.5 else f
    return math.copysign(r, x)


def mround_array(x):
    """MATLAB round, elementwise (halves away from zero)."""
    x = np.asarray(x, dtype=float)
    return np.sign(x) * np.floor(np.abs(x) + 0.5)


def sech(x):
    return 1.0 / np.cosh(x)


def mgradient(f, x):
    f = np.asarray(f, dtype=float).ravel()
    x = np.asarray(x, dtype=float).ravel()
    n = f.size
    g = np.empty(n)
    g[0] = (f[1] - f[0]) / (x[1] - x[0])
    g[-1] = (f[-1] - f[-2]) / (x[-1] - x[-2])
    if n > 2:
        g[1:-1] = (f[2:] - f[:-2]) / (x[2:] - x[:-2])
    return g


def mtrapz(x, y):
    """MATLAB trapz(x, y) for vectors (y may also be a matrix: columns are integrated, as in MATLAB)."""
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y)
    if y.ndim == 1:
        return np.dot(np.diff(x), y[:-1] + y[1:]) / 2
    return np.dot(np.diff(x), y[:-1, :] + y[1:, :]) / 2


def mcumtrapz(x, y):
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    dt = np.diff(x) / 2
    return np.concatenate([[0.0], np.cumsum(dt * (y[:-1] + y[1:]))])


def interp1(x, y, xq):
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    q = np.asarray(xq, dtype=float)
    out = np.interp(q, x, y)
    out = np.where((q < x[0]) | (q > x[-1]) | np.isnan(q), np.nan, out)
    return float(out) if out.ndim == 0 else out


def interp1_extrap(x, y, xq):
    x = np.asarray(x, dtype=float).ravel()
    y = np.asarray(y, dtype=float).ravel()
    xq = float(xq)
    if x[0] <= xq <= x[-1]:
        return float(np.interp(xq, x, y))
    j = 0 if xq < x[0] else x.size - 2
    return float(y[j] + (y[j + 1] - y[j]) * (xq - x[j]) / (x[j + 1] - x[j]))


def conv_same(u, v):
    u = np.asarray(u, dtype=float).ravel()
    v = np.asarray(v, dtype=float).ravel()
    full = np.convolve(u, v, mode="full")
    start = (v.size - 1) // 2 if v.size % 2 == 1 else v.size // 2
    return full[start:start + u.size]


def first_true(mask):
    """Index of the first True (MATLAB find(mask, 1, 'first')), or None."""
    idx = np.flatnonzero(np.asarray(mask))
    return int(idx[0]) if idx.size else None


_NONFINITE = re.compile(r"(?<![A-Za-z])(-?)(nan|inf)(?![A-Za-z])")


def mfmt(fmt, *args):
    """fmt % args, with NaN/Inf spelled as MATLAB's sprintf spells them (Python writes nan/inf)."""
    s = fmt % args
    return _NONFINITE.sub(lambda m: m.group(1) + ("NaN" if m.group(2) == "nan" else "Inf"), s)
