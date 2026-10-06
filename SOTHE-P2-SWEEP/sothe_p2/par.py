"""Process-level parallelism for independent tasks.

Every task is a module-level function run in its own process with single-threaded BLAS/FFT, so a
parallel run gives the same numbers as a serial one.  The pool is created once (spawn start method on
every OS, so Linux and Windows behave alike) and reused by all stages of a run."""
import multiprocessing as mp
import sys
import time
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures.process import BrokenProcessPool

_EX = None
_WORKERS = 1
_BROKEN = False              # set when a pool broke: later calls run serially instead of starting another pool


def set_workers(n):
    global _WORKERS
    _WORKERS = max(1, int(n))


def workers():
    return _WORKERS


def _report(label, done, n, t0):
    el = time.time() - t0
    print("  %s: %d/%d done (%.0f s elapsed, ~%.0f s remaining)" % (label, done, n, el, el / done * (n - done)))
    sys.stdout.flush()


def pmap(fn, arglist, progress=None):
    """[fn(*a) for a in arglist], in parallel when more than one worker is configured.
    PROGRESS (a label) prints a progress line about every tenth of the tasks."""
    global _EX, _BROKEN
    arglist = [a if isinstance(a, tuple) else (a,) for a in arglist]
    n = len(arglist)
    step = max(1, n // 10)
    t0 = time.time()
    if _WORKERS <= 1 or n <= 1 or _BROKEN:
        out = []
        for i, a in enumerate(arglist):
            out.append(fn(*a))
            if progress and ((i + 1) % step == 0 or i + 1 == n):
                _report(progress, i + 1, n, t0)
        return out
    if _EX is None:
        _EX = ProcessPoolExecutor(max_workers=_WORKERS, mp_context=mp.get_context("spawn"))
    try:
        futures = [_EX.submit(fn, *a) for a in arglist]
        if progress:
            done = 0
            for fu in as_completed(futures):
                if fu.exception() is not None and isinstance(fu.exception(), BrokenProcessPool):
                    raise fu.exception()
                done += 1
                if done % step == 0 or done == n:
                    _report(progress, done, n, t0)
        return [f.result() for f in futures]
    except BrokenProcessPool:
        # a worker process died (out of memory, or a main module that spawned workers cannot import, such as
        # `python -` from a pipe): record it and finish the tasks serially in this process
        print("  WARNING: the worker pool broke; running the %d tasks (and all later ones) serially" % n)
        sys.stdout.flush()
        shutdown()
        _BROKEN = True
        return [fn(*a) for a in arglist]


def shutdown():
    global _EX
    if _EX is not None:
        _EX.shutdown(wait=True)
        _EX = None
