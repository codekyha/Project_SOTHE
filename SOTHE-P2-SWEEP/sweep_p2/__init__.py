"""sweep_p2 -- SOTHE-P2-SWEEP: one full parameter sweep of the Paper-2 computations, in Python (NumPy, SciPy, Matplotlib).

The sweep is a single job.  Its tasks run in one process pool, in three phases:
  A  independent tasks:
       B1 kinematic layer over (N, k_op, delta3, kappa_g)        B3 radiative loss and drift (split step)
       B5 BdG about the stationary TOD soliton phi               B6 Raman / self-steepening extension (RK4IP)
       B7 mode-conversion scattering with the kinematic flow     B8 Ref. [25] sweep variants
       B9 extraction-bias parity tables (Supplement Tables 2, 3)
  B  tasks that need phase A:  B4 comoving channel census and Hawking window (drift epochs from B3)
  C  reductions:  B2 entanglement layer (closed form), B10 anchor products at (N, k_op) = (1, 1.5), tables, figures,
     the oracle ledger (ORACLES.json), SUMMARY.md, MANIFEST.sha256, the results tarball.

Numerics reused unchanged from SOTHE-P2 1.0.0 (vendored in sothe_p2/): kit.splitstep_loss, kit.tod_soliton,
kit.bdg_spectrum_bg, kit.channel_census, suite.kinematic_kappa and the thermal layer, ref25.src (Ref. [25] ports).
New numerics: sweep_p2.b_raman (GNLSE with Raman and self-steepening), sweep_p2.b_scatter (scattering solve).

Status of every number: RECORDED Python output (HR-3).  Nothing here edits a manuscript."""
import os

__version__ = "1.0.0"
PACKAGE = "SOTHE-P2-SWEEP"
PACK_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
RESULTS_PREFIX = "SOTHE-P2-SWEEP_results_"
