"""sothe_p2 -- SOTHE-P2: the Paper-2 computations of the SOTHE project in Python (NumPy + Matplotlib).

A standalone Python rewrite of the MATLAB Altay pack P2_altay_pack 1.0.0 (68 .m files, docs/PORT_MAP.md), for
UHeM Altay (SLURM) or any computer.  Stages (python p2.py help):
  S0  runtime probe                                   stages.stage_S0
  S1  C-05   gate G4a                                 g4a.py      run_g4a_canonical.m 0.1.1
  S2  C-05b  G4c loss scan, T-10 + T-12 tokens        g4c.py      run_g4c_loss_scan.m 0.1.2 (+ sub-grid peak)
  S3  G-12   provenance of entropy_trajectory.csv     g12.py      stage_S3_g12.m (Ref. [25] solver)
  S4  T-14   oracle ledgers v0.1.0 and v0.2.0         oracles.py  run_oracles.m, run_oracles_v020.m
  S5  T-13, T-18 manuscript figures                   figures.py  make_paper2_figs.m 0.2.1
  S6  phase-2 data products (10 CSV, JSON, 6 figs)    phase2.py   run_phase2_all.m v0.1.0 (+ suite_figs.py)
  S7  Ref. [25] robustness sweep (300 runs)           ref25/src.py  run_robustness_sweep.m
Shared numerics: suite.py (sothe_phase2_matlab v0.1.0), kit.py (counssug g4a_kit), ref25/ (SOTHE_pkg v1.0.0).

Status of every number this package produces: RECORDED (HR-3: only MATLAB R2025b output is canonical).
"""
import os

__version__ = "1.0.0"
PACKAGE = "SOTHE-P2"
PACK_ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
PORT_OF = ("P2_altay_pack 1.0.0 (MATLAB, 68 files): g4a_kit (run_g4a_canonical 0.1.1, run_g4c_loss_scan 0.1.2, "
           "p2_splitstep_loss 0.1.1 + sub-grid peak position (not in the MATLAB source)), sothe_phase2_matlab v0.1.0 "
           "(run_phase2_all, six figure scripts) + v0.2.0 overlay, make_paper2_figs 0.2.1, pack stage drivers; "
           "Ref. [25] = SOTHE_pkg v1.0.0: src/*.m (ported) and validation/gnlse_dimensionless.py (unchanged)")
