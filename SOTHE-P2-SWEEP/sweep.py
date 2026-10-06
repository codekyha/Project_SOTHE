#!/usr/bin/env python3
"""sweep.py -- single entry point of SOTHE-P2-SWEEP: one full parameter sweep of the Paper-2 computations in Python,
run as ONE SLURM job on UHeM Altay (or on any computer).  Commands: python sweep.py help"""
import os
import sys

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sweep_p2.cli import main  # noqa: E402  (after the path insert)

if __name__ == "__main__":
    main()
