#!/usr/bin/env python3
"""SOTHE-P2 entry point:  python p2.py <command> [options]      (python p2.py help)

Thin wrapper around sothe_p2/cli.py, kept at the top of the pack so that the commands work from the pack folder
without installing anything (python -m sothe_p2 <command> does the same)."""
import os
import sys

sys.dont_write_bytecode = True
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))

from sothe_p2.cli import main  # noqa: E402  (after the path set-up)

if __name__ == "__main__":
    main()
