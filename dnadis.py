#!/usr/bin/env python3
"""
Compatibility shim for dnadis.

This file maintains backwards compatibility for:
    python dnadis.py [args]

All functionality is now in the dnadis/ package.
"""
from __future__ import annotations

# Re-export everything from the package
from dnadis import *  # noqa: F401, F403
from dnadis import main

if __name__ == "__main__":
    main()
