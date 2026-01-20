#!/usr/bin/env python3
"""
Compatibility shim for final_finalizer.

This file maintains backwards compatibility for:
    python final_finalizer.py [args]

All functionality is now in the final_finalizer/ package.
"""
from __future__ import annotations

# Re-export everything from the package
from final_finalizer import *  # noqa: F401, F403
from final_finalizer import main

if __name__ == "__main__":
    main()
