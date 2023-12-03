#!/usr/bin/env python3

"""
Convenience wrapper for running censo directly from the source tree.
"""
import sys
from src.censo import main

if __name__ == "__main__":
    sys.exit(main())
