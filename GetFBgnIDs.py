#!/usr/bin/env python3
"""
Canonical entrypoint for FBgnID conversion.

Usage:
    python GetFBgnIDs.py <input_directory> [gene_column]
"""

from HelperScripts.GetFBgnIDs import main


if __name__ == "__main__":
    raise SystemExit(main())
