#!/usr/bin/env python3
"""
Prepare metadata for WGBS Atlas analysis
Select only immune cells
Author: Alice Balard
"""

import argparse
import pandas as pd
import os

parser = argparse.ArgumentParser(description="Prepare metadata for WGBS Atlas")
parser.add_argument("--meta", required=True, help="Input metadata CSV file")
parser.add_argument("--output", required=False, help="Output file path")
args = parser.parse_args()

# Read input metadata
df = pd.read_csv(args.meta)

# Add Analysis group
df["Analysis group"] = df["Source Tissue"].astype(str) + " - " + df["Cell type"].astype(str)

# Keep only rows where "Immune?" is True
df_filtered = df[df["Immune?"] == True]

# Determine output path
meta_out = args.output if args.output else os.path.splitext(args.meta)[0] + "_9_immuneOnly.csv"

# Save
df_filtered.to_csv(meta_out, index=False)
print(f"✅ Wrote a modified metadata to: {meta_out}")
