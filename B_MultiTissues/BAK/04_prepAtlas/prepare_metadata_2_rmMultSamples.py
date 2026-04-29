#!/usr/bin/env python3
"""
Prepare metadata for WGBS Atlas analysis
Rm samples if multiple per indivdual
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

# Keep only the first occurrence of each PatientID
df_unique = df.drop_duplicates(subset="PatientID", keep="first")

# Determine output path
meta_out = args.output if args.output else os.path.splitext(args.meta)[0] + "_2_rmMultSamples.csv"

# Save
df_unique.to_csv(meta_out, index=False)
print(f"✅ Wrote a modified metadata to: {meta_out}")

