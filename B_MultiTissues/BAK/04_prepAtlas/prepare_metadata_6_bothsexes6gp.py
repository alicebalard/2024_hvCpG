#!/usr/bin/env python3
"""
Prepare metadata for WGBS Atlas analysis
Select all but only males groups
Author: Alice Balard
"""

import argparse
import pandas as pd
import os

RANDOM_SEED = 42  # for reproducibility

parser = argparse.ArgumentParser(description="Prepare metadata for WGBS Atlas")
parser.add_argument("--meta", required=True, help="Input metadata CSV file")
parser.add_argument("--output", required=False, help="Output file path")
args = parser.parse_args()

# Read input metadata
df = pd.read_csv(args.meta)

# Add Analysis group
df["Analysis group"] = df["Source Tissue"].astype(str) + " - " + df["Cell type"].astype(str)

# Ensure consistent formatting (optional)
df["sex"] = df["sex"].str.strip()

# Keep groups where sexes are mixed (at least one M and one F)
df_mixed = df.groupby("Analysis group").filter(
    lambda g: ("M" in g["sex"].values) and ("F" in g["sex"].values)
)

# Among those, keep only groups with at least 3 individuals
counts = df_mixed.groupby("Analysis group").size()
eligible_groups = counts[counts >= 3].index

# Sample 6 groups (or fewer if not enough eligible)
n_to_sample = min(6, len(eligible_groups))
sampled_groups = (
    pd.Series(eligible_groups)
    .sample(n=n_to_sample, random_state=RANDOM_SEED, replace=False)
    .tolist()
)

# Keep only rows from sampled groups
df_filtered = df_mixed[df_mixed["Analysis group"].isin(sampled_groups)]

# Determine output path
meta_out = args.output if args.output else os.path.splitext(args.meta)[0] + "_6_bothsexes6gp.csv"

# Save
df_filtered.to_csv(meta_out, index=False)
print(f"✅ Wrote a modified metadata to: {meta_out}")
