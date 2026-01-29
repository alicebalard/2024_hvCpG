#!/usr/bin/env python3
"""
Prepare metadata for WGBS Atlas analysis
Select all but immune cells
Author: Alice Balard
"""

import argparse
import pandas as pd
import os
import numpy as np

parser = argparse.ArgumentParser(description="Prepare metadata for WGBS Atlas")
parser.add_argument("--meta", required=True, help="Input metadata CSV file")
parser.add_argument("--output", required=False, help="Output file path")
args = parser.parse_args()

# Read input metadata
df = pd.read_csv(args.meta)

# Add Analysis group
df["Analysis group"] = df["Source Tissue"].astype(str) + " - " + df["Cell type"].astype(str)

# Keep only rows where "Immune?" is False
df_filtered = df[df["Immune?"] == False]

# set a seed for reproducibility
RANDOM_SEED = 42

# Restrict to groups with at least 3 samples
group_counts = df_filtered.groupby("Analysis group").size()
eligible_groups = group_counts[group_counts >= 3].index

# Sample up to 11 groups (or fewer if not enough eligible)
n_to_sample = min(11, len(eligible_groups))
sampled_groups = (
    pd.Series(eligible_groups)
    .sample(n=n_to_sample, random_state=RANDOM_SEED, replace=False)
    .tolist()
)

# Keep only rows belonging to the sampled groups
df_sampled = df_filtered[df_filtered["Analysis group"].isin(sampled_groups)]

# Determine output path
meta_out = args.output if args.output else os.path.splitext(args.meta)[0] + "_11_noImmune_sample11groups.csv"

# Save
df_sampled.to_csv(meta_out, index=False)
print(f"✅ Wrote a modified metadata to: {meta_out}")
