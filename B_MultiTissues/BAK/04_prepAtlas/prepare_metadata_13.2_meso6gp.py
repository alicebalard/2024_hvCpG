#!/usr/bin/env python3
"""
Prepare metadata for WGBS Atlas analysis
Adds 'Analysis group' column and writes a new CSV
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

# Select only rows where Germ layer = Meso
df = df[df["Germ layer"] == "Meso"]

# Among those groups, keep only groups with at least 3 individuals (rows)
counts = df.groupby("Analysis group").size()
eligible_groups = counts[counts >= 3].index

RANDOM_SEED = 42  # for reproducibility

# Select 6 samples
n_to_sample = min(6, len(eligible_groups))
sampled_groups = (
    pd.Series(eligible_groups)
    .sample(n=n_to_sample, random_state=RANDOM_SEED, replace=False)
    .tolist()
)

df_sampled = df[df["Analysis group"].isin(sampled_groups)]

# Determine output path
meta_out = args.output if args.output else os.path.splitext(args.meta)[0] + "_13.2_meso6gp.csv"

# Save
df_sampled.to_csv(meta_out, index=False)
print(f"✅ Wrote a modified metadata to: {meta_out}")
