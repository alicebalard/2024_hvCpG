
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

# Validate columns
for col in ["Source Tissue", "Cell type"]:
    if col not in df.columns:
        raise ValueError(f"Missing column: {col}")

# Add Analysis group
df["Analysis group"] = df["Source Tissue"].astype(str) + " - " + df["Cell type"].astype(str)

# Select only rows where Germ layer = Ecto
df = df[df["Germ layer"] == "Ecto"]

# Determine output path
meta_out = args.output if args.output else os.path.splitext(args.meta)[0] + "_14_ecto.csv"

# Save
df.to_csv(meta_out, index=False)
print(f"✅ Wrote a modified metadata to: {meta_out}")
