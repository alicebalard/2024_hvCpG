#!/usr/bin/env python3
"""
Prepare metadata for WGBS Atlas analysis (Pairs: FF)
Select exactly 2 females per Analysis group (if possible).
Author: Alice Balard (amended)
"""

import argparse
import os
import pandas as pd
from numpy.random import default_rng

def normalize_sex(s):
    if pd.isna(s):
        return pd.NA
    s = str(s).strip().upper()
    if s in {"M", "MALE"}:
        return "M"
    if s in {"F", "FEMALE"}:
        return "F"
    return pd.NA

def main():
    parser = argparse.ArgumentParser(description="Prepare metadata pairs: FF")
    parser.add_argument("--meta", required=True, help="Input metadata CSV file")
    parser.add_argument("--output", required=False, help="Output CSV")
    parser.add_argument("--n_groups", type=int, default=None,
                        help="Max number of Analysis groups to include (default: all eligible)")
    parser.add_argument("--seed", type=int, default=42, help="Random seed")
    args = parser.parse_args()

    rng = default_rng(args.seed)

    df = pd.read_csv(args.meta)

    # Build Analysis group
    df["Analysis group"] = df["Source Tissue"].astype(str) + " - " + df["Cell type"].astype(str)

    # Normalize sex
    df["sex"] = df["sex"].map(normalize_sex)

    # Eligible groups: groups with at least 2 females
    selected = []
    eligible_groups = []
    for grp, g in df.groupby("Analysis group", dropna=False):
        females = g[g["sex"] == "F"]
        if len(females) >= 2:
            eligible_groups.append(grp)

    # Optionally downsample number of groups
    if args.n_groups is not None and args.n_groups < len(eligible_groups):
        eligible_groups = list(rng.choice(eligible_groups, size=args.n_groups, replace=False))

    # From each eligible group, pick exactly 2 random females
    for grp in eligible_groups:
        g = df[df["Analysis group"] == grp]
        females = g[g["sex"] == "F"]
        pick_idx = rng.choice(females.index, size=2, replace=False)
        selected.append(df.loc[pick_idx])

    if selected:
        out_df = pd.concat(selected, axis=0).sort_values(["Analysis group"])
    else:
        out_df = pd.DataFrame(columns=df.columns)

    # Default output name
    meta_out = args.output if args.output else os.path.splitext(args.meta)[0] + "_5_pairs_FF.csv"
    out_df.to_csv(meta_out, index=False)

    print(f"✅ Pairs FF: {out_df['Analysis group'].nunique()} groups, {len(out_df)} samples")
    print(f"✅ Wrote modified metadata to: {meta_out}")

if __name__ == "__main__":
    main()
