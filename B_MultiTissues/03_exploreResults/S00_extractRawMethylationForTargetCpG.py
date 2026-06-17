#!/usr/bin/env python3
"""
extract_cpg_methylation_byPatient.py — Extract methylation values per patient
==============================================================================
Same as extract_cpg_methylation.py but uses PatientID instead of Sample name.
One row per (cpg_site, patient_id, source_tissue_celltype).

Output columns:
    cpg_site | patient_id | source_tissue_celltype | methylation

Usage
-----
python extract_cpg_methylation_byPatient.py \\
    --cpg_list   /path/to/my_cpg_sites.txt \\
    --cpg_bed    /data/hg38/CpG.bed.gz \\
    --beta_files "/data/betaFiles/GSM*.hg38.beta" \\
    --meta       /path/to/SupTab1_Loyfer2023_amended.csv \\
    --output     /path/to/output.tsv \\
    --minCov     10

Author: Alice Balard
"""

import os
import sys
import glob
import argparse

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from cpg_utils import (
    load_cpg_names_from_bed,
    build_sample_to_path_map,
)

# ──────────────────────────────────────────────
#  Arguments
# ──────────────────────────────────────────────

parser = argparse.ArgumentParser(
    description="Extract methylation values for target CpGs, one row per patient.",
    formatter_class=argparse.RawDescriptionHelpFormatter,
    epilog=__doc__,
)
parser.add_argument("--cpg_list",    required=True)
parser.add_argument("--cpg_bed",     required=True)
parser.add_argument("--beta_files",  required=True)
parser.add_argument("--meta",        required=True)
parser.add_argument("--output",      required=True)
parser.add_argument("--minCov",      type=int, default=10)
parser.add_argument("--id_pattern",  default=r"-([A-Za-z0-9_]+)\.hg38\.beta$")
parser.add_argument("--sample_col",  default="Sample name")
parser.add_argument("--patient_col", default="PatientID")
parser.add_argument("--tissue_col",  default="Source Tissue")
parser.add_argument("--cell_col",    default="Cell type")
parser.add_argument("--id_sep",      default="-")

args = parser.parse_args()

# ──────────────────────────────────────────────
#  Load inputs
# ──────────────────────────────────────────────

print("\n" + "="*62)
print("  extract_cpg_methylation_byPatient.py")
print("="*62)

# 1. Target CpGs
with open(args.cpg_list) as f:
    target_cpgs = [l.strip() for l in f if l.strip() and not l.startswith("#")]
print(f"  Target CpGs : {len(target_cpgs):,}")

# 2. BED reference → index
cpg_names = load_cpg_names_from_bed(args.cpg_bed)
cpg_index = {name: i for i, name in enumerate(cpg_names)}

target_indices = {}
missing = []
for cpg in target_cpgs:
    idx = cpg_index.get(cpg)
    if idx is None:
        missing.append(cpg)
    else:
        target_indices[cpg] = idx

if missing:
    print(f"  Warning: {len(missing):,} CpG(s) not found in reference — skipped.")

sorted_cpgs    = sorted(target_indices.keys(), key=lambda c: target_indices[c])
sorted_indices = np.array([target_indices[c] for c in sorted_cpgs])
print(f"  Matched {len(sorted_cpgs):,} / {len(target_cpgs):,} CpGs.")

# 3. Beta files
beta_files = sorted(glob.glob(args.beta_files))
if not beta_files:
    raise FileNotFoundError(f"No beta files matched: {args.beta_files}")
sample_to_path = build_sample_to_path_map(beta_files, id_pattern=args.id_pattern)

# 4. Metadata — short_id → (patient_id, tissue_cell)
meta = pd.read_csv(args.meta)
for col in [args.sample_col, args.patient_col, args.tissue_col, args.cell_col]:
    if col not in meta.columns:
        raise ValueError(f"Metadata missing required column: '{col}'")

def shorten_id(name, sep):
    parts = str(name).split(sep)
    return parts[-1] if len(parts) > 1 else name

meta["_short_id"]    = meta[args.sample_col].apply(lambda x: shorten_id(x, args.id_sep))
meta["_tissue_cell"] = (meta[args.tissue_col].astype(str)
                        + " - " + meta[args.cell_col].astype(str))

id_to_patient = dict(zip(meta["_short_id"], meta[args.patient_col]))
id_to_label   = dict(zip(meta["_short_id"], meta["_tissue_cell"]))

# ──────────────────────────────────────────────
#  Extract: one row per (patient, tissue_cell, cpg)
# ──────────────────────────────────────────────

print(f"\n  Extracting from {len(sample_to_path):,} beta files...")

rows = []
for i, (short_id, path) in enumerate(sample_to_path.items()):
    if (i + 1) % 50 == 0:
        print(f"    {i+1:,} / {len(sample_to_path):,} processed...")

    patient_id  = id_to_patient.get(short_id)
    tissue_cell = id_to_label.get(short_id)

    if patient_id is None or tissue_cell is None:
        print(f"    Skipping {short_id}: not found in metadata.")
        continue

    mm   = np.memmap(path, dtype=np.uint8, mode="r").reshape(-1, 2)
    meth = mm[sorted_indices, 0].astype(np.float32)
    cov  = mm[sorted_indices, 1]
    beta = np.where(cov == 0, np.nan, meth / cov).astype(np.float32)
    beta[cov < args.minCov] = np.nan

    for cpg, b in zip(sorted_cpgs, beta):
        rows.append({
            "cpg_site"              : cpg,
            "patient_id"            : patient_id,
            "source_tissue_celltype": tissue_cell,
            "methylation"           : None if np.isnan(b) else round(float(b), 6),
        })

# ──────────────────────────────────────────────
#  Output
# ──────────────────────────────────────────────

df = pd.DataFrame(rows, columns=[
    "cpg_site", "patient_id", "source_tissue_celltype", "methylation"
])
df = df.sort_values(["cpg_site", "patient_id", "source_tissue_celltype"])

os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
df.to_csv(args.output, sep="\t", index=False, na_rep="NA")

print(f"\n  Output          : {args.output}")
print(f"  Rows            : {len(df):,}")
print(f"  Unique patients : {df['patient_id'].nunique():,}")
print(f"  Unique tissues  : {df['source_tissue_celltype'].nunique():,}")
print(f"  NA              : {df['methylation'].isna().sum():,} "
      f"({df['methylation'].isna().mean()*100:.1f}%)")
print("\n🎉 Done!")
