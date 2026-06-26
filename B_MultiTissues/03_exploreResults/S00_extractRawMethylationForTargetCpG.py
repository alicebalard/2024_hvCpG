#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
S00_extractRawMethylationForTargetCpG.py - Extract methylation values per patient
==================================================================================
Given a list of CpG sites (one chr_pos per line), extracts raw beta values
from all matching WGBS .beta files and returns a long-format dataframe:

    cpg_site | patient_id | source_tissue_celltype | germ_layer | methylation

Usage
-----
python S00_extractRawMethylationForTargetCpG.py \
    --cpg_list   /path/to/my_cpg_sites.txt \
    --cpg_bed    /path/to/hg38/CpG.bed.gz \
    --beta_files "/path/to/betaFiles/GSM*.hg38.beta" \
    --meta       /path/to/SupTab1_Loyfer2023_amended.csv \
    --output     /path/to/output.tsv \
    --minCov     10

Input CpG list format (one per line, chr_pos):
    chr1_17752162
    chr1_17752319
    chr1_17752421

Author: Alice Balard
"""

import os
import sys
import glob
import argparse

import numpy as np
import pandas as pd

sys.path.insert(0, "/SAN/ghlab/epigen/Alice/hvCpG_project/code/2024_hvCpG/B_MultiTissues/01_dataPrep")
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
parser.add_argument("--cpg_list",       required=True,
                    help="Text file with one CpG site per line in chr_pos format.")
parser.add_argument("--cpg_bed",        required=True,
                    help="CpG BED reference file (.bed or .bed.gz).")
parser.add_argument("--beta_files",     required=True,
                    help="Glob pattern for .beta files.")
parser.add_argument("--meta",           required=True,
                    help="Metadata CSV.")
parser.add_argument("--output",         required=True,
                    help="Output TSV path.")
parser.add_argument("--minCov",         type=int, default=10,
                    help="Minimum read coverage (default: 10). Below = NA.")
parser.add_argument("--id_pattern",     default=r"-([A-Za-z0-9_]+)\.hg38\.beta$",
                    help="Regex to extract sample ID from beta filename.")
parser.add_argument("--sample_col",     default="Sample name",
                    help="Metadata column for sample names (default: 'Sample name').")
parser.add_argument("--patient_col",    default="PatientID",
                    help="Metadata column for patient ID (default: 'PatientID').")
parser.add_argument("--tissue_col",     default="Source Tissue",
                    help="Metadata column for source tissue (default: 'Source Tissue').")
parser.add_argument("--cell_col",       default="Cell type",
                    help="Metadata column for cell type (default: 'Cell type').")
parser.add_argument("--germ_layer_col", default="Germ layer",
                    help="Metadata column for germ layer (default: 'Germ layer').")
parser.add_argument("--id_sep",         default="-",
                    help="Separator to shorten sample IDs (default: '-').")

args = parser.parse_args()

# ──────────────────────────────────────────────
#  Load inputs
# ──────────────────────────────────────────────

print("\n" + "="*62)
print("  S00_extractRawMethylationForTargetCpG.py")
print("="*62)

# 1. Target CpGs
with open(args.cpg_list) as f:
    target_cpgs = [l.strip() for l in f if l.strip() and not l.startswith("#")]
print(f"  Target CpGs : {len(target_cpgs):,}")

# 2. BED reference -> index
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
    print(f"  Warning: {len(missing):,} CpG(s) not found in reference -- skipped.")
    for m in missing[:10]:
        print(f"    {m}")
    if len(missing) > 10:
        print(f"    ... and {len(missing) - 10} more.")

sorted_cpgs    = sorted(target_indices.keys(), key=lambda c: target_indices[c])
sorted_indices = np.array([target_indices[c] for c in sorted_cpgs])
print(f"  Matched {len(sorted_cpgs):,} / {len(target_cpgs):,} CpGs.")

# 3. Beta files
beta_files = sorted(glob.glob(args.beta_files))
if not beta_files:
    raise FileNotFoundError(f"No beta files matched: {args.beta_files}")
print(f"  Beta files  : {len(beta_files):,}")
sample_to_path = build_sample_to_path_map(beta_files, id_pattern=args.id_pattern)

# 4. Metadata
meta = pd.read_csv(args.meta)
for col in [args.sample_col, args.patient_col, args.tissue_col,
            args.cell_col, args.germ_layer_col]:
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
id_to_germ    = dict(zip(meta["_short_id"], meta[args.germ_layer_col]))

# ──────────────────────────────────────────────
#  Extract: one row per (patient, tissue_cell, cpg)
# ──────────────────────────────────────────────

# print(f"\n  Extracting from {len(sample_to_path):,} beta files...")
# 
# rows = []
# for i, (short_id, path) in enumerate(sample_to_path.items()):
#     if (i + 1) % 50 == 0:
#         print(f"    {i+1:,} / {len(sample_to_path):,} processed...")
# 
#     patient_id  = id_to_patient.get(short_id)
#     tissue_cell = id_to_label.get(short_id)
# 
#     if patient_id is None or tissue_cell is None:
#         print(f"    Skipping {short_id}: not found in metadata.")
#         continue
# 
#     mm   = np.memmap(path, dtype=np.uint8, mode="r").reshape(-1, 2)
#     meth = mm[sorted_indices, 0].astype(np.float32)
#     cov  = mm[sorted_indices, 1]
#     beta = np.where(cov == 0, np.nan, meth / cov).astype(np.float32)
#     beta[cov < args.minCov] = np.nan
# 
#     for cpg, b in zip(sorted_cpgs, beta):
#         rows.append({
#             "cpg_site"              : cpg,
#             "patient_id"            : patient_id,
#             "source_tissue_celltype": tissue_cell,
#             "methylation"           : None if np.isnan(b) else round(float(b), 6),
#         })
# ── Replace the extraction loop with a batched, streaming version ─────────────

import csv

print(f"\n  Extracting from {len(sample_to_path):,} beta files...")

os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)

with open(args.output, "w", newline="") as fout:
    writer = csv.writer(fout, delimiter="\t")
    writer.writerow(["cpg_site", "patient_id", "source_tissue_celltype", "methylation"])

    for i, (short_id, path) in enumerate(sample_to_path.items()):
        if (i + 1) % 50 == 0:
            print(f"    {i+1:,} / {len(sample_to_path):,} processed...")

        patient_id  = id_to_patient.get(short_id)
        tissue_cell = id_to_label.get(short_id)

        if patient_id is None or tissue_cell is None:
            print(f"    Skipping {short_id}: not found in metadata.")
            continue

        try:
            # open memmap, extract only target rows, then close immediately
            mm   = np.memmap(path, dtype=np.uint8, mode="r")
            mm   = mm.reshape(-1, 2)
            meth = mm[sorted_indices, 0].astype(np.float32)
            cov  = mm[sorted_indices, 1]
            del mm   # release memmap immediately after extraction

            with np.errstate(invalid="ignore", divide="ignore"):
                beta = np.where(cov == 0, np.nan, meth / cov).astype(np.float32)
            beta[cov < args.minCov] = np.nan

            for cpg, b in zip(sorted_cpgs, beta):
                mval = None if np.isnan(b) else round(float(b), 6)
                writer.writerow([cpg, patient_id, tissue_cell, mval if mval is not None else "NA"])

        except OSError as e:
            print(f"    ERROR on {short_id} ({path}): {e} — skipping.")
            continue

print(f"\n  Output: {args.output}")
print("\n🎉 Done!")
