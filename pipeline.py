#!/usr/bin/env python3
"""
pipeline.py
General RNA-RNA sliding-window interaction pipeline using IntaRNA for both
fast scanning (replacement for RiSearch) and final rescoring.
"""

import os
import argparse
import subprocess
from pathlib import Path
import yaml
from Bio import SeqIO
import pandas as pd
from tqdm import tqdm

# ---------- Utilities ----------
def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)

def read_fasta_to_dict(fasta_path):
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(str(fasta_path), "fasta")}

def write_fasta(records, path):
    with open(path, "w") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n{seq}\n")

# ---------- Window generation ----------
def generate_windows(fasta_path, window, step, out_fa):
    seqs = read_fasta_to_dict(fasta_path)
    records = []
    for rid, seq in seqs.items():
        L = len(seq)
        for start in range(0, max(1, L - window + 1), step):
            s = start
            e = start + window
            subseq = seq[s:e]
            wid = f"{rid}|{s+1}-{e}"  # 1-based coords
            records.append((wid, subseq))
    write_fasta(records, out_fa)
    return out_fa

# ---------- Fast search (IntaRNA replacement for RiSearch) ----------
def call_intarna_fast(query_fa, target_fa, out_tsv, energy_cutoff=None, max_hits=None,
                      seedTQ=3, seedBP=8, threads=1, intarna_bin=None):
    """Run IntaRNA in 'fast scan' mode to produce candidate interactions."""
    if intarna_bin is None:
        intarna_bin = os.path.join(os.path.dirname(__file__), "bin", "IntaRNA")
    out_tsv = str(out_tsv)
    cmd = [
        intarna_bin,
        "--query", str(query_fa),
        "--target", str(target_fa),
        "--outMode", "C",
        "--out", out_tsv,
        "--threads", str(threads),
        "--seedTQ", str(seedTQ),
        "--seedBP", str(seedBP)
    ]
    if energy_cutoff is not None:
        cmd += ["--energyThreshold", str(energy_cutoff)]
    print("Running IntaRNA fast search:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return out_tsv

# ---------- Compatibility wrapper for pipeline_runner ----------
def call_vienna_fast(query_fa, target_fa, out_tsv, energy_cutoff=-6.0, max_hits=None, threads=1, vienna_bin=None):
    """
    Wrapper to call IntaRNA fast scan, compatible with previous RiSearch pipeline.
    """
    return call_intarna_fast(
        query_fa=query_fa,
        target_fa=target_fa,
        out_tsv=out_tsv,
        energy_cutoff=energy_cutoff,
        max_hits=max_hits,
        threads=threads,
        intarna_bin=vienna_bin
    )

# ---------- Extract flanking context ----------
def extract_contexts(target_fa, hits_df, flank, out_dir):
    ensure_dir(out_dir)
    seqs = read_fasta_to_dict(target_fa)
    rows = []
    for idx, row in hits_df.iterrows():
        tid = str(row['target_id'])
        start = int(row['t_start']) if pd.notna(row['t_start']) else None
        end = int(row['t_end']) if pd.notna(row['t_end']) else None
        if start is None or end is None:
            continue
        seq = seqs.get(tid)
        if seq is None:
            raise RuntimeError(f"target id {tid} not found in {target_fa}")
        L = len(seq)
        s = max(0, start - 1 - flank)
        e = min(L, end + flank)
        context_seq = seq[s:e]
        fid = f"{tid}|{start}-{end}|flank{flank}"
        outp = Path(out_dir) / f"{fid}.fa"
        write_fasta([(fid, context_seq)], outp)
        rows.append({
            'hit_index': idx,
            'target_id': tid,
            't_start': start,
            't_end': end,
            'context_fasta': str(outp),
        })
    return pd.DataFrame(rows)

# ---------- Rescore with IntaRNA ----------
def call_intarna(query_fa, target_context_fa, out_prefix, threads=1, intarna_bin=None):
    """Run IntaRNA for rescoring (full prediction)."""
    if intarna_bin is None:
        intarna_bin = os.path.join(os.path.dirname(__file__), "bin", "IntaRNA")
    out_prefix = str(out_prefix)
    out_csv = out_prefix + ".csv"
    cmd = [
        intarna_bin,
        "--query", str(query_fa),
        "--target", str(target_context_fa),
        "--outMode", "C",
        "--out", out_csv,
        "--threads", str(threads)
    ]
    subprocess.run(cmd, check=True)
    return out_csv

def parse_intarna_csv(csv_path):
    df = pd.read_csv(csv_path, comment="#")
    # canonicalize energy column
    if 'E' in df.columns and 'energy' not in df.columns:
        df = df.rename(columns={'E':'energy'})
    if 'Energy' in df.columns and 'energy' not in df.columns:
        df = df.rename(columns={'Energy':'energy'})
    return df

# ---------- Aggregation & ranking ----------
def aggregate_and_rank(intarna_results, hits_df, out_tsv):
    rows = []
    for hit_index, csvp in intarna_results:
        df = parse_intarna_csv(csvp)
        if df.shape[0] == 0:
            continue
        energy_col = 'energy' if 'energy' in df.columns else df.columns[0]
        best_row = df.loc[df[energy_col].idxmin()].to_dict()
        merged = hits_df.loc[hit_index].to_dict()
        merged.update({f"intarna_{k}": v for k, v in best_row.items()})
        merged['intarna_csv'] = csvp
        rows.append(merged)
    outdf = pd.DataFrame(rows)
    outdf.to_csv(out_tsv, index=False)
    return outdf
