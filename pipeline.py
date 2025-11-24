#!/usr/bin/env python3
"""
pipeline.py
Pure Python RNA-RNA sliding-window interaction pipeline.
Includes logging/debugging to track window generation, candidate scoring,
and sequence extraction.
"""

from pathlib import Path
import pandas as pd
from Bio import SeqIO

# ---------- Utilities ----------
def log(msg, callback=None):
    if callback:
        callback(msg)
    else:
        print(msg)

def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)

def read_fasta_to_dict(fasta_path):
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(str(fasta_path), "fasta")}

def write_fasta(records, path):
    with open(path, "w") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n{seq}\n")

# ---------- Window generation ----------
def generate_windows(fasta_path, window, step, log_callback=None):
    seqs = read_fasta_to_dict(fasta_path)
    records = []
    log(f"Generating windows of size {window}, step {step}", log_callback)
    for rid, seq in seqs.items():
        L = len(seq)
        for start in range(0, max(1, L - window + 1), step):
            end = start + window
            wid = f"{rid}|{start+1}-{end}"  # 1-based coords
            records.append({
                'query_id': wid,
                'sequence': seq[start:end],
                'q_start': start + 1,
                'q_end': end
            })
    log(f"Generated {len(records)} windows across {len(seqs)} sequences", log_callback)
    return pd.DataFrame(records)

# ---------- Candidate scoring ----------
def simple_complementarity_score(seqA, seqB):
    """
    Simple score: count of complementary base pairs (A-U, G-C) in alignment.
    Assumes equal-length sequences.
    """
    comp = {'A':'U','U':'A','G':'C','C':'G'}
    score = sum(1 for a,b in zip(seqA, seqB) if comp.get(a.upper(),'')==b.upper())
    return score

def score_candidates(df_queries, seqB_dict, log_callback=None):
    """
    Score each query window against full target sequences in seqB_dict.
    Returns dataframe with candidate hits.
    """
    rows = []
    log(f"Scoring candidates against target sequences ({len(seqB_dict)} sequences)...", log_callback)
    for idx, row in df_queries.iterrows():
        query_seq = row['sequence']
        for tid, target_seq in seqB_dict.items():
            # Sliding along target sequence
            for start in range(0, len(target_seq) - len(query_seq) + 1):
                end = start + len(query_seq)
                target_slice = target_seq[start:end]
                score = simple_complementarity_score(query_seq, target_slice)
                rows.append({
                    'query_id': row['query_id'],
                    'target_id': tid,
                    'q_start': row['q_start'],
                    'q_end': row['q_end'],
                    't_start': start + 1,
                    't_end': end,
                    'score': score
                })
        if idx % 10 == 0:
            log(f"Scored {idx+1}/{len(df_queries)} query windows", log_callback)
    log(f"Completed scoring {len(df_queries)} windows", log_callback)
    return pd.DataFrame(rows)
