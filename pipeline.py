#!/usr/bin/env python3
"""
pipeline.py
Pure Python RNA-RNA sliding-window interaction pipeline.
Tracks window coordinates for Sequence A and B, and dynamically extracts sequences.
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
    """Read FASTA and return dictionary {id: sequence}."""
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(str(fasta_path), "fasta")}

# ---------- Window generation ----------
def generate_windows(fasta_path, window, step, log_callback=None):
    """Generate sliding windows with nucleotide coordinates."""
    seqs = read_fasta_to_dict(fasta_path)
    records = []
    log(f"Generating windows of size {window}, step {step}", log_callback)
    for rid, seq in seqs.items():
        L = len(seq)
        for start in range(0, max(1, L - window + 1), step):
            end = start + window
            wid = f"{rid}|{start+1}-{end}"  # 1-based coordinates
            records.append({
                'query_id': wid,
                'seq_id': rid,
                'q_start': start + 1,
                'q_end': end
            })
    log(f"Generated {len(records)} windows across {len(seqs)} sequences", log_callback)
    return pd.DataFrame(records)

# ---------- Candidate scoring ----------
def simple_complementarity_score(seqA, seqB):
    """
    Simple complementarity score: count of complementary base pairs (A-U, G-C).
    Assumes equal-length sequences.
    """
    comp = {'A':'U','U':'A','G':'C','C':'G'}
    return sum(1 for a,b in zip(seqA, seqB) if comp.get(a.upper(),'')==b.upper())

def score_candidates(df_queries, seqA_dict, seqB_dict, log_callback=None, flank=0):
    """
    Score each query window against full target sequences.
    Keeps track of nucleotide positions for both sequences.
    Returns a dataframe with columns: query_id, target_id, q_start, q_end, t_start, t_end, score.
    """
    rows = []
    log(f"Scoring {len(df_queries)} query windows against {len(seqB_dict)} target sequences...", log_callback)

    for idx, row in df_queries.iterrows():
        query_seq_full = seqA_dict[row['seq_id']]
        q_start = max(0, row['q_start'] - 1 - flank)  # 0-based start including flank
        q_end = min(len(query_seq_full), row['q_end'] + flank)  # end with flank
        query_seq = query_seq_full[q_start:q_end]

        for tid, target_seq_full in seqB_dict.items():
            L_target = len(target_seq_full)
            window_len = len(query_seq)
            for t_start in range(0, L_target - window_len + 1):
                t_end = t_start + window_len
                target_seq = target_seq_full[t_start:t_end]
                score = simple_complementarity_score(query_seq, target_seq)
                rows.append({
                    'query_id': row['query_id'],
                    'seq_id': row['seq_id'],
                    'target_id': tid,
                    'q_start': row['q_start'],
                    'q_end': row['q_end'],
                    't_start': t_start + 1,
                    't_end': t_end,
                    'score': score
                })

        if (idx+1) % 10 == 0 or (idx+1) == len(df_queries):
            log(f"Scored {idx+1}/{len(df_queries)} query windows", log_callback)

    scored_df = pd.DataFrame(rows)
    log(f"Completed scoring. Total candidate hits: {len(scored_df)}", log_callback)
    return scored_df
