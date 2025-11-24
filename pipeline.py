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
        for start in range(0, len(seq) - window + 1, step):
            end = start + window
            wid = f"{rid}|{start+1}-{end}"  # 1-based coordinates
            records.append({
                'query_id': wid,
                'seq_id': rid,
                'sequence': seq[start:end],
                'q_start': start + 1,
                'q_end': end
            })
    log(f"Generated {len(records)} windows across {len(seqs)} sequences", log_callback)
    return pd.DataFrame(records)

# ---------- Scoring ----------
def simple_complementarity_score(seqA, seqB):
    """Count complementary base pairs (A-U, G-C)"""
    comp = {'A':'U','U':'A','G':'C','C':'G'}
    return sum(1 for a,b in zip(seqA, seqB) if comp.get(a.upper(),'')==b.upper())

def quick_filter(df_queries, seqB_dict, top_k=10, log_callback=None):
    """Quick first-pass: score all windows and keep top_k per query window"""
    log(f"Running quick filter (top {top_k} hits per query window)...", log_callback)
    candidate_rows = []
    for idx, row in df_queries.iterrows():
        query_seq = row['sequence']
        for tid, target_seq in seqB_dict.items():
            for t_start in range(0, len(target_seq) - len(query_seq) + 1):
                t_end = t_start + len(query_seq)
                slice_target = target_seq[t_start:t_end]
                score = simple_complementarity_score(query_seq, slice_target)
                candidate_rows.append({
                    'query_id': row['query_id'],
                    'seq_id': row['seq_id'],
                    'target_id': tid,
                    'q_start': row['q_start'],
                    'q_end': row['q_end'],
                    't_start': t_start+1,
                    't_end': t_end,
                    'score': score
                })
        if (idx+1) % 10 == 0 or (idx+1) == len(df_queries):
            log(f"Quick-scored {idx+1}/{len(df_queries)} windows", log_callback)
    
    all_hits = pd.DataFrame(candidate_rows)
    # Keep top_k per query
    top_hits = all_hits.groupby('query_id').apply(lambda x: x.nlargest(top_k, 'score')).reset_index(drop=True)
    log(f"Quick filtering done. Kept {len(top_hits)} top candidates", log_callback)
    return top_hits

def full_rescore(top_hits_df, seqA_dict, seqB_dict, log_callback=None):
    """Rescore top candidates using full sequences (same as simple complementarity here)"""
    log("Starting full rescoring on top candidates...", log_callback)
    rescored_rows = []
    for idx, row in top_hits_df.iterrows():
        query_seq_full = seqA_dict[row['seq_id']]
        target_seq_full = seqB_dict[row['target_id']]

        # extract sequences based on nucleotide positions
        q_seq = query_seq_full[row['q_start']-1:row['q_end']]
        t_seq = target_seq_full[row['t_start']-1:row['t_end']]
        score = simple_complementarity_score(q_seq, t_seq)
        rescored_rows.append({
            **row,
            'query_seq': q_seq,
            'target_seq': t_seq,
            'rescore_energy_B': score
        })
        if (idx+1) % 10 == 0 or (idx+1) == len(top_hits_df):
            log(f"Rescored {idx+1}/{len(top_hits_df)} candidates", log_callback)
    df_rescored = pd.DataFrame(rescored_rows)
    log("Full rescoring done.", log_callback)
    return df_rescored
