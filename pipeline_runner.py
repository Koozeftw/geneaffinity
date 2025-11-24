#!/usr/bin/env python3
"""
pipeline_runner.py
Pure Python RNA-RNA sliding-window pipeline runner for Streamlit.
Exposes run_pipeline() at top level.
"""

import pandas as pd
from pipeline import generate_windows, score_candidates, read_fasta_to_dict

def run_pipeline(
    geneA_path,
    geneB_path,
    window_sizes=[35],
    step=1,
    top_k_per_window=10,
    flank=0,
    log_callback=None
):
    """
    Run sliding-window RNA-RNA pipeline (pure Python).
    Tracks nucleotide positions, and returns sequences for each window.
    
    Parameters:
    - geneA_path: FASTA path for query sequence(s)
    - geneB_path: FASTA path for target sequence(s)
    - window_sizes: list of integers
    - step: int, sliding window step
    - top_k_per_window: int, how many top hits to keep per window
    - flank: int, optional flanking region length for scoring
    - log_callback: optional function for logging
    """
    def log(msg):
        if log_callback:
            log_callback(msg)
        else:
            print(msg)

    log("Reading sequences...")
    seqA_dict = read_fasta_to_dict(geneA_path)
    seqB_dict = read_fasta_to_dict(geneB_path)

    all_results = []

    for w in window_sizes:
        log(f"Generating windows of size {w}...")
        df_windows = generate_windows(geneA_path, w, step, log_callback)

        log(f"Scoring candidates for window size {w}...")
        scored_df = score_candidates(df_windows, seqA_dict, seqB_dict, log_callback, flank=flank)

        log(f"Selecting top {top_k_per_window} hits per query...")
        top_hits = scored_df.sort_values('score', ascending=False).groupby('query_id').head(top_k_per_window).reset_index(drop=True)

        # Add actual sequences for query and target based on positions
        log("Extracting sequences for each hit...")
        def extract_query_seq(row):
            full_seq = seqA_dict[row['seq_id']]
            return full_seq[row['q_start']-1 : row['q_end']]

        def extract_target_seq(row):
            full_seq = seqB_dict[row['target_id']]
            return full_seq[row['t_start']-1 : row['t_end']]

        top_hits['query_seq'] = top_hits.apply(extract_query_seq, axis=1)
        top_hits['target_seq'] = top_hits.apply(extract_target_seq, axis=1)

        all_results.append(top_hits)

    final_df = pd.concat(all_results, ignore_index=True)
    log("Pipeline finished successfully!")
    return final_df
