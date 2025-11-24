#!/usr/bin/env python3
"""
pipeline_runner.py
Runs the pure Python RNA-RNA pipeline using sliding windows.
"""

from pathlib import Path
import pandas as pd
from pipeline import generate_windows, score_candidates, log, read_fasta_to_dict

def run_pipeline(
    geneA_path,
    geneB_path,
    window_sizes=[35],
    step=5,
    flank=100,
    top_k_per_window=100,
    log_callback=None
):
    """
    Run the sliding-window RNA-RNA pipeline.
    Returns a dataframe with all scored candidates.
    """
    seqA_dict = read_fasta_to_dict(geneA_path)
    seqB_dict = read_fasta_to_dict(geneB_path)

    all_results = []

    for w in window_sizes:
        log(f"=== Window size {w} ===", log_callback)

        # Generate sliding windows for Sequence A
        df_windows = generate_windows(geneA_path, w, step, log_callback)
        log(f"Generated {len(df_windows)} query windows", log_callback)

        # Score candidates against Sequence B
        scored_df = score_candidates(df_windows, seqA_dict, seqB_dict, log_callback, flank)
        log(f"Scored candidates for window size {w}", log_callback)

        # Keep top K hits per query
        top_df = scored_df.sort_values('score', ascending=False).groupby('query_id').head(top_k_per_window).reset_index(drop=True)
        log(f"Selected top {top_k_per_window} hits per query", log_callback)

        all_results.append(top_df)

    final_df = pd.concat(all_results, ignore_index=True)
    log(f"Pipeline finished. Total hits: {len(final_df)}", log_callback)
    return final_df
