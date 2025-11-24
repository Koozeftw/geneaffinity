#!/usr/bin/env python3
"""
pipeline_runner.py
Runs the pure Python RNA-RNA sliding-window pipeline with coordinate tracking.
"""

import pandas as pd
from pathlib import Path
from pipeline import generate_windows, score_candidates, log, read_fasta_to_dict

def run_pipeline(
    geneA_path,
    geneB_path,
    window_sizes=[35],
    step=5,
    flank=0,
    top_k_per_window=100,
    thermo_mode="B",  # B=fast, C=full later
    topN_modeC=None,
    log_callback=None
):
    """
    Run sliding-window RNA-RNA pipeline.
    Returns a dataframe with all scored candidate interactions.
    """
    # Read sequences into dictionaries
    seqA_dict = read_fasta_to_dict(geneA_path)
    seqB_dict = read_fasta_to_dict(geneB_path)

    all_results = []

    for window in window_sizes:
        log(f"=== Window size {window} ===", log_callback)

        # 1️⃣ Generate windows
        df_windows = generate_windows(geneA_path, window, step, log_callback)
        log(f"Generated {len(df_windows)} windows for Sequence A.", log_callback)

        # 2️⃣ Score candidates
        scored_df = score_candidates(df_windows, seqA_dict, seqB_dict, log_callback, flank=flank)
        log(f"Scoring complete. Total candidates: {len(scored_df)}", log_callback)

        # 3️⃣ Select top hits per query window
        top_df = scored_df.sort_values("score", ascending=False).groupby("query_id").head(top_k_per_window).reset_index(drop=True)
        log(f"Selected top {top_k_per_window} hits per query window.", log_callback)

        # 4️⃣ Extract sequences dynamically for display
        top_df['query_seq'] = top_df.apply(lambda row: seqA_dict[row['seq_id']][row['q_start']-1 : row['q_end']], axis=1)
        top_df['target_seq'] = top_df.apply(lambda row: seqB_dict[row['target_id']][row['t_start']-1 : row['t_end']], axis=1)

        all_results.append(top_df)

    # 5️⃣ Combine results from all windows
    final_df = pd.concat(all_results, ignore_index=True)
    log(f"Pipeline finished. Total final hits: {len(final_df)}", log_callback)
    return final_df
