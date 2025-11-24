from pipeline import generate_windows, score_candidates, log, read_fasta_to_dict

def run_pipeline(
    geneA_path,
    geneB_path,
    window_sizes=[35],
    step=5,
    flank=0,
    energy_cutoff_fast=-6.0,
    top_k_per_window=100,
    thermo_mode="B",
    topN_modeC=None,
    log_callback=None
):
    """
    Run pure Python RNA-RNA sliding-window pipeline with optional logging.
    """
    log("Loading target sequences...", log_callback)
    seqB_dict = read_fasta_to_dict(geneB_path)

    all_results = []
    for w in window_sizes:
        log(f"=== Processing window size {w} ===", log_callback)
        df_windows = generate_windows(geneA_path, w, step, log_callback)
        log(f"Scoring candidates for window size {w}...", log_callback)
        scored_df = score_candidates(df_windows, seqB_dict, log_callback)

        # Keep top K per query
        scored_df = scored_df.sort_values('score', ascending=False).groupby('query_id').head(top_k_per_window).reset_index(drop=True)
        log(f"Top {top_k_per_window} hits selected per query window", log_callback)

        # Optional Mode C: no real thermodynamics, placeholder for demonstration
        if thermo_mode=="C" and topN_modeC:
            log(f"Mode C: further rescoring top {topN_modeC} hits per query...", log_callback)
            # Here you could add a more expensive scoring function
            # For now, we just keep the topN_modeC hits
            scored_df = scored_df.groupby('query_id').head(topN_modeC).reset_index(drop=True)

        all_results.append(scored_df)

    final_df = pd.concat(all_results, ignore_index=True)
    log(f"Pipeline completed. Total hits: {len(final_df)}", log_callback)
    return final_df
