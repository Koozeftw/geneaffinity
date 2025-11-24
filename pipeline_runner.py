from pipeline import generate_windows, quick_filter, full_rescore, read_fasta_to_dict
import pandas as pd

def run_pipeline(
    geneA_path,
    geneB_path,
    window_sizes=[35],
    step=1,
    top_k_per_window=10,
    log_callback=None
):
    """Run sliding-window RNA-RNA pipeline (pure Python)"""
    seqA_dict = read_fasta_to_dict(geneA_path)
    seqB_dict = read_fasta_to_dict(geneB_path)

    all_results = []

    for w in window_sizes:
        # 1️⃣ Generate windows
        df_windows = generate_windows(geneA_path, w, step, log_callback)

        # 2️⃣ Quick filter
        top_hits = quick_filter(df_windows, seqB_dict, top_k=top_k_per_window, log_callback=log_callback)

        # 3️⃣ Full rescoring
        rescored_df = full_rescore(top_hits, seqA_dict, seqB_dict, log_callback=log_callback)

        all_results.append(rescored_df)

    final_df = pd.concat(all_results, ignore_index=True)
    return final_df
