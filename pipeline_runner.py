import pandas as pd
from pipeline import generate_windows, call_risearch, parse_risearch_tsv, extract_contexts, call_intarna, aggregate_and_rank
from pathlib import Path

def run_pipeline(
    geneA_path,
    geneB_path,
    window_sizes=[35],
    step=5,
    flank=100,
    energy_cutoff_fast=-6.0,
    top_k_per_window=100,
    threads=1,
    log_callback=print
):
    tmp_dir = Path("tmp")
    tmp_dir.mkdir(exist_ok=True)
    results_all = []

    for w in window_sizes:
        log_callback(f"=== Running window size {w} ===")
        # Generate windows
        win_fa = tmp_dir / f"windows_w{w}.fa"
        generate_windows(geneA_path, w, step, win_fa)

        # Fast search
        fast_out = tmp_dir / f"risearch_w{w}.tsv"
        call_risearch(win_fa, geneB_path, fast_out, energy_cutoff=energy_cutoff_fast, max_hits=100000)

        # Parse results
        hits_df = parse_risearch_tsv(fast_out)
        hits_df = hits_df.sort_values('energy').groupby('query_id').head(top_k_per_window).reset_index(drop=True)

        # Extract flanking context
        contexts_dir = tmp_dir / f"contexts_w{w}"
        contexts_dir.mkdir(exist_ok=True)
        contexts_df = extract_contexts(geneB_path, hits_df, flank, contexts_dir)

        # Rescore with IntaRNA
        intarna_results = []
        for idx, r in contexts_df.iterrows():
            qid = hits_df.loc[int(r['hit_index']), 'query_id']
            qtmp = contexts_dir / f"query_{qid}.fa"
            for rec in open(win_fa):
                pass  # you can use SeqIO if needed
            out_pref = contexts_dir / f"intarna_{idx}"
            csvp = call_intarna(qtmp, r['context_fasta'], out_pref, threads=threads)
            intarna_results.append((int(r['hit_index']), csvp))

        # Aggregate
        agg_out = tmp_dir / f"aggregated_w{w}.csv"
        agg_df = aggregate_and_rank(intarna_results, hits_df, agg_out)
        results_all.append(agg_df)

    return pd.concat(results_all, ignore_index=True)

