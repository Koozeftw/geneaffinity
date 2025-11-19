import os
from pathlib import Path
import pandas as pd
from pipeline import (
    generate_windows,
    call_intarna_fast,     # use IntaRNA for fast scan
    extract_contexts,
    call_intarna,
    aggregate_and_rank,
)
from Bio import SeqIO

def run_pipeline(
    geneA_path,
    geneB_path,
    window_sizes=[35],
    step=5,
    flank=100,
    energy_cutoff_fast=-6.0,
    top_k_per_window=100,
    threads=1,
    intarna_bin=None,
    log_callback=None,
):
    """
    Run the RNA-RNA sliding-window pipeline using IntaRNA
    for both fast scanning and full rescoring.
    """

    def log(msg):
        if log_callback:
            log_callback(msg)
        else:
            print(msg)

    if intarna_bin is None:
        intarna_bin = "/Users/colekuznitz/miniforge3/bin/IntaRNA"

    all_results = []

    for w in window_sizes:
        log(f"=== Window size {w} ===")

        # 1️⃣ Generate sliding windows
        win_fa = Path(f"tmp_windows_w{w}.fa")
        generate_windows(geneA_path, w, step, win_fa)
        log(f"Generated {len(list(SeqIO.parse(win_fa, 'fasta')))} windows for geneA.")

        # 2️⃣ Fast search with IntaRNA
        fast_out = Path(f"tmp_intarna_fast_w{w}.tsv")
        call_intarna_fast(win_fa, geneB_path, fast_out, energy_cutoff=energy_cutoff_fast,
                          threads=threads, intarna_bin=intarna_bin)
        log(f"IntaRNA fast scan completed. Output: {fast_out}")

        # 3️⃣ Parse hits and select top K per query
        hits_df = pd.read_csv(fast_out, sep="\t", comment="#")
        hits_df = hits_df.sort_values("energy").groupby("query_id").head(top_k_per_window).reset_index(drop=True)
        log(f"Selected top {top_k_per_window} hits per query.")

        # 4️⃣ Extract contexts with flanking regions
        contexts_dir = Path(f"tmp_contexts_w{w}")
        contexts_dir.mkdir(parents=True, exist_ok=True)
        contexts_df = extract_contexts(geneB_path, hits_df, flank, contexts_dir)
        log(f"Extracted {len(contexts_df)} context sequences with flank={flank}.")

        # 5️⃣ Rescore with IntaRNA
        intarna_results = []
        for idx, row in contexts_df.iterrows():
            qid = hits_df.loc[int(row["hit_index"]), "query_id"]
            tmp_query_fa = contexts_dir / f"tmp_query_{qid}.fa"

            # extract query sequence
            for rec in SeqIO.parse(win_fa, "fasta"):
                if rec.id == qid:
                    with open(tmp_query_fa, "w") as fh:
                        fh.write(f">{rec.id}\n{str(rec.seq)}\n")
                    break

            out_pref = contexts_dir / f"intarna_{idx}"
            csv_file = call_intarna(tmp_query_fa, row["context_fasta"], out_pref,
                                    threads=threads, intarna_bin=intarna_bin)
            intarna_results.append((int(row["hit_index"]), csv_file))
            log(f"IntaRNA rescoring done for hit {idx}")

            tmp_query_fa.unlink(missing_ok=True)

        # 6️⃣ Aggregate & rank results
        aggregated_df = aggregate_and_rank(intarna_results, hits_df, contexts_dir / f"aggregated_w{w}.csv")
        log(f"Aggregated results written for window {w}.")
        all_results.append(aggregated_df)

    final_df = pd.concat(all_results, ignore_index=True)
    log("Pipeline finished successfully.")
    return final_df
