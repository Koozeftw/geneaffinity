import os
from pathlib import Path
import pandas as pd
from pipeline import (
    generate_windows,
    call_vienna_fast,          # <-- NEW: replaces call_risearch
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
    vienna_bin=None,
    intarna_bin=None,
    log_callback=None,
):
    """
    Run the full RNA-RNA sliding-window pipeline using ViennaRNA
    for the fast screen instead of RiSearch.

    Parameters
    ----------
    geneA_path : str
        FASTA path for query sequences.
    geneB_path : str
        FASTA path for target sequences.
    window_sizes : list[int]
        Window sizes to use.
    step : int
        Step size for sliding window.
    flank : int
        Flanking sequence length for rescoring.
    energy_cutoff_fast : float
        Energy cutoff for the ViennaRNA fast screen.
    top_k_per_window : int
        Number of top hits per query to rescore.
    threads : int
        Number of threads for IntaRNA.
    vienna_bin : str
        Path to RNAduplex binary.
    intarna_bin : str
        Path to IntaRNA binary.
    log_callback : callable
        Optional callback for logging.
    """

    def log(msg):
        if log_callback:
            log_callback(msg)
        else:
            print(msg)

    if vienna_bin is None:
        vienna_bin = os.path.join(os.path.dirname(__file__), "RNAduplex")
    if intarna_bin is None:
        intarna_bin = os.path.join(os.path.dirname(__file__), "IntaRNA")

    all_results = []

    for w in window_sizes:
        log(f"=== Window size {w} ===")

        # 1️⃣ Generate sliding windows
        win_fa = Path(f"tmp_windows_w{w}.fa")
        generate_windows(geneA_path, w, step, win_fa)
