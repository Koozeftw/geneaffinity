# pipeline_runner.py
"""
pipeline_runner.py

High-level runner that orchestrates the pure-Python pipeline (Option B default,
optional Mode C re-evaluation). Designed to be called from Streamlit app.py.

Expected functions to exist in pipeline.py:
  - generate_windows(fasta_path, window, step, out_fa)
  - fast_scan_all(query_fa, target_fa, seed_len=..., min_seed_hits=..., ...)
  - rescore_candidates(candidates, query_fa, target_fa, mode='B', top_n_for_C=20)
  - read_fasta_to_dict(path)
  - write_fasta(records, path)      (optional convenience)
  - write_candidates_tsv(candidates, out_tsv)  (optional)
  - write_rescored_csv(results, out_csv)      (optional)

This runner focuses on wiring those pieces together and returning a pandas DataFrame.
"""

from pathlib import Path
import tempfile
import os
import pandas as pd
import csv
from typing import List, Dict, Optional

# Import pure-Python pipeline functions
from pipeline import (
    generate_windows,
    fast_scan_all,
    rescore_candidates,
    read_fasta_to_dict,
    write_fasta,
    write_candidates_tsv,
    write_rescored_csv,
)

def _log_maybe(callback, msg: str):
    if callback:
        try:
            callback(msg)
        except Exception:
            # fail silently for logging errors
            print("[logger error]", msg)
    else:
        print(msg)

def _select_top_k_per_query(candidates: List[Dict], k: int) -> List[Dict]:
    """
    candidates: list of dicts with 'query_id', 'target_id', and 'approx_energy' (or similar)
    returns filtered list keeping at most k candidates per query_id (best approx_energy)
    """
    if k is None or k <= 0:
        return candidates
    # group by query_id
    grouped = {}
    for c in candidates:
        q = c.get('query_id')
        grouped.setdefault(q, []).append(c)
    out = []
    for q, lst in grouped.items():
        # sort by approx_energy (more negative/better first). Fallback to 0 if not present.
        lst_sorted = sorted(lst, key=lambda x: x.get('approx_energy', 0.0))
        out.extend(lst_sorted[:k])
    return out

def run_pipeline(
    geneA_path: str,
    geneB_path: str,
    window_sizes: List[int] = [35],
    step: int = 5,
    flank: int = 100,
    energy_cutoff_fast: float = -6.0,
    top_k_per_window: int = 100,
    threads: int = 1,                # placeholder — threading not implemented in this runner
    thermo_mode: str = "B",          # "B" or "C"
    topN_modeC: int = 20,            # when thermo_mode == "C", re-evaluate top N from B with Mode C
    seed_len: int = 8,
    min_seed_hits: int = 1,
    max_hits_per_query: Optional[int] = None,
    log_callback=None
) -> pd.DataFrame:
    """
    High-level runner.

    Parameters
    ----------
    geneA_path : path to FASTA for gene A (will be windowed)
    geneB_path : path to FASTA for gene B (targets to scan against)
    window_sizes : list of window sizes to generate from geneA
    step : sliding-window step
    flank : (unused here for rescoring) kept for API compatibility
    energy_cutoff_fast : approximate energy threshold for including candidate seeds
    top_k_per_window : keep top-K candidates per query after fast scan
    thermo_mode : 'B' (default) or 'C' (run B, then re-evaluate topN with C)
    topN_modeC : number of top-B hits to re-evaluate with Mode C
    seed_len, min_seed_hits, max_hits_per_query : pass-through to fast scan
    log_callback : optional function to receive log messages (str)
    """
    _log_maybe(log_callback, "Starting pipeline (pure-Python mode).")
    all_results = []

    # Make a small temporary working directory so intermediate files don't clutter repo
    with tempfile.TemporaryDirectory() as workdir:
        workdir = Path(workdir)
        _log_maybe(log_callback, f"Working directory: {workdir}")

        # Ensure input files exist
        geneA_path = Path(geneA_path)
        geneB_path = Path(geneB_path)
        if not geneA_path.exists():
            raise FileNotFoundError(f"geneA_path not found: {geneA_path}")
        if not geneB_path.exists():
            raise FileNotFoundError(f"geneB_path not found: {geneB_path}")

        # Load gene B sequences once (used by rescoring)
        _log_maybe(log_callback, "Reading target FASTA (geneB)...")
        targets_dict = read_fasta_to_dict(str(geneB_path))

        for w in window_sizes:
            _log_maybe(log_callback, f"=== Window size {w} ===")

            # 1) Generate windows from geneA
            windows_fa = workdir / f"windows_w{w}.fa"
            _log_maybe(log_callback, f"Generating windows (window={w}, step={step})...")
            generate_windows(str(geneA_path), w, step, str(windows_fa))

            # Read how many windows created
            windows = read_fasta_to_dict(str(windows_fa))
            n_windows = len(windows)
            _log_maybe(log_callback, f"Generated {n_windows} windows (saved to {windows_fa}).")

            # 2) Fast seed-based scan (Option B fast stage)
            _log_maybe(log_callback, f"Running fast seed scan (seed_len={seed_len})...")
            candidates = fast_scan_all(
                query_fa=str(windows_fa),
                target_fa=str(geneB_path),
                seed_len=seed_len,
                min_seed_hits=min_seed_hits,
                max_hits_per_query=max_hits_per_query,
                approx_score_threshold=energy_cutoff_fast
            )

            _log_maybe(log_callback, f"Fast scan found {len(candidates)} candidate hits before filtering.")

            # 3) Keep top_k_per_window per query
            filtered = _select_top_k_per_query(candidates, top_k_per_window)
            _log_maybe(log_callback, f"Filtered to {len(filtered)} candidates after keeping top {top_k_per_window} per query.")

            # Optional: write candidates TSV for inspection (in workdir)
            cand_tsv = workdir / f"candidates_w{w}.tsv"
            try:
                write_candidates_tsv(filtered, str(cand_tsv))
                _log_maybe(log_callback, f"Wrote candidate TSV to {cand_tsv}")
            except Exception:
                # If write helper not present or fails, write a minimal TSV
                _log_maybe(log_callback, "Falling back to manual candidate TSV writer.")
                if filtered:
                    keys = sorted(filtered[0].keys())
                    with open(cand_tsv, 'w', newline='') as fh:
                        writer = csv.DictWriter(fh, fieldnames=keys, delimiter='\t')
                        writer.writeheader()
                        for r in filtered:
                            writer.writerow({k: r.get(k, '') for k in keys})
                    _log_maybe(log_callback, f"Wrote candidate TSV to {cand_tsv}")

            # 4) Rescore with Mode B (moderate DP) for all filtered candidates
            _log_maybe(log_callback, "Rescoring all candidates with Mode B (moderate DP)...")
            rescored_results = rescore_candidates(
                filtered,
                query_fa=str(windows_fa),
                target_fa=str(geneB_path),
                mode='B',
                top_n_for_C=topN_modeC
            )
            _log_maybe(log_callback, f"Mode B rescoring completed for {len(rescored_results)} candidates.")

            # 5) If thermo_mode == 'C', do Mode C on top-N results (rescore_candidates can embed C but we will handle merging)
            if thermo_mode.upper() == 'C':
                _log_maybe(log_callback, f"Mode C requested: re-evaluating top {topN_modeC} hits from Mode B with higher-accuracy DP...")
                # rescore_candidates already supports mode='C' and will annotate top-N; if it didn't, we could call affine_gap separately.
                # For safety, call rescore_candidates again with mode='C' so that 'rescore_energy_C' fields are present.
                rescored_results = rescore_candidates(
                    filtered,
                    query_fa=str(windows_fa),
                    target_fa=str(geneB_path),
                    mode='C',
                    top_n_for_C=topN_modeC
                )
                _log_maybe(log_callback, "Mode C re-evaluation completed.")

            # 6) Convert rescored_results (list of dict) into pandas DataFrame
            if rescored_results:
                df = pd.DataFrame(rescored_results)
            else:
                df = pd.DataFrame(columns=[
                    "query_id","target_id","q_pos","t_pos","approx_energy",
                    "rescore_energy_B","q_start_B","q_end_B","t_start_B","t_end_B",
                    "rescore_energy_C","q_start_C","q_end_C","t_start_C","t_end_C"
                ])

            # Save a CSV of rescored results for this window size
            out_csv = workdir / f"rescored_w{w}.csv"
            try:
                write_rescored_csv(rescored_results, str(out_csv))
                _log_maybe(log_callback, f"Wrote rescored CSV to {out_csv}")
            except Exception:
                # fallback to pandas to write
                df.to_csv(out_csv, index=False)
                _log_maybe(log_callback, f"Wrote rescored CSV to {out_csv} (pandas fallback)")

            all_results.append(df)

        # 7) concatenate across window sizes
        if all_results:
            final_df = pd.concat(all_results, ignore_index=True, sort=False)
        else:
            final_df = pd.DataFrame()

        _log_maybe(log_callback, "Pipeline finished successfully.")
        return final_df
