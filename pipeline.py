#!/usr/bin/env python3
"""
pipeline.py
General RNA-RNA sliding-window interaction pipeline.

Usage examples:
  # Full sweep (reads config.yaml)
  python pipeline.py tune --config config.yaml

  # Make windows only
  python pipeline.py make-windows --seq data/MALAT1.fa --window 35 --step 5 --out malat1_windows.fa

  # Create a config YAML
  python pipeline.py create-config --geneA data/MALAT1.fa --geneB data/targets.fa --out config.yaml
"""

import os
import argparse
import subprocess
from pathlib import Path
import yaml
from Bio import SeqIO
import pandas as pd
from tqdm import tqdm

# ---------- Utilities ----------
def ensure_dir(path):
    Path(path).mkdir(parents=True, exist_ok=True)

def read_fasta_to_dict(fasta_path):
    return {rec.id: str(rec.seq) for rec in SeqIO.parse(str(fasta_path), "fasta")}

def write_fasta(records, path):
    with open(path, "w") as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n{seq}\n")

# ---------- Window generation ----------
def generate_windows(fasta_path, window, step, out_fa):
    seqs = read_fasta_to_dict(fasta_path)
    records = []
    for rid, seq in seqs.items():
        L = len(seq)
        for start in range(0, max(1, L - window + 1), step):
            s = start
            e = start + window
            subseq = seq[s:e]
            wid = f"{rid}|{s+1}-{e}"  # 1-based coords
            records.append((wid, subseq))
    write_fasta(records, out_fa)
    return out_fa

# ---------- Fast search (risearch wrapper) ----------
def call_risearch(query_fa, target_fa, out_tsv, energy_cutoff=-6.0, max_hits=1000, risearch_bin=None):
    """Run RiSearch2 using the provided binary path or default 'risearch2' in the same folder."""
    if risearch_bin is None:
        risearch_bin = os.path.join(os.path.dirname(__file__),"bin", "risearch2")
    cmd = [
        risearch_bin,
        "-q", str(query_fa),
        "-t", str(target_fa),
        "--energy-cutoff", str(energy_cutoff),
        "--max-hits", str(max_hits),
        "-o", str(out_tsv)
    ]
    print("Running:", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True)
    except Exception as e:
        raise RuntimeError(f"risearch2 failed: {e}")
    return out_tsv

def parse_risearch_tsv(tsv_path):
    df = pd.read_csv(tsv_path, sep="\t", comment="#", dtype=str)
    expected_cols = ['query_id', 'target_id', 'q_start', 'q_end', 't_start', 't_end', 'energy']
    if not all(c in df.columns for c in expected_cols):
        raise RuntimeError(f"risearch2 output columns do not include expected: {df.columns.tolist()}")
    for c in ['q_start', 'q_end', 't_start', 't_end', 'energy']:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    return df

# ---------- Extract flanking context ----------
def extract_contexts(target_fa, hits_df, flank, out_dir):
    ensure_dir(out_dir)
    seqs = read_fasta_to_dict(target_fa)
    rows = []
    for idx, row in hits_df.iterrows():
        tid = str(row['target_id'])
        start = int(row['t_start'])
        end = int(row['t_end'])
        seq = seqs.get(tid)
        if seq is None:
            raise RuntimeError(f"target id {tid} not found in {target_fa}")
        L = len(seq)
        s = max(0, start - 1 - flank)
        e = min(L, end + flank)
        context_seq = seq[s:e]
        fid = f"{tid}|{start}-{end}|flank{flank}"
        outp = Path(out_dir) / f"{fid}.fa"
        write_fasta([(fid, context_seq)], outp)
        rows.append({
            'hit_index': idx,
            'target_id': tid,
            't_start': start,
            't_end': end,
            'context_fasta': str(outp),
        })
    return pd.DataFrame(rows)

# ---------- Rescore with IntaRNA ----------
def call_intarna(query_fa, target_context_fa, out_prefix, threads=1, intarna_bin=None):
    """Run IntaRNA using the provided binary path or default 'IntaRNA' in the same folder."""
    if intarna_bin is None:
        intarna_bin = os.path.join(os.path.dirname(__file__), "bin", "IntaRNA")
    out_prefix = str(out_prefix)
    out_csv = out_prefix + ".csv"
    cmd = [
        intarna_bin,
        "--query", str(query_fa),
        "--target", str(target_context_fa),
        "--outMode", "C",
        "--out", out_csv,
        "--threads", str(threads)
    ]
    try:
        subprocess.run(cmd, check=True)
    except Exception as e:
        raise RuntimeError(f"IntaRNA failed: {e}")
    return out_csv

def parse_intarna_csv(csv_path):
    df = pd.read_csv(csv_path, comment="#")
    if 'E' not in df.columns and 'Energy' in df.columns:
        df = df.rename(columns={'Energy':'E'})
    return df

# ---------- Aggregation & ranking ----------
def aggregate_and_rank(intarna_results, hits_df, out_tsv):
    rows = []
    for hit_index, csvp in intarna_results:
        try:
            df = parse_intarna_csv(csvp)
        except Exception:
            continue
        if df.shape[0] == 0:
            continue
        best_row = df.loc[df['E'].idxmin()].to_dict()
        merged = hits_df.loc[hit_index].to_dict()
        merged.update({f"intarna_{k}": v for k, v in best_row.items()})
        merged['intarna_csv'] = csvp
        rows.append(merged)
    if not rows:
        print("No intarna results parsed.")
        return None
    outdf = pd.DataFrame(rows)
    outdf = outdf.sort_values(by='intarna_E')
    outdf.to_csv(out_tsv, index=False)
    return outdf

# ---------- High-level pipeline steps ----------
def run_stage_make_windows(args):
    generate_windows(args.seq, args.window, args.step, args.out)
    print("Wrote windows to", args.out)

def run_stage_fast_search(args):
    call_risearch(args.query, args.target, args.out, energy_cutoff=args.energy, max_hits=args.max_hits)
    print("Fast search complete; output at", args.out)

def run_stage_extract_contexts(args):
    df = pd.read_csv(args.hits, sep="\t", comment="#")
    outdf = extract_contexts(args.target, df, args.flank, args.out_dir)
    outdf.to_csv(Path(args.out_dir)/"contexts_summary.tsv", index=False)
    print("Contexts written to", args.out_dir)

def run_stage_rescore(args):
    hits_df = pd.read_csv(args.hits)
    if args.top_k is not None:
        hits_df = hits_df.sort_values('energy').groupby('query_id').head(args.top_k).reset_index(drop=True)
    contexts_df = extract_contexts(args.target, hits_df, args.flank, args.context_dir)
    intarna_results = []
    for idx, r in tqdm(contexts_df.iterrows(), total=len(contexts_df)):
        qid = hits_df.loc[int(r['hit_index']), 'query_id']
        query_seq = None
        if not args.query_windows:
            raise RuntimeError("query_windows required for rescoring stage")
        for rec in SeqIO.parse(args.query_windows, "fasta"):
            if rec.id == qid:
                tmpq = Path(args.context_dir) / f"tmp_query_{qid}.fa"
                write_fasta([(qid, str(rec.seq))], tmpq)
                query_seq = tmpq
                break
        if query_seq is None:
            raise RuntimeError(f"Query id {qid} not found in {args.query_windows}")
        out_pref = Path(args.context_dir) / f"intarna_{idx}"
        csvp = call_intarna(query_seq, r['context_fasta'], out_pref, threads=args.intarna_threads)
        intarna_results.append((int(r['hit_index']), csvp))
        try:
            os.remove(query_seq)
        except OSError:
            pass
    agg = aggregate_and_rank(intarna_results, hits_df, args.out)
    print("Rescoring done. Aggregated output:", args.out)
    return agg

# ---------- Tuning ----------
def tune_windows(config):
    geneA = config['geneA']
    geneB = config['geneB']
    window_sizes = config['window_sizes']
    step = config.get('step', 5)
    flank = config.get('flank', 100)
    energy_cutoff_fast = config.get('energy_cutoff_fast', -6)
    top_k = config.get('top_k_per_window', 100)
    intarna_threads = config.get('intarna_threads', 2)
    tmp_dir = Path(config.get('tmp_dir', 'tmp'))
    out_dir = Path(config.get('out_dir', 'results'))
    ensure_dir(tmp_dir)
    ensure_dir(out_dir)
    summary_rows = []
    for w in window_sizes:
        print(f"=== Window size {w} ===")
        win_fa = tmp_dir / f"windows_w{w}.fa"
        generate_windows(geneA, w, step, win_fa)
        fast_out = tmp_dir / f"risearch_w{w}.tsv"
        # use local binary
        call_risearch(win_fa, geneB, fast_out, energy_cutoff=energy_cutoff_fast, max_hits=100000)
        try:
            hits_df = parse_risearch_tsv(fast_out)
        except Exception as e:
            print("Failed to parse risearch output:", e)
            continue
        hits_df = hits_df.sort_values('energy').groupby('query_id').head(top_k).reset_index(drop=True)
        contexts_dir = tmp_dir / f"contexts_w{w}"
        contexts_dir.mkdir(parents=True, exist_ok=True)
        contexts_df = extract_contexts(geneB, hits_df, flank, contexts_dir)
        intarna_results = []
        for idx, r in contexts_df.iterrows():
            qid = hits_df.loc[int(r['hit_index']), 'query_id']
            qtmp = contexts_dir / f"query_{qid}.fa"
            for rec in SeqIO.parse(str(win_fa), "fasta"):
                if rec.id == qid:
                    write_fasta([(qid, str(rec.seq))], qtmp)
                    break
            out_pref = contexts_dir / f"intarna_{idx}"
            try:
                csvp = call_intarna(qtmp, r['context_fasta'], out_pref, threads=intarna_threads)
            except Exception as e:
                print("IntaRNA failed for", qtmp, r['context_fasta'], e)
                continue
            intarna_results.append((int(r['hit_index']), csvp))
            try:
                os.remove(qtmp)
            except Exception:
                pass
        agg_out = out_dir / f"aggregated_w{w}.csv"
        agg_df = aggregate_and_rank(intarna_results, hits_df, agg_out)
        if agg_df is None:
            continue
        n_hits = len(agg_df)
        mean_E = float(agg_df['intarna_E'].mean())
        best_E = float(agg_df['intarna_E'].min())
        summary_rows.append({'window': w, 'n_hits': n_hits, 'mean_intarna_E': mean_E, 'best_intarna_E': best_E, 'agg_csv': str(agg_out)})
        print(f"Window {w}: hits={n_hits}, meanE={mean_E:.2f}, bestE={best_E:.2f}")
    summary_df = pd.DataFrame(summary_rows)
    summary_df.to_csv(out_dir / "tuning_summary.csv", index=False)
    print("Tuning complete. Summary at", out_dir / "tuning_summary.csv")
    return summary_df

# ---------- Create config ----------
def create_config(args):
    cfg = {
        'geneA': args.geneA,
        'geneB': args.geneB,
        'window_sizes': args.window_sizes,
        'step': args.step,
        'flank': args.flank,
        'energy_cutoff_fast': args.energy_cutoff_fast,
        'top_k_per_window': args.top_k_per_window,
        'threads': args.threads,
        'intarna_threads': args.intarna_threads,
        'tmp_dir': args.tmp_dir,
        'out_dir': args.out_dir
    }
    with open(args.out, 'w') as fh:
        yaml.dump(cfg, fh)
    print(f"Config YAML written to {args.out}")

# ---------- CLI ----------
def main():
    parser = argparse.ArgumentParser()
    sub = parser.add_subparsers(dest='cmd')

    # create-config
    p = sub.add_parser('create-config')
    p.add_argument('--geneA', required=True)
    p.add_argument('--geneB', required=True)
    p.add_argument('--out', default='config.yaml')
    p.add_argument('--window_sizes', nargs='+', type=int, default=[25,35,45])
    p.add_argument('--step', type=int, default=5)
    p.add_argument('--flank', type=int, default=100)
    p.add_argument('--energy_cutoff_fast', type=float, default=-6)
    p.add_argument('--top_k_per_window', type=int, default=50)
    p.add_argument('--threads', type=int, default=4)
    p.add_argument('--intarna_threads', type=int, default=2)
    p.add_argument('--tmp_dir', default='tmp')
    p.add_argument('--out_dir', default='results')

    # make-windows
    p = sub.add_parser('make-windows')
    p.add_argument('--seq', required=True)
    p.add_argument('--window', type=int, required=True)
    p.add_argument('--step', type=int, default=5)
    p.add_argument('--out', required=True)

    # fast-search
    p = sub.add_parser('fast-search')
    p.add_argument('--query', required=True)
    p.add_argument('--target', required=True)
    p.add_argument('--out', required=True)
    p.add_argument('--energy', type=float, default=-6.0)
    p.add_argument('--max-hits', type=int, default=1000)

    # extract-contexts
    p = sub.add_parser('extract-contexts')
    p.add_argument('--target', required=True)
    p.add_argument('--hits', required=True)
    p.add_argument('--flank', type=int, required=True)
    p.add_argument('--out_dir', required=True)

    # rescore
    p = sub.add_parser('rescore')
    p.add_argument('--query_windows', required=True)
    p.add_argument('--target', required=True)
    p.add_argument('--hits', required=True)
    p.add_argument('--flank', type=int, default=100)
    p.add_argument('--context_dir', required=True)
    p.add_argument('--intarna_threads', type=int, default=2)
    p.add_argument('--top_k', type=int, default=None)
    p.add_argument('--out', required=True)

    # tune
    p = sub.add_parser('tune')
    p.add_argument('--config', required=True)

    args = parser.parse_args()

    if args.cmd == 'create-config':
        create_config(args)
    elif args.cmd == 'make-windows':
        run_stage_make_windows(args)
    elif args.cmd == 'fast-search':
        run_stage_fast_search(args)
    elif args.cmd == 'extract-contexts':
        run_stage_extract_contexts(args)
    elif args.cmd == 'rescore':
        run_stage_rescore(args)
    elif args.cmd == 'tune':
        with open(args.config) as fh:
            cfg = yaml.safe_load(fh)
        tune_windows(cfg)
    else:
        parser.print_help()

if __name__ == "__main__":
    main()
