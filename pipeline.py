#!/usr/bin/env python3
"""
pipeline.py
General RNA-RNA sliding-window interaction pipeline using IntaRNA for both
fast scanning (replacement for RiSearch) and final rescoring.

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

# ---------- Fast search (IntaRNA replacement for RiSearch) ----------
def call_intarna_fast(query_fa, target_fa, out_tsv, energy_cutoff=None, max_hits=None,
                      seedTQ=3, seedBP=8, threads=1, intarna_bin=None):
    """Run IntaRNA in 'fast scan' mode to produce candidate interactions.

    Uses: --outMode C and per-user seed params (seedTQ, seedBP).
    The 'moderate' preset requested is seedTQ=3, seedBP=8 by default.
    """
    if intarna_bin is None:
        intarna_bin = os.path.join(os.path.dirname(__file__), "bin", "IntaRNA")
    out_tsv = str(out_tsv)
    cmd = [
        intarna_bin,
        "--query", str(query_fa),
        "--target", str(target_fa),
        "--outMode", "C",
        "--out", out_tsv,
        "--threads", str(threads),
        "--seedTQ", str(seedTQ),
        "--seedBP", str(seedBP)
    ]
    # optional filters (IntaRNA supports --energyThreshold but name can vary by version)
    if energy_cutoff is not None:
        # IntaRNA uses --E or --energyThreshold depending on version. We'll attempt --energyThreshold.
        cmd += ["--energyThreshold", str(energy_cutoff)]
    # Run
    print("Running IntaRNA fast search:", " ".join(cmd))
    try:
        subprocess.run(cmd, check=True)
    except Exception as e:
        raise RuntimeError(f"IntaRNA (fast) failed: {e}")
    return out_tsv

def parse_intarna_fast_output(tsv_path):
    """
    Parse IntaRNA --outMode C output and return a DataFrame with at least:
      ['query_id','target_id','q_start','q_end','t_start','t_end','energy']
    IntaRNA column names vary; we detect common names and normalize.
    """
    # Try reading as tab or comma separated, skipping commented lines
    with open(tsv_path, "r") as fh:
        # try to detect delimiter by first non-comment line
        first_line = None
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            first_line = line
            break
    if first_line is None:
        return pd.DataFrame(columns=['query_id','target_id','q_start','q_end','t_start','t_end','energy'])
    delimiter = "\t" if "\t" in first_line else ","
    df = pd.read_csv(tsv_path, sep=delimiter, comment="#", dtype=str)
    # Common column name variants mapping to canonical names
    col_map = {}
    cols = [c.lower() for c in df.columns]
    for c in df.columns:
        lc = c.lower()
        if 'query' in lc and 'id' in lc:
            col_map[c] = 'query_id'
        elif lc in ('qstart', 'q_start', 'q_startpos', 'qstartpos') or ('qstart' in lc and 'start' in lc):
            col_map[c] = 'q_start'
        elif lc in ('qend', 'q_end', 'q_endpos', 'qendpos') or ('qend' in lc and 'end' in lc):
            col_map[c] = 'q_end'
        elif 'target' in lc and 'id' in lc:
            col_map[c] = 'target_id'
        elif lc in ('tstart', 't_start', 't_startpos', 'tstartpos') or ('tstart' in lc and 'start' in lc):
            col_map[c] = 't_start'
        elif lc in ('tend', 't_end', 't_endpos', 'tendpos') or ('tend' in lc and 'end' in lc):
            col_map[c] = 't_end'
        elif lc in ('e', 'energy', 'freeenergy', 'interactionenergy') or 'energy' in lc:
            col_map[c] = 'energy'
        # else ignore for now
    # rename found columns
    df = df.rename(columns=col_map)
    # If some canonical columns are absent but present under other patterns, attempt more heuristics:
    # e.g., IntaRNA sometimes uses columns: qAcc, qStart, qEnd, tAcc, tStart, tEnd, E
    heuristic_map = {}
    for c in df.columns:
        if c.lower().startswith('q') and 'acc' in c.lower() and 'query_id' not in df.columns:
            heuristic_map[c] = 'query_id'
        if c.lower().startswith('t') and 'acc' in c.lower() and 'target_id' not in df.columns:
            heuristic_map[c] = 'target_id'
    if heuristic_map:
        df = df.rename(columns=heuristic_map)
    # Ensure canonical columns exist (create if missing with NaNs)
    for rc in ['query_id','target_id','q_start','q_end','t_start','t_end','energy']:
        if rc not in df.columns:
            df[rc] = pd.NA
    # coerce numerics
    for c in ['q_start','q_end','t_start','t_end','energy']:
        df[c] = pd.to_numeric(df[c], errors='coerce')
    # reorder to canonical
    df = df[['query_id','target_id','q_start','q_end','t_start','t_end','energy'] + [c for c in df.columns if c not in
                                                                                       ['query_id','target_id','q_start','q_end','t_start','t_end','energy']]]
    return df

# ---------- Extract flanking context ----------
def extract_contexts(target_fa, hits_df, flank, out_dir):
    ensure_dir(out_dir)
    seqs = read_fasta_to_dict(target_fa)
    rows = []
    for idx, row in hits_df.iterrows():
        tid = str(row['target_id'])
        start = int(row['t_start']) if pd.notna(row['t_start']) else None
        end = int(row['t_end']) if pd.notna(row['t_end']) else None
        if start is None or end is None:
            # skip hits without coordinates
            continue
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
    """Run IntaRNA for rescoring (full IntaRNA prediction)."""
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
    # canonicalize energy column if named differently
    if 'E' in df.columns and 'Energy' not in df.columns and 'energy' not in df.columns:
        df = df.rename(columns={'E':'energy'})
    if 'Energy' in df.columns and 'energy' not in df.columns:
        df = df.rename(columns={'Energy':'energy'})
    # keep df as-is; the aggregator will look for 'energy' or 'E' entries under intarna_
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
        # determine energy column
        energy_col = None
        for candidate in ['energy', 'E', 'Energy']:
            if candidate in df.columns:
                energy_col = candidate
                break
        if energy_col is None:
            # skip if no energy info
            continue
        best_row = df.loc[df[energy_col].idxmin()].to_dict()
        merged = hits_df.loc[hit_index].to_dict()
        merged.update({f"intarna_{k}": v for k, v in best_row.items()})
        merged['intarna_csv'] = csvp
        rows.append(merged)
    if not rows:
        print("No intarna results parsed.")
        return None
    outdf = pd.DataFrame(rows)
    # find intarna energy column name present
    intarna_energy_cols = [c for c in outdf.columns if c.lower().startswith('intarna_') and ('energy' in c.lower() or c.lower().endswith('_e') or c.lower().endswith('_energy') or c.lower().endswith('e'))]
    # fallback: search for 'intarna_E' or 'intarna_energy'
    if 'intarna_E' in outdf.columns:
        sort_col = 'intarna_E'
    elif 'intarna_energy' in outdf.columns:
        sort_col = 'intarna_energy'
    elif intarna_energy_cols:
        sort_col = intarna_energy_cols[0]
    else:
        # try to find any column that holds numeric intarna energies by name
        candidates = [c for c in outdf.columns if 'intarna' in c.lower() and c.lower().endswith(('e','energy'))]
        sort_col = candidates[0] if candidates else None
    if sort_col:
        outdf = outdf.sort_values(by=sort_col)
    outdf.to_csv(out_tsv, index=False)
    # try to produce Excel as well
    xlsx_path = Path(out_tsv).with_suffix('.xlsx')
    try:
        df_for_xl = outdf.copy()
        df_for_xl.to_excel(xlsx_path, index=False)
        print("Also wrote Excel:", xlsx_path)
    except Exception as e:
        print("Could not write Excel file (openpyxl not installed or other issue):", e)
        print("CSV remains available at", out_tsv)
    return outdf

# ---------- High-level pipeline steps ----------
def run_stage_make_windows(args):
    generate_windows(args.seq, args.window, args.step, args.out)
    print("Wrote windows to", args.out)

def run_stage_fast_search(args):
    # uses IntaRNA as fast-search with moderate seed (seedTQ=3, seedBP=8)
    call_intarna_fast(args.query, args.target, args.out, energy_cutoff=args.energy, max_hits=args.max_hits,
                      seedTQ=3, seedBP=8, threads=1, intarna_bin=getattr(args, 'intarna_bin', None))
    print("Fast search (IntaRNA) complete; output at", args.out)

def run_stage_extract_contexts(args):
    df = parse_intarna_fast_output(args.hits)
    outdf = extract_contexts(args.target, df, args.flank, args.out_dir)
    outdf.to_csv(Path(args.out_dir)/"contexts_summary.tsv", index=False)
    print("Contexts written to", args.out_dir)

def run_stage_rescore(args):
    hits_df = parse_intarna_fast_output(args.hits) if args.hits.endswith(('.tsv','.txt','.csv')) else pd.read_csv(args.hits)
    if args.top_k is not None:
        if 'energy' in hits_df.columns:
            hits_df = hits_df.sort_values('energy').groupby('query_id').head(args.top_k).reset_index(drop=True)
        else:
            hits_df = hits_df.groupby('query_id').head(args.top_k).reset_index(drop=True)
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
        csvp = call_intarna(query_seq, r['context_fasta'], out_pref, threads=args.intarna_threads, intarna_bin=getattr(args,'intarna_bin',None))
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
    energy_cutoff_fast = config.get('energy_cutoff_fast', None)
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
        fast_out = tmp_dir / f"intarna_fast_w{w}.tsv"
        # use IntaRNA fast scan (moderate seeds: seedTQ=3 seedBP=8)
        call_intarna_fast(win_fa, geneB, fast_out, energy_cutoff=energy_cutoff_fast, max_hits=100000,
                          seedTQ=3, seedBP=8, threads=1)
        try:
            hits_df = parse_intarna_fast_output(fast_out)
        except Exception as e:
            print("Failed to parse IntaRNA fast output:", e)
            continue
        # filter top_k per query by reported (fast) energy if present
        if 'energy' in hits_df.columns:
            hits_df = hits_df.sort_values('energy').groupby('query_id').head(top_k).reset_index(drop=True)
        else:
            hits_df = hits_df.groupby('query_id').head(top_k).reset_index(drop=True)
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
        # attempt to find intarna energy column for statistics
        ecols = [c for c in agg_df.columns if 'intarna' in c.lower() and 'energy' in c.lower()]
        if ecols:
            mean_E = float(pd.to_numeric(agg_df[ecols[0]], errors='coerce').mean())
            best_E = float(pd.to_numeric(agg_df[ecols[0]], errors='coerce').min())
        else:
            mean_E = float('nan')
            best_E = float('nan')
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
    p.add_argument('--energy_cutoff_fast', type=float, default=None)
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
    p.add_argument('--energy', type=float, default=None)
    p.add_argument('--max-hits', type=int, default=1000)
    p.add_argument('--intarna_bin', required=False, default=None)

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
    p.add_argument('--intarna_bin', required=False, default=None)

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
