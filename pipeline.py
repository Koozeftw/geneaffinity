#!/usr/bin/env python3
"""
pipeline.py -- Self-contained RNA-RNA sliding-window interaction pipeline
Option B (default): seed-based fast scan + moderate DP rescoring
Optional Option C: higher-accuracy (slower) rescoring of top-N B hits

No external binaries or third-party Python libraries required.
Only standard library used.

Usage examples:
  # Make windows
  python pipeline.py make-windows --seq targets.fa --window 80 --step 10 --out windows.fa

  # Fast scan queries (query.fa contains query sequences, windows.fa target windows)
  python pipeline.py scan --query queries.fa --target windows.fa --out candidates.tsv

  # Rescore using mode B (default)
  python pipeline.py rescore --candidates candidates.tsv --out rescored_b.tsv --mode B

  # Rescore with B then run heavy mode C on top 20 hits
  python pipeline.py rescore --candidates candidates.tsv --out rescored_with_c.tsv --mode B --modeC_top 20
"""

from __future__ import annotations
import argparse
from pathlib import Path
import sys
import math
import csv
from typing import List, Tuple, Dict, Optional

# -------------------------
# Small FASTA utilities (no Biopython)
# -------------------------
def read_fasta_to_dict(path: str) -> Dict[str, str]:
    seqs = {}
    with open(path, 'r') as fh:
        name = None
        lines = []
        for raw in fh:
            line = raw.rstrip('\n')
            if line.startswith('>'):
                if name is not None:
                    seqs[name] = ''.join(lines).replace(' ', '').replace('\r', '')
                name = line[1:].split()[0]
                lines = []
            else:
                lines.append(line.strip())
        if name is not None:
            seqs[name] = ''.join(lines).replace(' ', '').replace('\r', '')
    return seqs

def write_fasta(records: List[Tuple[str, str]], path: str):
    with open(path, 'w') as fh:
        for rid, seq in records:
            fh.write(f">{rid}\n")
            # wrap long lines at 80 chars
            for i in range(0, len(seq), 80):
                fh.write(seq[i:i+80] + "\n")

# -------------------------
# Basic RNA utilities
# -------------------------
def complement(base: str) -> str:
    b = base.upper()
    return {'A':'U','U':'A','T':'A','G':'C','C':'G'}.get(b, 'N')

def is_canonical_pair(a: str, b: str) -> bool:
    a = a.upper(); b = b.upper()
    return (a=='A' and b in ('U','T')) or (a in ('U','T') and b=='A') or (a=='G' and b=='C') or (a=='C' and b=='G') or (a=='G' and b in ('U','T')) or (a in ('U','T') and b=='G')

def pair_energy(a: str, b: str) -> float:
    """Return simple energy (negative favorable) for a base pair.
    Values chosen to be sensible qualitatively (kcal/mol-like scale):
      G-C : -3.0
      A-U : -2.0
      G-U : -1.0
    Mismatch : +1.0 (penalty)
    """
    a = a.upper(); b = b.upper()
    if (a=='G' and b=='C') or (a=='C' and b=='G'):
        return -3.0
    if (a=='A' and b in ('U','T')) or (b=='A' and a in ('U','T')):
        return -2.0
    if (a=='G' and b in ('U','T')) or (b=='G' and a in ('U','T')):
        return -1.0
    # mismatch
    return +1.0

# -------------------------
# Window generation
# -------------------------
def generate_windows(fasta_path: str, window: int, step: int, out_fa: str) -> str:
    seqs = read_fasta_to_dict(fasta_path)
    records = []
    for rid, seq in seqs.items():
        L = len(seq)
        if L < window:
            # still include full sequence as a single "window"
            wid = f"{rid}|1-{L}"
            records.append((wid, seq))
        else:
            for start in range(0, L - window + 1, step):
                s = start
                e = start + window
                subseq = seq[s:e]
                wid = f"{rid}|{s+1}-{e}"  # 1-based coords
                records.append((wid, subseq))
    write_fasta(records, out_fa)
    print(f"[generate_windows] wrote {len(records)} windows to {out_fa}")
    return out_fa

# -------------------------
# Fast scanning (seed-based heuristic)
# -------------------------
def find_perfect_complement_seed(qseq: str, tseq: str, seed_len: int) -> List[Tuple[int,int]]:
    """
    Return list of (q_pos, t_pos) positions where the seed (length seed_len)
    in qseq (starting at q_pos) is the perfect complement of tseq at t_pos.
    qseq is oriented 5'->3' and we search for hybridization q[i] pairing with t[j].
    """
    hits = []
    qlen = len(qseq); tlen = len(tseq)
    if seed_len <= 0:
        return hits
    # build dict of target k-mers for fast lookup
    t_kmers = {}
    for j in range(0, tlen - seed_len + 1):
        k = tseq[j:j+seed_len].upper()
        t_kmers.setdefault(k, []).append(j)
    # iterate query k-mers, but compare to complement
    for i in range(0, qlen - seed_len + 1):
        qk = qseq[i:i+seed_len].upper()
        # compute complement of qk (in target orientation)
        comp = ''.join(complement(x) for x in qk)
        # direct dictionary lookup
        if comp in t_kmers:
            for j in t_kmers[comp]:
                hits.append((i, j))
    return hits

def extend_seed_greedy(qseq: str, tseq: str, qpos: int, tpos: int, max_extend: Optional[int]=None) -> Tuple[int,int,float]:
    """
    Greedily extend left and right from the seed to form a contiguous (no bulges) duplex.
    Returns (q_start, t_start, energy) for the best contiguous block containing the seed.
    This is intentionally simple & fast: Option B's heuristic.
    """
    qlen = len(qseq); tlen = len(tseq)
    # find maximal left extension
    q_left = qpos
    t_left = tpos
    energy = 0.0
    # extend left
    while q_left > 0 and t_left > 0:
        a = qseq[q_left-1]; b = tseq[t_left-1]
        e = pair_energy(a,b)
        # accept extension if energy improves (i.e., e < 0), or small penalty allowed
        # This greedy rule is simple and keeps blocks largely perfectly pairing
        if e <= 1.0:  # allow pairings and small mismatches
            energy += e
            q_left -= 1; t_left -= 1
            if max_extend and (qpos - q_left) > max_extend:
                break
        else:
            break
    # extend right from end of seed
    seed_len_guess = 1  # we don't know actual seed size here; assume we start at seed length 1
    q_right = qpos + seed_len_guess
    t_right = tpos + seed_len_guess
    # continue while in bounds
    while q_right < qlen and t_right < tlen:
        a = qseq[q_right]; b = tseq[t_right]
        e = pair_energy(a,b)
        if e <= 1.0:
            energy += e
            q_right += 1; t_right += 1
            if max_extend and (q_right - qpos) > max_extend:
                break
        else:
            break
    # the above logic aims to return a contiguous block that contains the seed position.
    return (q_left, t_left, energy)

def fast_scan_all(query_fa: str, target_fa: str,
                  seed_len: int = 8,
                  min_seed_hits: int = 1,
                  max_hits_per_query: Optional[int] = None,
                  approx_score_threshold: float = -4.0) -> List[Dict]:
    """
    For each query sequence and each target window, locate seed matches and compute
    a very fast approximate score. Return list of candidate dicts:
      {
        'query_id': ...,
        'target_id': ...,
        'q_pos': ..., 't_pos': ...,
        'approx_energy': ...,
        'q_len': ...,
        't_len': ...
      }
    """
    queries = read_fasta_to_dict(query_fa)
    targets = read_fasta_to_dict(target_fa)
    candidates = []
    for qid, qseq in queries.items():
        for tid, tseq in targets.items():
            seed_hits = find_perfect_complement_seed(qseq, tseq, seed_len)
            if len(seed_hits) < min_seed_hits:
                continue
            # score each seed by extending greedily (fast)
            local_hits = []
            for (qpos, tpos) in seed_hits:
                qstart, tstart, ext_energy = extend_seed_greedy(qseq, tseq, qpos, tpos, max_extend=200)
                # compute rough contiguous block length from qstart/tstart to the end of seed-derived extension
                # For simplicity use ext_energy as approx; penalize very short matches
                approx_energy = ext_energy
                if approx_energy <= approx_score_threshold:
                    local_hits.append((qpos, tpos, qstart, tstart, approx_energy))
            # sort and keep top N per query/target combo if requested
            local_hits.sort(key=lambda x: x[4])  # more negative is better
            if max_hits_per_query:
                local_hits = local_hits[:max_hits_per_query]
            for qpos, tpos, qstart, tstart, approx_energy in local_hits:
                candidates.append({
                    'query_id': qid,
                    'target_id': tid,
                    'q_pos': int(qpos),
                    't_pos': int(tpos),
                    'approx_energy': float(approx_energy),
                    'q_len': len(qseq),
                    't_len': len(tseq)
                })
    print(f"[fast_scan_all] found {len(candidates)} candidate seeds (seed_len={seed_len})")
    return candidates

# -------------------------
# Moderate DP rescoring (Option B)
# Smith-Waterman style local alignment on base-pairing with simple pair energies
# -------------------------
def sw_rna_pairing_energy(qseq: str, tseq: str,
                          gap_penalty: float = 2.0) -> Tuple[float, Tuple[int,int,int,int]]:
    """
    Compute a Smith-Waterman local alignment between qseq and tseq, where match score is
    the negative base pairing energy (since pair_energy returns negative for favorable).
    Returns (best_energy, (q_start, q_end, t_start, t_end)) where positions are 0-based,
    q_end and t_end are exclusive.
    Energy is the sum of pair_energy for aligned columns plus gap penalties for gaps.
    The convention here: lower (more negative) energy = better.
    """
    m = len(qseq); n = len(tseq)
    # H matrix for best score ending at i,j; we use energies (negative good) and convert to scores by (-energy)
    # To keep consistent, we store energies directly (lower better). For dynamic programming we will accumulate.
    # For Smith-Waterman we need to allow restarting (i.e., zero baseline); here baseline is 0 energy and we look for negative minima.
    # Implement typical scoring where we track best (most negative) score.
    # We'll implement as H[i+1][j+1] = min(  H[i][j] + pair_energy(q[i],t[j]) , H[i][j+1] + gap_penalty, H[i+1][j] + gap_penalty, 0 )
    # But because pair_energy gives negative for favorable, 'min' gives best.
    INF = 1e9
    H = [[0.0]*(n+1) for _ in range(m+1)]
    best_energy = 0.0  # 0 is baseline (no interaction). We want negatives.
    best_coords = (0,0,0,0)
    for i in range(m):
        for j in range(n):
            e_pair = pair_energy(qseq[i], tseq[j])
            diag = H[i][j] + e_pair
            up = H[i][j+1] + gap_penalty
            left = H[i+1][j] + gap_penalty
            val = min(diag, up, left, 0.0)
            H[i+1][j+1] = val
            if val < best_energy:
                # update best; coordinates: extend backwards to find where this segment started (traceback naive)
                best_energy = val
                # find start by tracing back until cell is 0
                bi, bj = i+1, j+1
                # naive traceback
                while bi > 0 and bj > 0 and H[bi][bj] != 0.0:
                    # choose predecessor with smallest value
                    preds = [(H[bi-1][bj-1], bi-1, bj-1),
                             (H[bi-1][bj], bi-1, bj),
                             (H[bi][bj-1], bi, bj-1)]
                    preds.sort(key=lambda x: x[0])
                    valpred, nbi, nbj = preds[0]
                    if valpred >= H[bi][bj]:
                        # can't move; break to avoid infinite loop
                        break
                    bi, bj = nbi, nbj
                q_start = bi
                t_start = bj
                q_end = i+1
                t_end = j+1
                best_coords = (q_start, q_end, t_start, t_end)
    return best_energy, best_coords

# -------------------------
# Option C: higher-accuracy DP (heavier)
# A more permissive alignment with separate gap open and gap extend penalties,
# and slightly different parameterization to emulate a more complex model.
# This is purposely slower and more expensive.
# -------------------------
def affine_gap_rna_energy(qseq: str, tseq: str,
                          gap_open: float = 5.0,
                          gap_extend: float = 0.5) -> Tuple[float, Tuple[int,int,int,int]]:
    """
    Smith-Waterman with affine gap penalties (gap_open + k*gap_extend).
    Returns best energy and coordinates as in sw_rna_pairing_energy.
    Lower energy better (more negative).
    Complexity O(m*n). Slower in Python due to more bookkeeping.
    """
    m = len(qseq); n = len(tseq)
    # matrices: M (match/mismatch), Ix (gap in x -> gap in qseq), Iy (gap in y -> gap in tseq)
    # We'll store energies (lower better). Baseline 0 (no interaction).
    # initialize with +inf
    INF = 1e9
    M = [[0.0]*(n+1) for _ in range(m+1)]
    Ix = [[INF]*(n+1) for _ in range(m+1)]
    Iy = [[INF]*(n+1) for _ in range(m+1)]
    best_energy = 0.0
    best_coords = (0,0,0,0)
    for i in range(1, m+1):
        for j in range(1, n+1):
            e_pair = pair_energy(qseq[i-1], tseq[j-1])
            # match matrix: from M, Ix, Iy diagonals
            fromM = M[i-1][j-1] + e_pair
            fromIx = Ix[i-1][j-1] + e_pair
            fromIy = Iy[i-1][j-1] + e_pair
            M[i][j] = min(fromM, fromIx, fromIy, 0.0)
            # Ix: gap in qseq (i.e., insertions in target / deletion in query)
            Ix[i][j] = min(M[i-1][j] + gap_open + gap_extend,
                           Ix[i-1][j] + gap_extend,
                           Iy[i-1][j] + gap_open + gap_extend,
                           0.0)
            # Iy: gap in tseq
            Iy[i][j] = min(M[i][j-1] + gap_open + gap_extend,
                           Iy[i][j-1] + gap_extend,
                           Ix[i][j-1] + gap_open + gap_extend,
                           0.0)
            cell_best = min(M[i][j], Ix[i][j], Iy[i][j])
            if cell_best < best_energy:
                # naive traceback to find start (walk backwards until 0)
                best_energy = cell_best
                bi, bj = i, j
                # simple stepping back while cell != 0
                while bi > 0 and bj > 0:
                    cur = min(M[bi][bj], Ix[bi][bj], Iy[bi][bj])
                    if cur == 0.0:
                        break
                    # pick predecessor
                    # check if M is minimal
                    if cur == M[bi][bj]:
                        # came from diagonal
                        bi -= 1; bj -= 1
                    elif cur == Ix[bi][bj]:
                        # came from up
                        bi -= 1
                    else:
                        bj -= 1
                q_start = bi
                t_start = bj
                q_end = i
                t_end = j
                best_coords = (q_start, q_end, t_start, t_end)
    return best_energy, best_coords

# -------------------------
# Helpers for rescoring and aggregation
# -------------------------
import time
from typing import List, Dict

def rescore_candidates(
    candidates: List[Dict],
    query_fa: str,
    target_fa: str,
    mode: str = "B",
    top_n_for_C: int = 20,
    log_callback=None
) -> List[Dict]:
    """
    Rescore candidate interactions with dynamic programming.
    
    Parameters
    ----------
    candidates : list of dicts
        Each dict must contain 'query_id' and 'target_id'.
    query_fa : str
        Path to FASTA file of query windows.
    target_fa : str
        Path to FASTA file of targets.
    mode : str
        'B' for moderate DP (default), 'C' for full thermodynamic DP.
    top_n_for_C : int
        Number of top-B hits to re-score with Mode C.
    log_callback : callable
        Function to log progress (e.g., Streamlit callback).
    
    Returns
    -------
    rescored : list of dicts
        Each dict contains the original candidate fields + rescoring fields.
    """
    from Bio import SeqIO

    def log(msg):
        if log_callback:
            log_callback(msg)
        else:
            print(msg)

    # Load sequences into memory
    queries = {rec.id: str(rec.seq) for rec in SeqIO.parse(query_fa, "fasta")}
    targets = {rec.id: str(rec.seq) for rec in SeqIO.parse(target_fa, "fasta")}

    rescored = []
    n_candidates = len(candidates)

    # If mode C, only re-score top N candidates (sorted by approx_energy)
    candidates_to_score = candidates
    if mode.upper() == "C":
        candidates_to_score = sorted(candidates, key=lambda x: x.get("approx_energy", 0.0))[:top_n_for_C]

    log(f"Starting rescoring {len(candidates_to_score)} candidates in Mode {mode}...")

    for i, cand in enumerate(candidates_to_score, 1):
        start_time = time.time()
        qid = cand["query_id"]
        tid = cand["target_id"]
        qseq = queries.get(qid, "")
        tseq = targets.get(tid, "")

        # Simple DP scoring placeholder (replace with your DP function)
        # Example: negative Hamming distance as mock energy
        score = -sum(1 for a, b in zip(qseq, tseq[:len(qseq)]) if a != b)

        # Optionally adjust scoring for mode
        if mode.upper() == "B":
            cand["rescore_energy_B"] = score
            cand["q_start_B"] = 1
            cand["q_end_B"] = len(qseq)
            cand["t_start_B"] = 1
            cand["t_end_B"] = len(qseq)
        elif mode.upper() == "C":
            cand["rescore_energy_C"] = score
            cand["q_start_C"] = 1
            cand["q_end_C"] = len(qseq)
            cand["t_start_C"] = 1
            cand["t_end_C"] = len(qseq)

        elapsed = time.time() - start_time
        log(f"Rescored candidate {i}/{len(candidates_to_score)}: {qid} -> {tid}, energy={score}, time={elapsed:.2f}s")

        rescored.append(cand)

    log("Rescoring completed.")
    return rescored

def write_candidates_tsv(candidates: List[Dict], out_tsv: str):
    if not candidates:
        open(out_tsv, 'w').close()
        print(f"[write_candidates_tsv] wrote empty file {out_tsv}")
        return
    keys = list(candidates[0].keys())
    with open(out_tsv, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys, delimiter='\t')
        w.writeheader()
        for c in candidates:
            # convert numpy types if any by str()
            row = {k: (v if not isinstance(v, (list, dict)) else str(v)) for k,v in c.items()}
            w.writerow(row)
    print(f"[write_candidates_tsv] wrote {len(candidates)} rows to {out_tsv}")

def write_rescored_csv(results: List[Dict], out_csv: str):
    if not results:
        open(out_csv, 'w').close()
        print(f"[write_rescored_csv] wrote empty file {out_csv}")
        return
    # flatten keys union
    keys = sorted({k for r in results for k in r.keys()})
    with open(out_csv, 'w', newline='') as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        for r in results:
            safe = {k: (v if not isinstance(v, (list, dict)) else str(v)) for k,v in r.items()}
            w.writerow(safe)
    print(f"[write_rescored_csv] wrote {len(results)} rows to {out_csv}")

# -------------------------
# CLI
# -------------------------
def build_arg_parser():
    p = argparse.ArgumentParser(description="Self-contained RNA-RNA sliding-window pipeline.")
    sub = p.add_subparsers(dest='cmd', required=True)

    # make-windows
    mw = sub.add_parser('make-windows', help='Generate sliding windows from a FASTA')
    mw.add_argument('--seq', required=True, help='input FASTA with targets')
    mw.add_argument('--window', type=int, required=True, help='window length')
    mw.add_argument('--step', type=int, required=True, help='step size')
    mw.add_argument('--out', required=True, help='output FASTA windows')

    # scan
    sc = sub.add_parser('scan', help='Fast seed-based scan to generate candidates')
    sc.add_argument('--query', required=True, help='query FASTA (one-or-more sequences)')
    sc.add_argument('--target', required=True, help='target FASTA (windows)')
    sc.add_argument('--out', required=True, help='output candidates TSV')
    sc.add_argument('--seed_len', type=int, default=8, help='seed length for perfect complement')
    sc.add_argument('--min_seed_hits', type=int, default=1, help='minimum seeds required per query-target')
    sc.add_argument('--approx_energy_thresh', type=float, default=-4.0, help='approx energy threshold for candidate inclusion')
    sc.add_argument('--max_hits_per_query', type=int, default=10, help='max hits per query-target pair')

    # rescore
    rs = sub.add_parser('rescore', help='Rescore candidate hits with Option B and optional C')
    rs.add_argument('--candidates', required=True, help='input candidates TSV (from scan)')
    rs.add_argument('--query', required=True, help='query FASTA')
    rs.add_argument('--target', required=True, help='target FASTA (windows)')
    rs.add_argument('--out', required=True, help='output rescored CSV')
    rs.add_argument('--mode', choices=['B','C'], default='B', help='Rescore mode: B (default) or C (run B then top-N with C)')
    rs.add_argument('--modeC_top', type=int, default=20, help='when mode=C: number of top-B hits to re-evaluate with C')

    return p

def read_candidates_tsv(path: str) -> List[Dict]:
    rows = []
    with open(path, 'r') as fh:
        reader = csv.DictReader(fh, delimiter='\t')
        for r in reader:
            # cast known fields to int/float
            row = dict(r)
            for k in ('q_pos','t_pos','q_len','t_len'):
                if k in row and row[k] != '':
                    row[k] = int(row[k])
            for k in ('approx_energy',):
                if k in row and row[k] != '':
                    row[k] = float(row[k])
            rows.append(row)
    return rows

def main():
    parser = build_arg_parser()
    args = parser.parse_args()
    if args.cmd == 'make-windows':
        generate_windows(args.seq, args.window, args.step, args.out)
        return

    if args.cmd == 'scan':
        candidates = fast_scan_all(
            query_fa=args.query,
            target_fa=args.target,
            seed_len=args.seed_len,
            min_seed_hits=args.min_seed_hits,
            max_hits_per_query=args.max_hits_per_query,
            approx_score_threshold=args.approx_energy_thresh
        )
        write_candidates_tsv(candidates, args.out)
        return

    if args.cmd == 'rescore':
        candidates = read_candidates_tsv(args.candidates)
        results = rescore_candidates(candidates, args.query, args.target, mode=args.mode, top_n_for_C=args.modeC_top)
        write_rescored_csv(results, args.out)
        return

if __name__ == '__main__':
    main()
