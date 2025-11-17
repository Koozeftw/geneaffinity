import os
import subprocess
import pandas as pd
from pathlib import Path
from pipeline import generate_windows, parse_risearch_tsv, extract_contexts, call_intarna, aggregate_and_rank
import shutil
import urllib.request
import stat

# ============================
# Helper: download binaries
# ============================

def download_binary(url, dest_path):
    if dest_path.exists():
        return
    print(f"Downloading {url} → {dest_path}")
    urllib.request.urlretrieve(url, dest_path)
    # Make it executable
    st = os.stat(dest_path)
    os.chmod(dest_path, st.st_mode | stat.S_IEXEC)

# ============================
# Check and setup binaries
# ============================

def setup_risearch2(tmp_dir, log_callback=print):
    risearch_path = tmp_dir / "risearch2"
    if shutil.which("risearch2"):
        log_callback("✅ Found system risearch2")
        return shutil.which("risearch2")
    
    # Otherwise, attempt download
    log_callback("⚠️ risearch2 not found. Attempting to download precompiled binary...")
    # Example URL (replace with actual working precompiled binary URL)
    url = "https://example.com/risearch2_linux_x86_64"
    try:
        download_binary(url, risearch_path)
        return str(risearch_path)
    except Exception as e:
        log_callback(f"❌ Failed to get risearch2: {e}")
        return None

def setup_intarna(tmp_dir, log_callback=print):
    intarna_path = tmp_dir / "IntaRNA"
    if shutil.which("IntaRNA"):
        log_callback("✅ Found system IntaRNA")
        return shutil.which("IntaRNA")
    
    # Otherwise, attempt download
    log_callback("⚠️ IntaRNA not found. Attempting to download precompiled binary...")
    # Example URL (replace with actual working precompiled binary URL)
    url = "https://example.com/IntaRNA_linux_x86_64"
    try:
        download_binary(url, intarna_path)
        return str(intarna_path)
    except Exception as e:
        log_callback(f"❌ Failed to get IntaRNA: {e}")
        return None

# ============================
# Main pipeline runner
# ============================

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

    # Setup binaries
    risearch2_bin = setup_risearch2(tmp_dir, log_callback)
    intarna_bin = setup_intarna(tmp_dir, log_callback)

    results_all = []

    for w in window_sizes:
        log_callback(f"=== Running window size {w} ===")
        # 1️⃣ Generate windows
        win_fa = tmp_dir / f"windows_w{w}.fa"
        generate_windows(geneA_path, w, step, win_fa)

        # 2️⃣ Fast search
        fast_out = tmp_dir / f"risearch_w{w}.tsv"
        if risearch2_bin:
            try:
                cmd = [
                    risearch2_bin,
                    "-q", str(win_fa),
                    "-t", str(geneB_path),
                    "--energy-cutoff", str(energy_cutoff_fast),
                    "--max-hits", str(100000),
                    "-o", str(fast_out)
                ]
                log_callback(f"Running risearch2: {cmd}")
                subprocess.run(cmd, check=True)
            except Exception as e:
                log_callback(f"❌ risearch2 failed: {e}")
                continue
        else:
            log_callback("⚠️ Skipping risearch2 (binary not available). Creating empty hits file.")
            pd.DataFrame(columns=['query_id','target_id','q_start','q_end','t_start','t_end','energy']).to_csv(fast_out, sep="\t", index=False)

        # 3️⃣ Parse results
        hits_df = parse_risearch_tsv(fast_out)
        hits_df = hits_df.sort_values('energy').groupby('query_id').head(top_k_per_window).reset_index(drop=True)

        # 4️⃣ Extract context
        contexts_dir = tmp_dir / f"contexts_w{w}"
        contexts_dir.mkdir(exist_ok=True)
        contexts_df = extract_contexts(geneB_path, hits_df, flank, contexts_dir)

        # 5️⃣ Rescore with IntaRNA
        intarna_results = []
        if intarna_bin:
            for idx, r in contexts_df.iterrows():
                qid = hits_df.loc[int(r['hit_index']), 'query_id']
                qtmp = contexts_dir / f"query_{qid}.fa"
                # write single window query
                from Bio import SeqIO
                for rec in SeqIO.parse(win_fa, "fasta"):
                    if rec.id == qid:
                        with open(qtmp, "w") as fh:
                            fh.write(f">{qid}\n{str(rec.seq)}\n")
                        break
                out_pref = contexts_dir / f"intarna_{idx}"
                try:
                    csvp = call_intarna(qtmp, r['context_fasta'], out_pref, threads=threads)
                    intarna_results.append((int(r['hit_index']), csvp))
                except Exception as e:
                    log_callback(f"❌ IntaRNA failed for {qtmp}: {e}")
        else:
            log_callback("⚠️ Skipping IntaRNA (binary not available). Using empty results.")
            intarna_results = []

        # 6️⃣ Aggregate
        agg_out = tmp_dir / f"aggregated_w{w}.csv"
        agg_df = aggregate_and_rank(intarna_results, hits_df, agg_out)
        if agg_df is not None:
            results_all.append(agg_df)

    if results_all:
        return pd.concat(results_all, ignore_index=True)
    else:
        return pd.DataFrame()
