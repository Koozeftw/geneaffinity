import streamlit as st
import os
import tempfile
from datetime import datetime
from pipeline_runner import run_pipeline  # updated runner with local binaries

# --- Streamlit page config ---
st.set_page_config(page_title="FASTA Binding Pipeline", layout="wide")
st.title("FASTA Binding Pipeline 🔬")
st.write("Upload your sequences and run the pipeline to find potential binding regions.")

# --- File uploads ---
fileA = st.file_uploader("Upload Sequence A (FASTA)", type=["fa", "fasta"])
fileB = st.file_uploader("Upload Sequence B (FASTA)", type=["fa", "fasta"])

# --- Pipeline parameters ---
st.sidebar.header("Pipeline Parameters")
window = st.sidebar.number_input("Window Size", min_value=1, value=35)
step = st.sidebar.number_input("Step Size", min_value=1, value=5)
flank = st.sidebar.number_input("Flank Size", min_value=0, value=100)
energy_cutoff_fast = st.sidebar.number_input("Energy Cutoff (Fast)", value=-6.0, format="%.2f")
top_k = st.sidebar.number_input("Top Hits per Query", min_value=1, value=100)
threads = st.sidebar.number_input("Threads", min_value=1, value=1)

# --- Run button ---
if st.button("Run Pipeline"):
    if not fileA or not fileB:
        st.error("Please upload both Sequence A and Sequence B.")
    else:
        with tempfile.TemporaryDirectory() as tmpdir:
            # Save uploaded files
            pathA = os.path.join(tmpdir, fileA.name)
            pathB = os.path.join(tmpdir, fileB.name)
            with open(pathA, "wb") as f:
                f.write(fileA.getbuffer())
            with open(pathB, "wb") as f:
                f.write(fileB.getbuffer())

            st.info("Running pipeline...")
            log_window = st.empty()  # live log placeholder

            # Define local binaries
            risearch_bin = os.path.join(os.path.dirname(__file__), "risearch2")
            intarna_bin = os.path.join(os.path.dirname(__file__), "IntaRNA")

            # Run pipeline
            try:
                results_df = run_pipeline(
                    geneA_path=pathA,
                    geneB_path=pathB,
                    window_sizes=[window],
                    step=step,
                    flank=flank,
                    energy_cutoff_fast=energy_cutoff_fast,
                    top_k_per_window=top_k,
                    threads=threads,
                    risearch_bin=risearch_bin,
                    intarna_bin=intarna_bin,
                    log_callback=lambda msg: log_window.text(msg)  # live logging
                )

                st.success("Pipeline finished!")

                # Download results as CSV
                csv_data = results_df.to_csv(index=False)
                st.download_button(
                    label="⬇️ Download Results CSV",
                    data=csv_data,
                    file_name=f"aggregated_results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                    mime="text/csv"
                )

            except Exception as e:
                st.error(f"Pipeline failed: {e}")
