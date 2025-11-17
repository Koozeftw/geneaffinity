import streamlit as st
import subprocess
import os
from datetime import datetime
import tempfile

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
energy_cutoff = st.sidebar.number_input("Energy Cutoff", value=-10.0, format="%.2f")
threads = st.sidebar.number_input("Threads", min_value=1, value=1)

# --- Run button ---
if st.button("Run Pipeline"):
    if not fileA or not fileB:
        st.error("Please upload both Sequence A and Sequence B.")
    else:
        with tempfile.TemporaryDirectory() as tmpdir:
            pathA = os.path.join(tmpdir, fileA.name)
            pathB = os.path.join(tmpdir, fileB.name)
            # Save uploaded files
            with open(pathA, "wb") as f:
                f.write(fileA.getbuffer())
            with open(pathB, "wb") as f:
                f.write(fileB.getbuffer())

            st.info("Running pipeline...")
            log_window = st.empty()  # Placeholder for live logs

            # Local binaries
            risearch_bin = os.path.join(os.path.dirname(__file__), "risearch2")
            intarna_bin = os.path.join(os.path.dirname(__file__), "IntaRNA")

            cmd = [
                "python3",
                "pipeline.py",
                "fast-search",
                "--query", pathA,
                "--target", pathB,
                "--out", os.path.join(tmpdir, "results_risearch.tsv"),
                "--energy", str(energy_cutoff),
                "--max-hits", "1000",
            ]

            # Export environment variables for local binaries
            env = os.environ.copy()
            env["RISEARCH_BIN"] = risearch_bin
            env["INTARNA_BIN"] = intarna_bin

            process = subprocess.Popen(
                cmd,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
                text=True,
                bufsize=1,
                env=env
            )

            full_log = ""
            for line in process.stdout:
                full_log += line
                log_window.code(full_log)

            process.wait()

            results_path = os.path.join(tmpdir, "results_risearch.tsv")
            if os.path.exists(results_path):
                st.success("Pipeline finished!")
                with open(results_path, "r") as f:
                    st.download_button(
                        label="⬇️ Download Results",
                        data=f.read(),
                        file_name=f"results_{}.tsv".format(datetime.now().strftime("%Y%m%d_%H%M%S")),
                        mime="text/plain",
                    )
            else:
                st.error("Pipeline finished but results were not found. Check logs above.")
