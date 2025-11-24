import streamlit as st
import tempfile, os
from datetime import datetime
from Bio import SeqIO
import pandas as pd
import plotly.express as px
from pipeline_runner import run_pipeline

st.set_page_config(page_title="FASTA Binding Pipeline", layout="wide")
st.title("FASTA Binding Pipeline 🔬")
st.write("Upload two FASTA sequences and run the pipeline.")

fileA = st.file_uploader("Sequence A (FASTA)", type=["fa","fasta"])
fileB = st.file_uploader("Sequence B (FASTA)", type=["fa","fasta"])

st.sidebar.header("Parameters")
window = st.sidebar.number_input("Window size", min_value=1, value=35)
step = st.sidebar.number_input("Step size", min_value=1, value=1)
top_k = st.sidebar.number_input("Top hits per window", min_value=1, value=20)

if st.button("Run Pipeline") and fileA and fileB:
    with tempfile.TemporaryDirectory() as tmpdir:
        pathA = os.path.join(tmpdir, fileA.name)
        pathB = os.path.join(tmpdir, fileB.name)
        with open(pathA, "wb") as f: f.write(fileA.getbuffer())
        with open(pathB, "wb") as f: f.write(fileB.getbuffer())

        log_window = st.empty()
        def logger(msg): log_window.text(msg)

        results_df = run_pipeline(
            geneA_path=pathA,
            geneB_path=pathB,
            window_sizes=[window],
            step=step,
            top_k_per_window=top_k,
            log_callback=logger
        )

        st.success("Pipeline finished!")
        st.dataframe(results_df)

        # CSV download
        csv_data = results_df.to_csv(index=False)
        st.download_button("Download CSV", csv_data, f"results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv", "text/csv")

        # Plotly visualization
        st.header("Interactive visualization")
        # Optional: show top 3 hits per query for clarity
        top_df = results_df.groupby('query_id').apply(lambda x: x.nlargest(3, 'rescore_energy_B')).reset_index(drop=True)
        top_df['label'] = top_df.apply(lambda r: f"Query: {r['query_id']}<br>Target: {r['target_id']}<br>Score: {r['rescore_energy_B']}<br>Query seq: {r['query_seq']}<br>Target seq: {r['target_seq']}", axis=1)

        fig = px.scatter(
            top_df,
            x='t_start',
            y='query_id',
            size=[r['t_end']-r['t_start'] for _,r in top_df.iterrows()],
            color='rescore_energy_B',
            color_continuous_scale='RdYlBu_r',
            hover_name='label',
            labels={'query_id':'Query window','t_start':'Target start'}
        )
        fig.update_traces(marker=dict(line=dict(width=1, color='black')))
        st.plotly_chart(fig, use_container_width=True)
