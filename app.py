import streamlit as st
import os
import tempfile
from datetime import datetime
from pipeline_runner import run_pipeline
import pandas as pd
import plotly.express as px
from Bio import SeqIO

st.set_page_config(page_title="FASTA Binding Pipeline", layout="wide")
st.title("FASTA Binding Pipeline 🔬")
st.write("Upload two FASTA sequences and run a fully self-contained Python pipeline.")

# --- File uploads ---
fileA = st.file_uploader("Upload Sequence A (FASTA)", type=["fa", "fasta"])
fileB = st.file_uploader("Upload Sequence B (FASTA)", type=["fa", "fasta"])

# --- Parameters ---
st.sidebar.header("Pipeline Parameters")
window = st.sidebar.number_input("Window Size", min_value=1, value=35)
step = st.sidebar.number_input("Step Size", min_value=1, value=5)
flank = st.sidebar.number_input("Flank Size", min_value=0, value=100)
top_k = st.sidebar.number_input("Top Hits per Window", min_value=1, value=100)
mode = st.sidebar.radio("Thermodynamic Mode", ["B", "C"])
topN_C = None
if mode=="C":
    topN_C = st.sidebar.number_input("Top N for Mode C", min_value=1, value=20)

# --- Run button ---
if st.button("Run Pipeline"):
    if not fileA or not fileB:
        st.error("Please upload both sequences.")
    else:
        with tempfile.TemporaryDirectory() as tmpdir:
            pathA = os.path.join(tmpdir, fileA.name)
            pathB = os.path.join(tmpdir, fileB.name)
            with open(pathA, "wb") as f: f.write(fileA.getbuffer())
            with open(pathB, "wb") as f: f.write(fileB.getbuffer())

            log_window = st.empty()
            def log(msg): log_window.text(msg)

            try:
                results_df = run_pipeline(
                    geneA_path=pathA,
                    geneB_path=pathB,
                    window_sizes=[window],
                    step=step,
                    flank=flank,
                    top_k_per_window=top_k,
                    thermo_mode=mode,
                    topN_modeC=topN_C,
                    log_callback=log
                )

                st.success("Pipeline finished!")

                # --- Add actual sequences dynamically ---
                seqA_dict = {rec.id: str(rec.seq) for rec in SeqIO.parse(pathA, "fasta")}
                seqB_dict = {rec.id: str(rec.seq) for rec in SeqIO.parse(pathB, "fasta")}

                # Detect coordinate columns
                q_start = next(c for c in results_df.columns if c.startswith('q_start'))
                q_end   = next(c for c in results_df.columns if c.startswith('q_end'))
                t_start = next(c for c in results_df.columns if c.startswith('t_start'))
                t_end   = next(c for c in results_df.columns if c.startswith('t_end'))

                for c in [q_start,q_end,t_start,t_end]:
                    results_df[c] = results_df[c].astype(int)

                results_df['query_seq'] = results_df.apply(
                    lambda row: seqA_dict.get(row['query_id'].split('|')[0],'')[row[q_start]-1:row[q_end]],
                    axis=1
                )
                results_df['target_seq'] = results_df.apply(
                    lambda row: seqB_dict.get(row['target_id'], '')[row[t_start]-1:row[t_end]],
                    axis=1
                )

                st.dataframe(results_df)

                # --- Downloads ---
                csv_data = results_df.to_csv(index=False)
                st.download_button("Download CSV", csv_data,
                                   file_name=f"results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                                   mime="text/csv")

                try:
                    import io
                    excel_buffer = io.BytesIO()
                    results_df.to_excel(excel_buffer, index=False, engine="openpyxl")
                    st.download_button("Download Excel", excel_buffer.getvalue(),
                                       file_name=f"results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx",
                                       mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
                except Exception as e:
                    st.warning(f"Excel export failed: {e}")

                # --- Interactive Plot ---
                top_df = results_df.groupby('query_id').apply(lambda x: x.nlargest(3,'score')).reset_index(drop=True)
                top_df['label'] = top_df.apply(lambda r: f"Q:{r['query_id']} T:{r['target_id']} Score:{r['score']}<br>Query_seq:{r['query_seq']}<br>Target_seq:{r['target_seq']}", axis=1)

                fig = px.scatter(
                    top_df,
                    x='t_start',
                    y='query_id',
                    size=top_df[t_end]-top_df[t_start],
                    color='score',
                    color_continuous_scale='RdYlBu_r',
                    hover_name='label'
                )
                fig.update_traces(marker=dict(line=dict(width=1,color='black')))
                fig.update_layout(
                    xaxis_title='Target Sequence B Position',
                    yaxis_title='Query Windows',
                    coloraxis_colorbar=dict(title='Score')
                )
                st.plotly_chart(fig, use_container_width=True)

            except Exception as e:
                st.error(f"Pipeline failed: {e}")
