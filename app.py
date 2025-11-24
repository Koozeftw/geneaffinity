import streamlit as st
import os
import tempfile
from datetime import datetime
from pipeline_runner import run_pipeline  # pure Python version
import pandas as pd
import plotly.express as px

# --- Streamlit page config ---
st.set_page_config(page_title="FASTA Binding Pipeline", layout="wide")
st.title("FASTA Binding Pipeline 🔬")
st.write("Upload two FASTA sequences and run a fully self-contained Python pipeline to find potential binding regions.")

# --- File uploads ---
fileA = st.file_uploader("Upload Sequence A (FASTA)", type=["fa", "fasta"])
fileB = st.file_uploader("Upload Sequence B (FASTA)", type=["fa", "fasta"])

# --- Pipeline parameters ---
st.sidebar.header("Pipeline Parameters")
window = st.sidebar.number_input("Window Size", min_value=1, value=35)
step = st.sidebar.number_input("Step Size", min_value=1, value=5)
flank = st.sidebar.number_input("Flank Size", min_value=0, value=100)
energy_cutoff_fast = st.sidebar.number_input(
    "Fast Scan: Minimum Complementarity Score",
    value=-6.0,
    format="%.2f"
)
top_k = st.sidebar.number_input("Top Hits per Window to Keep", min_value=1, value=100)

# Mode B or C
mode = st.sidebar.radio(
    "Thermodynamic Mode",
    ["B (fast, simplified thermodynamics)", "C (full thermodynamics after B)"]
)
topN_C = None
if mode.startswith("C"):
    topN_C = st.sidebar.number_input(
        "Re-run full thermodynamics for top N hits",
        min_value=1,
        value=20
    )

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

            st.info("Running pipeline… this may take a moment.")
            log_window = st.empty()  # live logs

            try:
                # Run pipeline
                results_df = run_pipeline(
                    geneA_path=pathA,
                    geneB_path=pathB,
                    window_sizes=[window],
                    step=step,
                    flank=flank,
                    energy_cutoff_fast=energy_cutoff_fast,
                    top_k_per_window=top_k,
                    thermo_mode="B" if mode.startswith("B") else "C",
                    topN_modeC=topN_C,
                    log_callback=lambda msg: log_window.text(msg)
                )

                st.success("Pipeline finished! ✔️")

                # Display dataframe
                st.dataframe(results_df)

                # --- CSV download ---
                csv_data = results_df.to_csv(index=False)
                st.download_button(
                    label="⬇️ Download Results as CSV",
                    data=csv_data,
                    file_name=f"results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.csv",
                    mime="text/csv"
                )

                # --- Excel download ---
                try:
                    import io
                    excel_buffer = io.BytesIO()
                    results_df.to_excel(excel_buffer, index=False, engine="openpyxl")
                    st.download_button(
                        label="⬇️ Download Results as Excel",
                        data=excel_buffer.getvalue(),
                        file_name=f"results_{datetime.now().strftime('%Y%m%d_%H%M%S')}.xlsx",
                        mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet"
                    )
                except Exception as e:
                    st.warning(f"Could not generate Excel file: {e}")

                # --- Interactive Plotly visualization ---
                st.header("Interactive Binding Visualization")

                # Optional: limit to top 3 hits per query for clarity
                top_df = results_df.groupby('query_id').apply(lambda x: x.nsmallest(3, 'rescore_energy_B')).reset_index(drop=True)

                # Hover labels
                top_df['label'] = top_df.apply(
                    lambda row: f"Query: {row['query_id']}<br>Target: {row['target_id']}<br>Energy: {row['rescore_energy_B']}<br>Q: {row['q_start_B']}-{row['q_end_B']} T: {row['t_start_B']}-{row['t_end_B']}",
                    axis=1
                )

                # Plotly scatter
                fig = px.scatter(
                    top_df,
                    x='t_start_B',
                    y='query_id',
                    size=[row['t_end_B'] - row['t_start_B'] for _, row in top_df.iterrows()],
                    color='rescore_energy_B',
                    color_continuous_scale='RdYlBu_r',
                    hover_name='label',
                    labels={'query_id': 'Query Window', 't_start_B': 'Target Start'}
                )
                fig.update_traces(marker=dict(line=dict(width=1, color='black')))
                fig.update_layout(
                    xaxis_title='Target Sequence B Position',
                    yaxis_title='Query Windows (Sequence A)',
                    coloraxis_colorbar=dict(title='Mode B Energy')
                )
                st.plotly_chart(fig, use_container_width=True)

                # Dropdown for selecting a window
                selected_window = st.selectbox("Select a query window to see details:", top_df['query_id'].unique())
                window_data = top_df[top_df['query_id'] == selected_window]
                st.write(f"Details for {selected_window}:", window_data)

            except Exception as e:
                st.error(f"Pipeline failed: {e}")
