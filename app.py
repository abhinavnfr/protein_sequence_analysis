import os
import streamlit as st
import pandas as pd
from io import BytesIO
import time
import ingest_sequence as ingest
import generate_views as gv
import databricks_handler as dbh
import base64
import re


def main():
    st.set_page_config(
                        page_title="Protezard",
                        page_icon="ui_elements/logo.png"
                    )
    
    with open("ui_elements/dna_loop_animation.mp4", "rb") as video_file:
        encoded_string = base64.b64encode(video_file.read()).decode('utf-8')
    
    video_html = f"""
                    <style>
                    #myVideo {{
                        position: fixed;
                        right: 0;
                        bottom: 0;
                        min-width: 100%; 
                        min-height: 100%;
                    }}
                    .content {{
                        position: fixed;
                        bottom: 0;
                        background: rgba(0, 0, 0, 0.5);
                        color: #f1f1f1;
                        width: 100%;
                        padding: 20px;
                    }}
                    </style>
                    <video autoplay loop muted playsinline id="myVideo">
                        <source type="video/mp4" src="data:video/mp4;base64,{encoded_string}">
                    </video>
                """
    
    st.markdown(video_html, unsafe_allow_html=True)

    st.markdown("""
                    <style>
                    /* Make all labels and texts white */
                    /* Target labels with class 'css-1avcm0e' or similar */
                    label, p, .css-1cpxqw2, stFileUploaderFileName { 
                        color: #FFFFFF !important;
                    }
                    /* Specific for Streamlit labels and other texts */
                    .stText, .stMarkdown, .stButton, .stSuccess, .stWarning, .stError {
                        color: #FFFFFF !important;
                    }
                    </style>
                """,
                unsafe_allow_html=True
                )
    st.markdown("""
                    <style>
                    /* Change the spinner circle color to white */
                    .stSpinner > div > div {
                        border-top-color: #FFFFFF !important;
                    }
                    /* Change the spinner text color (including elapsed time) to white */
                    .stSpinner > div > span {
                        color: #FFFFFF !important;
                    }
                    </style>
                """,
                unsafe_allow_html=True
                )

    st.markdown("<h1 style='color: white;'>Protezard</h1>", unsafe_allow_html=True)
    st.markdown("<h2 style='color: white;'>A one-stop shop app for all your protein sequence analysis needs</h2><br>", unsafe_allow_html=True)

    # Step 1: Process input accession numbers file
    st.markdown("<p style='font-size: 24px; color: white;'>To get started, choose a text file containing accession numbers or a fasta file containing sequences</p>", unsafe_allow_html=True)
    input_file = st.file_uploader(label="Upload file", type=["txt", "fasta"])
    if input_file is not None:
        filename = input_file.name
        content = input_file.read().decode("utf-8")
        if filename.endswith(".txt"):
            accessions = [line.strip() for line in content.splitlines() if line.strip()]
            new_accessions = ingest.filter_new_sequences(accessions)
            new_accesions_count = len(new_accessions)
            ingest.add_new_accession_uc_table(new_accessions)
            ingest.add_fasta_uc_table()
        elif filename.endswith(".fasta"):
            fasta_records = re.split(r'(^>.*$)', content, flags=re.MULTILINE)
            seq_list = []
            accessions = []
            for i in range(1, len(fasta_records), 2):
                header = fasta_records[i].strip()
                sequence = fasta_records[i+1] if (i+1) < len(fasta_records) else ""
                fasta_seq = header + sequence
                accession_id = header[1:].split()[0]  # accession = first word after '>'
                seq_list.append((accession_id, fasta_seq.strip()))
                accessions.append(accession_id)
            st.write(accessions[0])
            ingest.add_fasta_batch_uc_table(seq_list)
    
        if st.button(label="Perform BLAST, EffectorP, PFAM Domain Search, and Molecular Weight Calculation end-to-end", type="primary"):
            # get sequences to BLAST
            seq_to_blast = ingest.get_seq_to_blast()
            seq_to_blast_count = len(seq_to_blast)

            for seq in seq_to_blast:
                try:
                    blasted_sequence = ingest.blast_sequence(seq[0], seq[1])
                    ingest.add_blast_uc_table(seq[0], blasted_sequence)
                except Exception as e:
                    st.error(f"Error processing {acc}: {e}")
                time.sleep(0.1)
            df_blast = gv.generate_view("workspace.curated.blast_sequence", accessions)
            results_blast = BytesIO()
            with pd.ExcelWriter(results_blast, engine='xlsxwriter') as writer:
                df_blast.to_excel(writer, index=False, sheet_name='results_blast')
            results_blast.seek(0)
            
            time.sleep(1)
            ingest.predict_effectorp()
            df_effectorp = gv.generate_view("workspace.curated.effectorp", accessions)
            results_effectorp = BytesIO()
            with pd.ExcelWriter(results_effectorp, engine='xlsxwriter') as writer:
                df_effectorp.to_excel(writer, index=False, sheet_name='results_effectorp')
            results_effectorp.seek(0)
            
            time.sleep(1)
            ingest.pfam_domain_search()
            df_pfam = gv.generate_view("workspace.curated.pfam", accessions)
            results_pfam = BytesIO()
            with pd.ExcelWriter(results_pfam, engine='xlsxwriter') as writer:
                df_pfam.to_excel(writer, index=False, sheet_name='results_pfam')
            results_pfam.seek(0)
            
            time.sleep(1)
            ingest.calculate_molecular_weight_kda()
            df_mw = gv.generate_view("workspace.curated.molecularweight", accessions)
            results_mw = BytesIO()
            with pd.ExcelWriter(results_mw, engine='xlsxwriter') as writer:
                df_mw.to_excel(writer, index=False, sheet_name='results_molecularweight')
            results_mw.seek(0)
            
            st.download_button(label="Download BLAST results", 
                            type="primary", 
                            icon=":material/download:",
                            data=results_blast,
                            file_name="results_blast.xlsx",
                            on_click="ignore",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
            st.download_button(label="Download EffectorP results", 
                            type="primary", 
                            icon=":material/download:",
                            data=results_effectorp,
                            file_name="results_effectorp.xlsx",
                            on_click="ignore",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
            st.download_button(label="Download InterProScan PFAM Domain Search results", 
                            type="primary", 
                            icon=":material/download:",
                            data=results_pfam,
                            file_name="results_pfam.xlsx",
                            on_click="ignore",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
            st.download_button(label="Download Molecular Weight results", 
                            type="primary", 
                            icon=":material/download:",
                            data=results_mw,
                            file_name="results_molecularweight.xlsx",
                            on_click="ignore",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

        elif st.button(label="Perform only BLAST", type="primary"):
            # get sequences to BLAST
            seq_to_blast = ingest.get_seq_to_blast()
            seq_to_blast_count = len(seq_to_blast)

            for seq in seq_to_blast:
                try:
                    blasted_sequence = ingest.blast_sequence(seq[0], seq[1])
                    ingest.add_blast_uc_table(seq[0], blasted_sequence)
                except Exception as e:
                    st.error(f"Error processing {acc}: {e}")
                time.sleep(0.1)
            df_blast = gv.generate_view("workspace.curated.blast_sequence", accessions)
            results_blast = BytesIO()
            with pd.ExcelWriter(results_blast, engine='xlsxwriter') as writer:
                df_blast.to_excel(writer, index=False, sheet_name='results_blast')
            results_blast.seek(0)
            st.download_button(label="Download BLAST results", 
                            type="primary", 
                            icon=":material/download:",
                            data=results_blast,
                            file_name="results_blast.xlsx",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
            
        elif st.button(label="Perform only EffectorP", type="primary"):
            time.sleep(1)
            ingest.predict_effectorp()
            df_effectorp = gv.generate_view("workspace.curated.effectorp", accessions)
            results_effectorp = BytesIO()
            with pd.ExcelWriter(results_effectorp, engine='xlsxwriter') as writer:
                df_effectorp.to_excel(writer, index=False, sheet_name='results_effectorp')
            results_effectorp.seek(0)
            st.download_button(label="Download EffectorP results", 
                            type="primary", 
                            icon=":material/download:",
                            data=results_effectorp,
                            file_name="results_effectorp.xlsx",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        
        elif st.button(label="Perform only InterproScan PFAM Domain Search", type="primary"):
            time.sleep(1)
            ingest.pfam_domain_search()
            df_pfam = gv.generate_view("workspace.curated.pfam", accessions)
            results_pfam = BytesIO()
            with pd.ExcelWriter(results_pfam, engine='xlsxwriter') as writer:
                df_pfam.to_excel(writer, index=False, sheet_name='results_pfam')
            results_pfam.seek(0)
            st.download_button(label="Download InterProScan PFAM Domain Search results", 
                            type="primary", 
                            icon=":material/download:",
                            data=results_pfam,
                            file_name="results_pfam.xlsx",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
        
        elif st.button(label="Perform only Molecular Weight Calculation", type="primary"):
            time.sleep(1)
            ingest.calculate_molecular_weight_kda()
            df_mw = gv.generate_view("workspace.curated.molecularweight", accessions)
            results_mw = BytesIO()
            with pd.ExcelWriter(results_mw, engine='xlsxwriter') as writer:
                df_mw.to_excel(writer, index=False, sheet_name='results_molecularweight')
            results_mw.seek(0)
            st.download_button(label="Download Molecular Weight results", 
                            type="primary", 
                            icon=":material/download:",
                            data=results_mw,
                            file_name="results_molecularweight.xlsx",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")


if __name__ == "__main__":
    main()
