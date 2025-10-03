import os
import streamlit as st
import pandas as pd
from io import BytesIO
import time
import ingest_sequence as ingest
import generate_views as gv
import databricks_handler as dbh
import base64


def main():
    st.set_page_config(
                        page_title="Protezard",
                        page_icon="🧬"
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
    st.markdown("""
                    <style>
                    /* Style for secondary download buttons */
                    button[kind="secondary"] {
                        background-color: #FFFFFF !important;  /* White background */
                        color: #000000 !important;  /* Dark text */
                        border: 1px solid #CCCCCC !important;  /* Optional border for visibility */
                    }
                    button[kind="secondary"]:hover {
                        background-color: #EEEEEE !important;  /* Slightly gray on hover */
                        color: #000000 !important;  /* Keep text dark on hover */
                    }
                    </style>
                """,
                unsafe_allow_html=True
                )

    st.markdown("<h1 style='color: white;'>Protezard</h1>", unsafe_allow_html=True)
    st.markdown("<h2 style='color: white;'>A one-stop shop app for all your protein sequence analysis needs</h2><br>", unsafe_allow_html=True)

    # Step 1: Process input accession numbers file
    st.markdown("<p style='font-size: 24px; color: white;'>To get started, choose a text file containing accession numbers</p>", unsafe_allow_html=True)
    input_file = st.file_uploader(label="Upload file", type=["txt"])
    if input_file is not None:
        accessions = [line.strip() for line in input_file.read().decode("utf-8").splitlines()]
        new_accessions = ingest.filter_new_sequences(accessions)
        new_accesions_count = len(new_accessions)
        ingest.add_new_accession_uc_table(new_accessions)
        ingest.add_fasta_uc_table()
    
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
                            type="secondary", 
                            icon=":material/download:",
                            data=results_blast,
                            file_name="results_blast.xlsx",
                            on_click="ignore",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
            st.download_button(label="Download EffectorP results", 
                            type="secondary", 
                            icon=":material/download:",
                            data=results_effectorp,
                            file_name="results_effectorp.xlsx",
                            on_click="ignore",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
            st.download_button(label="Download InterProScan PFAM Domain Search results", 
                            type="secondary", 
                            icon=":material/download:",
                            data=results_pfam,
                            file_name="results_pfam.xlsx",
                            on_click="ignore",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
            st.download_button(label="Download Molecular Weight results", 
                            type="secondary", 
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
                            type="secondary", 
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
                            type="secondary", 
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
                            type="secondary", 
                            icon=":material/download:",
                            data=results_mw,
                            file_name="results_molecularweight.xlsx",
                            mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

        
    # Step 2: Generate optional FASTA file
    # st.markdown("<br><p style='font-size: 24px;'>Following files are now available for download:</p>", unsafe_allow_html=True)

    # if st.button(label="Yes", type="primary"):
    #     df_fasta = dbh.
    #     fasta_file_content, total_accessions, results = fs.generate_fasta_file(accessions)
    #     st.write(results)

    #     # Save to session state
    #     st.session_state["fasta_file_content"] = fasta_file_content
    #     st.session_state["total_accessions"] = total_accessions

    #     # Save locally for BLAST
    #     with open("sequences.fasta", "w") as f:
    #         f.write(fasta_file_content.getvalue().decode("utf-8"))

    # # Retrieve from session state
    # fasta_file_content = st.session_state.get("fasta_file_content", None)
    # total_accessions = st.session_state.get("total_accessions", 0)

    # if fasta_file_content:
    #     st.download_button(label="Download FASTA file",
    #                        type="primary",
    #                        data=fasta_file_content.getvalue(),
    #                        file_name='out_sequences.fasta',
    #                        mime='text/plain')

    # # # Step 3: BLAST sequences
    # # st.markdown("<br><p style='font-size: 24px;'>Step 3: Perform BLAST on retrieved FASTA sequences to get top hits</p>", unsafe_allow_html=True)

    # # option = st.selectbox("How many sequences from the generated FASTA file would you like to process?",
    # #                       ["None", "First n sequences", "All sequences"],
    # #                       index=0, placeholder="Choose an option")

    # # st.markdown("<p style='color: #FF4B4B;'>(Note: The more sequences you select, the more time it will take to BLAST)</p><br>", unsafe_allow_html=True)

    # # if option == "First n sequences":
    # #     num_seq = st.number_input(label="Enter the number of first n sequences required to BLAST", value=0, step=1)
    # # elif option == "All sequences":
    # #     num_seq = total_accessions
    # # else:
    # #     num_seq = 0

    # # num_hits = st.number_input(label="Enter the number of top hits required for each sequence", value=0, step=1)

    # # if st.button(label="Click to BLAST the FASTA sequences", type="primary"):
    # #     if fasta_file_content:
    # #         with st.spinner("Running BLAST... Please wait"):
    # #             start_time = time.time()
    
    # #             # Perform BLAST
    # #             df = bs.generate_blast_dataframe("sequences.fasta", num_seq, num_hits)
    
    # #             execution_time = time.time() - start_time
    # #             st.success(f"BLAST completed in {round(execution_time/60, 2)} minutes")
                
    # #             st.write(df)
    
    # #             # Save df to Excel
    # #             output = BytesIO()
    # #             with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
    # #                 df.to_excel(writer, index=False, sheet_name='BLAST Results')
    # #             output.seek(0)
    
    # #             st.download_button(label="Download BLAST Results as Excel",
    # #                                data=output,
    # #                                file_name="blast_file.xlsx",
    # #                                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

    # #     else:
    # #         st.error("FASTA file not generated. Please generate the FASTA file first.")


    # # # Step 4: PFAM domain search via interpro scan
    # # st.markdown("<br><p style='font-size: 24px;'>Step 4: Perform PFAM domain search on retrieved FASTA sequences via InterProScan</p>", unsafe_allow_html=True)

    # # if st.button(label="Click to search PFAM domains of the FASTA sequences", type="primary"):
    # #     if fasta_file_content:
    # #         with st.spinner("Searching PFAM domains... Please wait"):
    # #             start_time = time.time()
    
    # #             # PFAM domain search
    # #             df = pf.generate_pfam_dataframe("sequences.fasta")

    # #             execution_time = time.time() - start_time
    # #             st.success(f"PFAM domain search completed in {round(execution_time/60, 2)} minutes")

    # #             st.write(df)

    # #             # Save df to Excel
    # #             output = BytesIO()
    # #             with pd.ExcelWriter(output, engine='xlsxwriter') as writer:
    # #                 df.to_excel(writer, index=False, sheet_name='BLAST Results')
    # #             output.seek(0)
    
    # #             st.download_button(label="Download PFAM Domain Search Results as Excel",
    # #                                data=output,
    # #                                file_name="pfam_file.xlsx",
    # #                                mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")
            
    
    # # Reset
    # st.markdown("<br><p style='font-size: 20px;'>Reset</p>", unsafe_allow_html=True)
    # if st.button("Click to reset and start again", type="secondary"):
    #     st.session_state.clear()


if __name__ == "__main__":
    main()
