import os
import streamlit as st
import pandas as pd
from io import BytesIO
import time
import ingest_sequence as ingest
# import blast_sequence as bs
# import pfam_domains_interpro_scan as pf
import databricks_handler as dbh


def main():
    st.markdown("<h1>Protein Sequence Analysis App</h1><br>", unsafe_allow_html=True)
    st.markdown("<p style='text-align: right; color: #FF4B4B;'>by Abhinav Rana</p>", unsafe_allow_html=True)

    # Step 1: Process input accession numbers file
    st.markdown("<br><p style='font-size: 24px;'>To get started, choose a text file containing accession numbers</p><br>", unsafe_allow_html=True)
    input_file = st.file_uploader(label="Upload file", type=["txt"])
    if input_file is not None:
        accessions = [line.strip() for line in input_file.read().decode("utf-8").splitlines()]
        new_accessions = ingest.filter_new_sequences(accessions)
        new_accesions_count = len(new_accessions)

        if st.button(label="Process sequence end-to-end, i.e., perform BLAST, InterProScan, EffectorP", type="primary"):
            # initiate progress bar
            progress_bar = st.progress(0)
            status_text = st.empty()
            
            for i, acc in enumerate(new_accessions):
                try:
                    status_text.text(f"Processing {i+1}/{new_accesions_count}: {acc}")
                    fasta_sequence = ingest.fetch_fasta_sequence(acc, acc)
                    blasted_sequence = ingest.blast_sequence(acc, fasta_sequence)
                    ingest.add_blast_uc_table(acc, blasted_sequence)
                except Exception as e:
                    st.error(f"Error processing {acc}: {e}")
                
                # update progress bar
                progress_percent = int(((i+1) / new_accesions_count) * 100)
                progress_bar.progress(progress_percent)
                time.sleep(0.1)
            
            time.sleep(1)
            ingest.predict_effectorp()
            time.sleep(1)
            ingest.pfam_domain_search()
            time.sleep(1)
            ingest.calculate_molecular_weight_kda()
            progress_bar.progress(100)
    
        elif st.button(label="Perform EffectorP alone", type="primary"):
            # initiate progress bar
            progress_bar = st.progress(0)
            status_text = st.empty()

            for i, acc in enumerate(new_accessions):
                try:
                    status_text.text(f"Processing {i+1}/{new_accesions_count}: {acc}")
                    fasta_sequence = ingest.fetch_fasta_sequence(acc, acc)
                    blasted_sequence = ingest.blast_sequence(acc, fasta_sequence)
                    ingest.add_blast_uc_table(acc, blasted_sequence)
                except Exception as e:
                    st.error(f"Error processing {acc}: {e}")
                
                # update progress bar
                progress_percent = int(((i+1) / new_accesions_count) * 100)
                progress_bar.progress(progress_percent)
                time.sleep(0.1)
            
            time.sleep(1)
            ingest.predict_effectorp()

    
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
            
    
    # Reset
    st.markdown("<br><p style='font-size: 20px;'>Reset</p>", unsafe_allow_html=True)
    if st.button("Click to reset and start again", type="secondary"):
        st.session_state.clear()


if __name__ == "__main__":
    main()
