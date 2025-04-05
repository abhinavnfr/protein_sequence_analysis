import os
import time
import streamlit as st
import pandas as pd
import fetch_fasta_sequence as fs
import blast_sequence as bs


def main():
    st.markdown("<h1>Protein Sequence Analysis App</h1><br>", unsafe_allow_html=True)
    st.markdown("<p style='text-align: right; color: #FF4B4B;'>by Abhinav Rana</p>", unsafe_allow_html=True)

    # Step 1: Upload accession number file
    st.markdown("<br><p style='font-size: 24px;'>Step 1: Upload a text file containing accession numbers</p><br>",
                unsafe_allow_html=True)
    input_file = st.file_uploader(label="Upload file", type=["txt"])

    # Step 2: Fetch FASTA Sequences
    st.markdown("<br><p style='font-size: 24px;'>Step 2: Fetch FASTA sequences from accession numbers</p>",
                unsafe_allow_html=True)
    fasta_file = None
    total_accessions = 0

    if st.button(label="Click to generate FASTA sequences", type="primary"):
        fasta_file, total_accessions = fs.generate_fasta_file(input_file)

    if fasta_file:
        st.download_button(label="Download FASTA file having sequences",
                           data=fasta_file,
                           file_name='sequences.fasta',
                           mime='text/plain')

    # Step 3: Perform BLAST
    st.markdown("<br><p style='font-size: 24px;'>Step 3: Perform BLAST on retrieved FASTA sequences to get top hits</p>",
                unsafe_allow_html=True)

    option = st.selectbox("How many sequences from the generated FASTA file would you like to process?",
                          ["None", "First n sequences", "All sequences"], index=0)

    st.markdown("<p style='color: #FF4B4B;'>(Note: The more sequences you select, the more time it will take to BLAST)</p><br>",
                unsafe_allow_html=True)

    if option == "First n sequences":
        num_seq = st.number_input(label="Enter the number of first n sequences required to BLAST", value=0, step=1)
    else:
        num_seq = total_accessions

    num_hits = st.number_input(label="Enter the number of top hits required for each sequence", value=0, step=1)

    if st.button(label="Click to BLAST the FASTA sequences", type="primary"):
        if fasta_file:
            with st.spinner("Running BLAST... Please wait"):
                start_time = time.time()

                # Save fasta_file to disk for BLAST input
                with open("sequences.fasta", "w") as f:
                    f.write(fasta_file.getvalue().decode())

                # Generate DataFrame
                df = bs.generate_blast_dataframe("sequences.fasta", num_seq, num_hits)

                execution_time = time.time() - start_time
                st.success(f"BLAST completed in {round(execution_time, 2)} seconds!")

                st.dataframe(df)

                # Export DataFrame to Excel
                blast_excel = "blast_file.xlsx"
                df.to_excel(blast_excel, index=False)

                with open(blast_excel, "rb") as f:
                    st.download_button(label="Download BLAST Excel File",
                                       data=f,
                                       file_name="blast_file.xlsx",
                                       mime="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet")

                # Auto delete temp files
                os.remove("sequences.fasta")
                os.remove(blast_excel)

        else:
            st.error("FASTA file not generated. Please generate the FASTA file first.")

    # Step 4: Reset
    st.markdown("<br><p style='font-size: 20px;'>Reset</p>", unsafe_allow_html=True)
    st.button("Click to reset and start again", type="secondary")


if __name__ == "__main__":
    main()
