import os
import streamlit as st
import fetch_fasta_sequence as fs
import blast_sequence as bs


def main():
    st.markdown("<h1 style='color: white;'>Protein Sequence Analysis App</h1><br>", unsafe_allow_html=True)
    st.markdown("<p style='text-align: right; color: #ADD8E6;'>by Abhinav Rana</p>", unsafe_allow_html=True)
    
    # file uploader
    st.markdown("<br><p style='font-size: 24px; color: white;'>Step 1: To get started, choose a text file containing accession numbers</p><br>", unsafe_allow_html=True)
    input_file = st.file_uploader(label="Upload file", type=["txt"])

    # fetch FASTA sequences
    st.markdown("<br><p style='font-size: 24px; color: white;'>Step 2: Fetch FASTA sequences from accession numbers</p>", unsafe_allow_html=True)
    fasta_file = None # initialize fasta_file
    total_accessions = 0 # initialize total_accessions
    if st.button(label="Click to generate FASTA sequences", type="primary"):
            fasta_file, total_accessions = fs.generate_fasta_file(input_file)
    if fasta_file:
            st.download_button(label="Download FASTA file having sequences", data=fasta_file, file_name='sequences.fasta', mime='text/plain')
                

    # perform BLAST on retrieved FASTA sequences to get top hits
    st.markdown("<br><p style='font-size: 24px; color: white;'>Step 3: Perform BLAST on retrieved FASTA sequences to get top hits</p>", unsafe_allow_html=True)
    option = st.selectbox("How many sequences from the generated FASTA file would you like to process?",["None", "First n sequences", "All sequences"], index=0, placeholder="Choose an option")
    st.markdown("<p style='color: #FF4B4B;'>(Note: The more sequences you select, the more time it will take to BLAST)</p><br>", unsafe_allow_html=True)
    
    if option == "First n sequences":
        num_seq = st.number_input(label="Enter the number of first n sequences required to BLAST", value=0, step=1)
    else:
        num_seq = total_accessions
        
    num_hits = st.number_input(label="Enter the number of top hits required for each sequence", value=0, step=1)
    
    blast_file = None  # initialize blast_file
    if st.button(label="Click to BLAST the FASTA sequences", type="primary"):
        if fasta_file:
            df = bs.generate_blast_dataframe(fasta_file, num_seq, num_hits)
            st.write(df)
        else:
            st.error("FASTA file not generated. Please generate the FASTA file first.")

    # reset
    st.markdown("<br><p style='font-size: 20px; color: white;'>Reset</p>", unsafe_allow_html=True)
    st.button("Click to reset and start again", type="secondary")


if __name__ == "__main__":
    main()
