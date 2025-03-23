import os
import streamlit as st
import fetch_fasta_sequence as fs


def main():
    st.markdown("<h1 style='color: black;'>Protein Sequence Analysis App</h1><br>", unsafe_allow_html=True)
    st.markdown("<p style='text-align: right; color: #00BFFF;'>by Abhinav Rana</p><br>", unsafe_allow_html=True)

    # file uploader
    st.markdown("<p style='font-size: 24px; color: black;'>Step 1: To get started, choose a text file containing accession numbers</p><br>", unsafe_allow_html=True)
    input_file = st.file_uploader(label="Upload file", type=["txt"], accept_multiple_files=False, label_visibility="visible")

    st.markdown("<p></p><br>", unsafe_allow_html=True)

    # fetch FASTA sequences
    st.markdown("<p style='font-size: 24px; color: black;'>Step 2: Fetch FASTA sequences from accession numbers</p>", unsafe_allow_html=True)
    fasta_file = False
    if st.button(label="Click to generate FASTA sequences", type="primary"):
            fasta_file = fs.generate_fasta_file(input_file)
    if fasta_file:
            with open(fasta_file, "rb") as f:
                    st.download_button(label="Download FASTA file having sequences", data=f, file_name='sequences.fasta', mime='text/plain')
                

    # perform BLAST on retrieved FASTA sequences to get top hits
    st.markdown("<br><p style='font-size: 24px; color: black;'>Step 3: Perform BLAST on retrieved FASTA sequences to get top hits</p>", unsafe_allow_html=True)


    # reset
    st.markdown("<p></p><br>", unsafe_allow_html=True)
    st.markdown("<p style='font-size: 20px; color: black;'>Reset</p>", unsafe_allow_html=True)
    st.button("Click to reset and start again", type="secondary")
    
    

if __name__ == "__main__":
    main()
