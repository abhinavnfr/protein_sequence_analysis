import os
import streamlit as st


def main():
    st.markdown("<h1 style='color: black;'>Protein Sequence Analysis App</h1><br>", unsafe_allow_html=True)
    st.markdown("<p style='text-align: right; color: #00BFFF;'>by Abhinav Rana</p><br>", unsafe_allow_html=True)

    # file uploader
    input_file = st.file_uploader(label="To get started, choose a text file containing accession numbers", type=["txt"], accept_multiple_files=False, label_visibility="visible")

    st.markdown("<p style='color: black;'>Step 1: Fetch FASTA sequences from accession numbers</p><br>")
    if st.button(label="Click to fetch FASTA sequences", type="primary"):
            st.write("Hello")
    st.button("Reset", type="secondary")             
                  

if __name__ == "__main__":
    main()
