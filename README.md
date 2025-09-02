# protein_sequence_analysis
This is a streamlit app developed to perform analysis of protein sequences using user provided accession numbers. 
It allows users to perform different types of analyses of protein sequences within a single centralized web application, without having the need to go to different research websites available in the field of Biotechnology.

1. Retrieves FASTA sequences using accession numbers provided by user using Biopython's Entrez module which searches NCBI's database to get FASTA sequences from accession numbers. 
2. Perform BLAST of the retrieved FASTA sequences to get top hits using NCBI's website https://blast.ncbi.nlm.nih.gov/Blast.cgi?PROGRAM=blastp.
3. Perform PFAM domain search

A Medallion Architecture has been implemented on Databricks to incrementally and progressively improve the data structure and quality as it flows through each layer of the architecture:
1. **raw layer** - stores the accession numbers input by user from streamlit web app's UI, their generated FASTA sequences, corresponding BLAST sequences, EffectorP Prediction of input and BLAST sequences, PFAM domain search results of input and BLAST sequences in a flat table structure.

Databricks Unity Catalog tables are implemented to store the data for processed sequences, 


