# protein_sequence_analysis

### This is a streamlit app developed to perform analysis of protein sequences using user provided accession numbers. 
It allows users to perform different types of analyses of protein sequences within a single centralized web application, without having the need to go to multiple different Bioinformatics tools available in the field of Biotechnology.
1. Retrieve FASTA sequences using accession numbers provided using https://www.ncbi.nlm.nih.gov/
2. Perform BLAST of the retrieved FASTA sequences to get top hits using https://blast.ncbi.nlm.nih.gov/Blast.cgi
3. Perform EffectorP search to predict fungal and oomycete effector proteins using https://effectorp.csiro.au/
4. Perform InterPro Scan PFAM domain search of protein sequences using https://www.ebi.ac.uk/interpro/search/sequence/

### Medallion Data Architecture and Databricks Unity Catalog Tables are used to store this high volume of data obtained from multiple different Bioinformatics tools in an organized manner.
1. Raw layer - holds raw data ingested from these different Bioinformatics web application.
2. Curated layer - holds cleaned, transformed and standardized data, curated from raw data.
3. Analytical layer - can hold customized views created based on user requirements from curated data, as and when required.
