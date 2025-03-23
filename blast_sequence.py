import os
import streamlit as st
import pandas as pd
import numpy as np
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO


def blast_sequence(fasta_file, idx, num_hits):
    """
    Performs a BLAST search for a given protein sequence from a FASTA file.
    
    Parameters:
        fasta_file (str): Path to the FASTA file containing protein sequences.
        idx (int): Index of the sequence in the FASTA file to be BLASTed.
        num_hits (int): The number of top BLAST hits to return based on percentage identity. The function returns hits with unique percentage identities.
    
    Returns:
        list: A list containing the queried sequence information along with details of the top 3 BLAST hits, including percentage identity, scientific name, and accession number.
    """
    st.write("opening fasta file")
    # read the protein sequence from the fasta file
    with open(fasta_file, 'r') as file:
        seq_record = list(SeqIO.parse(file, "fasta"))
    st.write("file parsed")
    if not seq_record:
        st.write("No protein sequences found in the file.")
        return

    st.write(f"Running BLAST for the sequence: {seq_record[idx].id}")

    # run BLAST against the NCBI database
    result_handle = NCBIWWW.qblast("blastp", "nr", str(seq_record[idx].seq))

    # save the result to an XML file
    result_file = "blast_results.xml"
    with open(result_file, "w") as out_handle:
        out_handle.write(result_handle.read())

    result_handle.close()

    # parse the BLAST results
    with open(result_file) as result_handle:
        blast_records = NCBIXML.parse(result_handle)

        # get the first BLAST record
        blast_record = next(blast_records)

        # list to store information about alignments
        hits = []

        # extract percentage identities and corresponding details from HSPs
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                percent_identity = (hsp.identities / hsp.align_length) * 100
                hits.append({
                    "percent_identity": percent_identity,
                    "scientific_name": alignment.hit_def.split(' >')[0],
                    "accession": alignment.accession
                })

        # sort hits by percentage identity in descending order
        hits = sorted(hits, key=lambda x: x["percent_identity"], reverse=True) # hits = hits[::-1]

        # get the top num_hits unique percentage identities
        top_hits = []
        unique_identities = set()

        for hit in hits:
            if len(top_hits) >= num_hits: # this decides how many top hits you want
                break
            if hit["percent_identity"] not in unique_identities:
                top_hits.append(hit)
                unique_identities.add(hit["percent_identity"])

    lst = [">" + seq_record[idx].id + " " + seq_record[idx].description + "\n" + str(seq_record[idx].seq),
           seq_record[idx].id, seq_record[idx].description, str(seq_record[idx].seq)]
    st.write(f"Top {num_hits} Percentage Identities with Details:")
    for i, hit in enumerate(top_hits, start=1):
        st.write(f"{i}: {hit['percent_identity']:.2f}% - Scientific Name: {hit['scientific_name']} - Accession: {hit['accession']}")
        lst.append(f"{hit['percent_identity']:.2f}%")
        lst.append(hit['scientific_name'])
        lst.append(hit['accession'])

    st.write("="*100)

    return lst

def generate_blast_dataframe(fasta_file, num_seq, num_hits):
    """
    Generates a BLAST results DataFrame for multiple protein sequences and saves it to an Excel file.
  
    Parameters:
        fasta_file (str): Path to the FASTA file containing protein sequences.
        num_seq (int): Number of sequences to process from the FASTA file.
        num_hits (int): Number of top BLAST hits to retrieve (currently not used in the function but can be integrated).
    
    Returns:
        str: The filename of the generated Excel file containing BLAST results.
    """
    st.write("Starting BLAST.....")
    df = pd.DataFrame(columns=["Complete_Input", "Input_Sequence_ID", "Input_Sequence_Description", "Input_Sequence",
                               "Top_1st_Percent_Identity", "Top_1st_Scientific_Name", "Top_1st_Accession",
                               "Top_2nd_Percent_Identity", "Top_2nd_Scientific_Name", "Top_2nd_Accession",
                               "Top_3rd_Percent_Identity", "Top_3rd_Scientific_Name", "Top_3rd_Accession"]
                      )

    for i in range(num_seq):
        lst = blast_sequence(fasta_file, i, num_hits)
        df = pd.concat([df, pd.DataFrame([lst], columns=df.columns)], ignore_index=True)
    
    # save dataframe to an excel file
    excel_file = "blast_results.xlsx"
    df.to_excel(excel_file, index=False)
    st.write("BLAST finished")
    return excel_file
