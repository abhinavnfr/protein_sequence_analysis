import os
import pandas as pd
import numpy as np
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import streamlit as st

# Function to BLAST a sequence
def blast_sequence(fasta_file, idx, num_hits):
    # Read the protein sequence from the fasta file
    with open(fasta_file, 'r') as file:
        seq_record = list(SeqIO.parse(file, "fasta"))

    if not seq_record:
        st.error("No protein sequences found in the file.")
        return

    st.write(f"Running BLAST for the sequence: {seq_record[idx].id}")

    # Run BLAST against the NCBI database
    result_handle = NCBIWWW.qblast("blastp", "nr", str(seq_record[idx].seq))

    # Save the result to an XML file
    result_file = "blast_results.xml"
    with open(result_file, "w") as out_handle:
        out_handle.write(result_handle.read())

    result_handle.close()

    # Parse the BLAST results
    with open(result_file) as result_handle:
        blast_records = NCBIXML.parse(result_handle)

        # Get the first BLAST record
        blast_record = next(blast_records)

        # List to store information about alignments
        hits = []

        # Extract percentage identities and corresponding details from HSPs
        for alignment in blast_record.alignments:
            for hsp in alignment.hsps:
                percent_identity = (hsp.identities / hsp.align_length) * 100
                hits.append({
                    "percent_identity": percent_identity,
                    "scientific_name": alignment.hit_def.split(' >')[0],
                    "accession": alignment.accession
                })

        # Sort hits by percentage identity in descending order
        hits = sorted(hits, key=lambda x: x["percent_identity"], reverse=True) # hits = hits[::-1]

        # Get the top 3 unique percentage identities
        top_hits = []
        unique_identities = set()

        for hit in hits:
            if len(top_hits) >= num_hits: # This decides how many top hits you want
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


def ordinal(n):
    # Returns ordinal suffix for a number e.g., 1 -> 1st, 2 -> 2nd, 3 -> 3rd, 4 -> 4th ...
    return f"{n}{'th' if 11<=n%100<=13 else {1:'st',2:'nd',3:'rd'}.get(n%10, 'th')}"


def generate_blast_dataframe(fasta_file, num_seq, num_hits):
    # Fixed columns
    base_columns = ["Complete_Input", "Input_Sequence_ID", "Input_Sequence_Description", "Input_Sequence"]

    # Dynamically generate columns based on num_hits
    hit_columns = []
    for i in range(1, num_hits + 1):
        suffix = ordinal(i)
        hit_columns.extend([
            f"Top_{suffix}_Percent_Identity",
            f"Top_{suffix}_Scientific_Name",
            f"Top_{suffix}_Accession"
        ])

    all_columns = base_columns + hit_columns

    # Create empty dataframe with dynamic columns
    df = pd.DataFrame(columns=all_columns)

    # BLAST sequences in a loop
    for i in range(num_seq):
        lst = blast_sequence(fasta_file, i, num_hits)
        df = pd.concat([df, pd.DataFrame([lst], columns=df.columns)], ignore_index=True)

    return df
