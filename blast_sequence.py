import os
import pandas as pd
import numpy as np
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO


def blast_sequence(fasta_file, idx, num_hits):
    """
    Perform BLAST search for a given protein sequence from the FASTA file.
    
    Parameters:
        fasta_file (BytesIO): In-memory FASTA file containing protein sequences.
        idx (int): Index of the sequence in the FASTA file to be BLASTed.
        num_hits (int): Number of top BLAST hits to return.
    
    Returns:
        list: Information about the queried sequence and top BLAST hits.
    """
    # Parse the in-memory FASTA file (no need to open a file on disk)
    fasta_file.seek(0)  # Make sure we're at the beginning of the BytesIO object
    seq_record = list(SeqIO.parse(fasta_file, "fasta"))

    if not seq_record:
        print("No protein sequences found in the file.")
        return []

    print(f"Running BLAST for the sequence: {seq_record[idx].id}")

    # Run BLAST against the NCBI database
    result_handle = NCBIWWW.qblast("blastp", "nr", str(seq_record[idx].seq))

    # Parse the BLAST results
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
    hits = sorted(hits, key=lambda x: x["percent_identity"], reverse=True)

    # Get the top `num_hits` unique percentage identities
    top_hits = []
    unique_identities = set()

    for hit in hits:
        if len(top_hits) >= num_hits:
            break
        if hit["percent_identity"] not in unique_identities:
            top_hits.append(hit)
            unique_identities.add(hit["percent_identity"])

    # Prepare the output list with sequence details and top hits
    lst = [">" + seq_record[idx].id + " " + seq_record[idx].description + "\n" + str(seq_record[idx].seq),
           seq_record[idx].id, seq_record[idx].description, str(seq_record[idx].seq)]

    for i, hit in enumerate(top_hits, start=1):
        lst.append(f"{hit['percent_identity']:.2f}%")
        lst.append(hit['scientific_name'])
        lst.append(hit['accession'])

    print("=" * 100)
    
    return lst


def generate_blast_dataframe(fasta_file, num_seq, num_hits):
    """
    Generates a dataframe of top BLAST hits for sequences from a FASTA file.
    
    Parameters:
        fasta_file (BytesIO): In-memory FASTA file containing protein sequences.
        num_seq (int): Number of sequences to process from the FASTA file.
        num_hits (int): Number of top BLAST hits to return for each sequence.
    
    Returns:
        pd.DataFrame: DataFrame containing top BLAST hits for each sequence.
    """
    df = pd.DataFrame(columns=["Complete_Input", "Input_Sequence_ID", "Input_Sequence_Description", "Input_Sequence",
                               "Top_1st_Percent_Identity", "Top_1st_Scientific_Name", "Top_1st_Accession",
                               "Top_2nd_Percent_Identity", "Top_2nd_Scientific_Name", "Top_2nd_Accession",
                               "Top_3rd_Percent_Identity", "Top_3rd_Scientific_Name", "Top_3rd_Accession"])

    for i in range(num_seq):
        lst = blast_sequence(fasta_file, i, num_hits)
        # Only add to dataframe if the list is not empty (in case of errors or empty sequences)
        if lst:
            df = pd.concat([df, pd.DataFrame([lst], columns=df.columns)], ignore_index=True)

    return df
