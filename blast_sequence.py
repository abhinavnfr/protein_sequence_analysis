import os
import pandas as pd
import numpy as np
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO

# Function to BLAST a sequence
def blast_sequence(fasta_file, idx, num_hits):
    # Read the protein sequence from the fasta file
    with open(fasta_file, 'r') as file:
        seq_record = list(SeqIO.parse(file, "fasta"))

    if not seq_record:
        print("No protein sequences found in the file.")
        return

    print(f"Running BLAST for the sequence: {seq_record[idx].id}")

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
    print("Top 3 Percentage Identities with Details:")
    for i, hit in enumerate(top_hits, start=1):
        print(f"{i}: {hit['percent_identity']:.2f}% - Scientific Name: {hit['scientific_name']} - Accession: {hit['accession']}")
        lst.append(f"{hit['percent_identity']:.2f}%")
        lst.append(hit['scientific_name'])
        lst.append(hit['accession'])

    print("="*100)

    return lst

def generate_blast_dataframe(fasta_file, num_seq, num_hits):
    # Create an empty dataframe to store the result
    df = pd.DataFrame(columns=["Complete_Input", "Input_Sequence_ID", "Input_Sequence_Description", "Input_Sequence",
                               "Top_1st_Percent_Identity", "Top_1st_Scientific_Name", "Top_1st_Accession",
                               "Top_2nd_Percent_Identity", "Top_2nd_Scientific_Name", "Top_2nd_Accession",
                               "Top_3rd_Percent_Identity", "Top_3rd_Scientific_Name", "Top_3rd_Accession"]
                      )

    # BLAST sequences on loop
    for i in range(num_seq):
      lst = blast_sequence(fasta_file, i, num_hits)
      df = pd.concat([df, pd.DataFrame([lst], columns=df.columns)], ignore_index=True)

    return df
