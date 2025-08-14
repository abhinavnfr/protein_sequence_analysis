import os
import pandas as pd
import numpy as np
from Bio.Blast import NCBIWWW, NCBIXML
from Bio import SeqIO
import streamlit as st



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
