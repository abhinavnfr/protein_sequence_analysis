from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import streamlit as st
from io import BytesIO
import databricks_handler as dbh
from datetime import datetime
import os
import pandas as pd
import numpy as np


# Set Entrez email from secrets
def get_entrez_email():
    return st.secrets["entrez"]["email"]


# filter only new sequences from input
def filter_new_sequences(accessions: list) -> list:
    uc_table = "workspace.raw.accession"
    try:
        conn = dbh.get_databricks_connection()
        cursor = conn.cursor()

        cursor.execute(f"SELECT id FROM {uc_table}")
        existing_ids = set(row[0] for row in cursor.fetchall())
        new_accessions = [acc for acc in accessions if acc not in existing_ids]
        return new_accessions
    except Exception as e:
        st.error(f"Error filtering new sequences: {str(e)}")
        return []


# fetch FASTA sequence from accession number
def fetch_fasta_sequence(accession: str) -> str:
    Entrez.email = get_entrez_email()
    try:
        with Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text") as handle:
            record = handle.read()
        return record
    except Exception as e:
        st.error(f"Failed to retrieve sequence for {accession}: {str(e)}")


# BLAST a FASTA sequence
def blast_sequence(accession, seq_record, num_hits=5):
    try:
        # Run BLAST against the NCBI database
        result_handle = NCBIWWW.qblast("blastp", "nr", str(seq_record))

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

        lst = [accession, seq_record]
        # st.write(f"Top {num_hits} Percentage Identities with Details:") # optional print statement
        for i, hit in enumerate(top_hits, start=1):
            lst.append(f"{hit['percent_identity']:.2f}%")
            lst.append(hit['scientific_name'])
            lst.append(hit['accession'])
        st.write(lst)
        return lst

    except Exception as e:
        st.error(f"Failed to BLAST sequence for {accession}: {str(e)}")


# # Update UC table workspace.raw.accession
# def update_uc_table_accession(blasted_accessions: list) -> list:
#     uc_table = "workspace.raw.accession"
#     try:
#         conn = dbh.get_databricks_connection()
#         cursor = conn.cursor()

#         progress_bar = st.progress(0)
#         percentage_text = st.empty()

#         if accessions:
#             for i, acc in enumerate(accessions):
#                 sequence = fetch_fasta_sequence(acc)
#                 cursor.execute(
#                     f"""INSERT INTO {uc_table} (id, sequence, record_create_ts) 
#                         VALUES (?, ?, CURRENT_TIMESTAMP)""",
#                     (acc, sequence)
#                 )
#                 progress_percentage = (i + 1) / new_accessions_count
#                 progress_bar.progress(progress_percentage)
#                 percentage_text.text(f"Updating sequences in database: {int(progress_percentage * 100)}%")

#             conn.commit()

#         cursor.close()
#         conn.close()

#         st.success(f"Already existing sequences in UC table {uc_table}: {len(accessions) - new_accessions_count}")
#         st.success(f"New sequences added to UC table {uc_table}: {new_accessions_count}")
#         st.success(f"Total sequences now in {uc_table}: {len(existing_ids) + new_accessions_count}")
#         return new_accessions

#     except Exception as e:
#         st.error(f"Error updating UC table {uc_table}: {str(e)}")
#         return []
   

# # generate FASTA file from accession numbers
# def generate_fasta_file(accessions: list):
#     uc_table = "workspace.curated.sequence"
#     results = {}
#     total_accessions = len(accessions)

#     try:
#         conn = dbh.get_databricks_connection()
#         cursor = conn.cursor()

#         # Retrieve sequences from UC table
#         for acc in accessions:
#             cursor.execute(f"SELECT sequence FROM {uc_table} WHERE id = ?", (acc,))
#             row = cursor.fetchone()
#             if row and row[0]:
#                 results[acc] = row[0]
        
#         cursor.close()
#         conn.close()
        
#         fasta_io = BytesIO()
#         for sequence in results.values():
#             fasta_io.write(f"{sequence.strip()}\n".encode('utf-8'))
        
#         fasta_io.seek(0)
#         st.success("FASTA file generated successfully!")
#         return fasta_io, total_accessions, results

#     except Exception as e:
#         st.error(f"Error generating FASTA file: {str(e)}")
#         return None, 0, {}






