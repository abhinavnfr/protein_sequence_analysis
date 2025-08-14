from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
import streamlit as st
from io import BytesIO
import databricks_handler as dbh
from datetime import datetime
import os
import pandas as pd
import numpy as np
import requests
import time


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
        st.success(f"New accessions found: {len(new_accessions)}")
        return new_accessions
    except Exception as e:
        st.error(f"Error filtering new sequences: {str(e)}")
        return []


# fetch FASTA sequence from accession number
def fetch_fasta_sequence(accession: str) -> str:
    Entrez.email = get_entrez_email()
    try:
        with st.spinner(f"Fetching FASTA sequence for accession: {accession}", show_time=True):
            with Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text") as handle:
                record = handle.read()
            st.success(f"Successfully fetched FASTA sequence for accession: {accession}")
            return record
    except Exception as e:
        st.error(f"Failed to retrieve sequence for {accession}: {str(e)}")


# BLAST a FASTA sequence
def blast_sequence(accession, fasta_sequence, num_hits=5):
    try:
        with st.spinner(f"BLASTing FASTA sequence for accession: {accession}", show_time=True):
            # Run BLAST against the NCBI database
            result_handle = NCBIWWW.qblast("blastp", "nr", str(fasta_sequence))

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

            blasted_sequence = [accession, fasta_sequence]
            for i, hit in enumerate(top_hits, start=1):
                blasted_sequence.append(f"{hit['percent_identity']:.2f}%")
                blasted_sequence.append(hit['scientific_name'])
                blasted_sequence.append(hit['accession'])
            st.success(f"BLASTed sequence for accession: {accession}")
            st.write(blasted_sequence)
            return blasted_sequence

    except Exception as e:
        st.error(f"Failed to BLAST sequence for {accession}: {str(e)}")


# Submit protein sequences to the InterProScan REST API.
def submit_to_interpro(blasted_sequence):
    url = st.secrets["interpro"]["url"]
    with open(fasta_file, "rb") as fasta:
        # Include both email and sequence file in the POST request
        files = {"sequence": fasta}
        data = {
            "email": get_entrez_email(),  # Replace with a valid email address
            "title": "InterProScan job",        # Optional title for your job
        }
        response = requests.post(url, files=files, data=data)

        # Check response status
        if response.status_code == 200:
            return response.text.strip()  # Return the job ID
        else:
            raise Exception(f"Error submitting job: {response.status_code} {response.text}")


def check_status(job_id):
    """
    Check the status of the InterProScan job.

    Parameters:
        job_id (str): The job ID.

    Returns:
        str: The job status (e.g., "RUNNING", "FINISHED", "ERROR").
    """
    url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/status/{job_id}"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        raise Exception(f"Error checking status: {response.text}")


def retrieve_results(job_id):
    """
    Retrieve results for a completed InterProScan job.

    Parameters:
        job_id (str): The job ID.

    Returns:
        str: Results in TSV format.
    """
    url = f"https://www.ebi.ac.uk/Tools/services/rest/iprscan5/result/{job_id}/tsv"
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        raise Exception(f"Error retrieving results: {response.text}")


def process_results(results, fasta_file):
    """
    Process the InterProScan results to extract PFAM domains with names,
    and dynamically parse Accession and Sequence_Name.

    Parameters:
        results (str): Results in TSV format.
        fasta_file (str): Path to the input FASTA file.

    Returns:
        pd.DataFrame: DataFrame with sequence details and PFAM domain names.
    """
    rows = []
    for line in results.strip().split("\n"):
        if line.startswith("#"):  # Ignore comment lines
            continue
        parts = line.split("\t")
        sequence_name = parts[0]
        database = parts[3]
        domain_acc = parts[4]  # Domain accession number (e.g., PF00931)
        domain_name = parts[5]  # Domain name (e.g., NB-ARC)

        # Filter only PFAM domains
        if database == "Pfam":
            rows.append((sequence_name, domain_acc, domain_name))

    # Aggregate domains for each sequence
    domain_dict = {}
    for sequence_name, domain_acc, domain_name in rows:
        if sequence_name not in domain_dict:
            domain_dict[sequence_name] = []
        domain_dict[sequence_name].append((domain_acc, domain_name))

    # Convert to DataFrame
    data = []
    for record in SeqIO.parse(fasta_file, "fasta"):
        # Split the description to separate Accession and Sequence_Name
        accession = record.id  # First part is the accession (e.g., XP_042375699.1)
        sequence_name = " ".join(record.description.split(" ")[1:])  # Everything after the accession
        sequence = str(record.seq)
        domains = domain_dict.get(record.id, [])
        domain_names = [f"{domain_acc} ({domain_name})" for domain_acc, domain_name in domains]
        row = [accession, sequence_name, sequence] + domain_names
        data.append(row)

    # Create column names
    max_domains = max(len(row) - 3 for row in data)  # Find max domains
    columns = ["Accession", "Sequence_Name", "Sequence"] + [f"PFAM_Domain_{i+1}" for i in range(max_domains)]
    return pd.DataFrame(data, columns=columns)


# def split_fasta(input_file, chunk_size=750):
#     output_prefix = input_file.replace(".fasta", "_chunk")
#     chunk_files = []

#     with open(input_file, "r") as f:
#         sequences = []
#         current_seq = ""
#         header = ""

#         for line in f:
#             line = line.strip()
#             if line.startswith(">"):
#                 if header and current_seq:
#                     sequences.append((header, current_seq))
#                 header = line
#                 current_seq = ""
#             else:
#                 current_seq += line

#         if header and current_seq:
#             sequences.append((header, current_seq))

#         for i in range(0, len(sequences), chunk_size):
#             chunk = sequences[i:i + chunk_size]
#             output_file = f"{output_prefix}_part{i//chunk_size + 1}.fasta"
#             chunk_files.append(output_file)
#             with open(output_file, "w") as out_f:
#                 for h, seq in chunk:
#                     out_f.write(f"{h}\n{seq}\n")

#     return chunk_files


def generate_pfam_dataframe(input_fasta_file):
    # chunk_files = split_fasta(input_fasta_file)

    for fasta_file in chunk_files:
        df_list = []
        # Step 1: Submit sequences to InterPro
        job_id = submit_to_interpro(fasta_file)
        st.write(f"Job submitted. Job ID: {job_id}")
    
        # Step 2: Check job status
        while True:
            status = check_status(job_id)
            # st.write(f"Job status: {status}") # optional print statement
            if status == "FINISHED":
                st.write(f"Job status: {status}")
                break
            elif status == "ERROR":
                raise Exception("Error occurred during InterProScan job.")
            time.sleep(30)  # Wait for 30 seconds before checking again
    
        # Step 3: Retrieve and process results
        results = retrieve_results(job_id)
        df = process_results(results, fasta_file)
        df_list.append(df)

    df_final = pd.concat(df_list, ignore_index=True)
    
    return df_final

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






