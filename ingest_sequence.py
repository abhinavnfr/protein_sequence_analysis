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
import re
import time
from bs4 import BeautifulSoup


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
        cursor.close()
        conn.close()
        return new_accessions
    except Exception as e:
        st.error(f"Error filtering new sequences: {str(e)}")
        return []


# fetch FASTA sequence from accession number
def fetch_fasta_sequence(accession, blast_accession):
    Entrez.email = get_entrez_email()
    try:
        with st.spinner(f"Fetching FASTA sequence for accession: {accession}", show_time=True):
            with Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text") as handle:
                record = handle.read()
            if accession == blast_accession:
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
                blasted_sequence.append(fetch_fasta_sequence(hit['accession'], accession))
            st.success(f"Successfully BLASTed sequence for accession: {accession}")
            return blasted_sequence

    except Exception as e:
        st.error(f"Failed to BLAST sequence for {accession}: {str(e)}")


# perform Effector P of protein sequence
def predict_effectorp(accession, blasted_sequence):
    with st.spinner(f"Predicting effectorp for accession: {accession}", show_time=True):
        try:
            # Submit sequence to the form endpoint. The HTML shows the textarea name is 'seq'
            submit_url = st.secrets["effectorp"]["submit_url"]
            
            effectorp_sequence = blasted_sequence
            
            for i in range(5, 22, 4):
                data = {
                    'seq': blasted_sequence[i],
                    'submit': 'Run EffectorP'
                }
                session = requests.Session()
                response = session.post(submit_url, data=data)
                if response.status_code != 200:
                    raise Exception(f"Form submission failed for accession: {blasted_sequence[0]}")
                
                # Regex pattern for UUID (dataset id)
                pattern = r'dataset ([0-9a-f\-]{36})'
                match = re.search(pattern, response.text)
                if match:
                    dataset_id = match.group(1)
                else:
                    raise Exception("Dataset ID not found")

                result_url = f"{st.secrets["effectorp"]["result_url"]}{dataset_id}"
                while True:
                    result_response = session.get(result_url)
                    result_soup = BeautifulSoup(result_response.text, "html.parser")
                    if "Summary table" in result_response.text:
                        break
                    # Wait and retry
                    time.sleep(5)
                
                # Extract the summary table from final result page
                table = result_soup.find("table")
                results = []
                for tr in table.find_all("tr"):
                    row = [td.get_text(strip=True) for td in tr.find_all(["td", "th"])]
                    if row:
                        results.append(row)
                effectorp_sequence += results[1][1:5]
            st.success(f"EffectorP predicted for accession: {accession}")
            return effectorp_sequence
        
        except Exception as e:
            st.error(f"Failed to predict EffectorP for accession: {accession}")


# submit protein sequences to the InterProScan REST API.
def submit_to_interpro(sequence):
    url = f'{st.secrets["interpro"]["url"]}run/'
    data = {
            "sequence": sequence.strip(),
            "email": get_entrez_email(), # Replace with a valid email address
            "title": "InterProScan job", # Optional title for your job
            }
    response = requests.post(url, data=data)

    # Check response status
    if response.status_code == 200:
        return response.text.strip()  # Return the job ID
    else:
        raise Exception(f"Error submitting job: {response.status_code} {response.text}")


# check the status of the InterProScan job.
def check_status(job_id):
    url = f'{st.secrets["interpro"]["url"]}status/{job_id}'
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        raise Exception(f"Error checking status: {response.text}")


# retrieve results for a completed InterProScan job.
def retrieve_results(job_id):
    url = f'{st.secrets["interpro"]["url"]}result/{job_id}/tsv'
    response = requests.get(url)
    if response.status_code == 200:
        return response.text
    else:
        raise Exception(f"Error retrieving results: {response.text}")


# perform Interpro Scan PFAM Domain
def pfam_domain_search(accession, effectorp_sequence):
    try:
        with st.spinner(f"Performing Interpro Scan PFAM Domain search for accession: {accession}", show_time=True):
            # Step 1: Submit sequence to InterPro
            job_id = submit_to_interpro(effectorp_sequence[1])
            st.write(f"Submitted Job ID for accession: {accession}: {job_id}")
    
            # Step 2: Check job status
            while True:
                status = check_status(job_id)
                st.write(f"Job status: {status}") # optional print statement
                if status == "FINISHED":
                    break
                elif status == "ERROR":
                    raise Exception(f"Error occurred during InterProScan job for accession: {accession}")
                time.sleep(30)  # Wait for 30 seconds before checking again
    
            # Step 3: Retrieve and process results
            results = retrieve_results(job_id)
            pfam_sequence = effectorp_sequence
            for line in results.strip().split("\n"):
                if line.startswith("#"):  # Ignore comment lines
                    continue
                parts = line.split("\t")
                database = parts[3]
                domain_acc = parts[4]  # Domain accession number (e.g., PF00931)
                domain_name = parts[5]  # Domain name (e.g., NB-ARC)

                # Filter only PFAM domains
                if database == "Pfam":
                    pfam_sequence.append(domain_acc)
                    pfam_sequence.append(domain_name)
            
            st.success(f"Completed Interpro Scan PFAM Domain search for accession: {accession}")
            return pfam_sequence
    except Exception as e:
        st.error(f"Failed to perform PFAM Domain search for {accession}: {str(e)}")


# update UC table raw.accession
def update_uc_table_accession(pfam_sequence: list) -> None:
    uc_table = "workspace.raw.accession"
    with st.spinner(f"Inserting processed sequence to UC table {uc_table} for accession: {pfam_sequence[0]}", show_time=True):
        try:
            conn = dbh.get_databricks_connection()
            cursor = conn.cursor()

            # Fetch column names for the table
            cursor.execute(f"DESCRIBE TABLE {uc_table}")
            columns_info = cursor.fetchall()
            table_columns = [row[0] for row in columns_info if row != ""]
 
            # Trim to only number of provided values
            insert_columns = table_columns[:len(pfam_sequence)] + ["record_create_ts"]

            # Prepare parameter placeholders (use ? for Databricks SQL)
            placeholders = ", ".join(["?"] * len(pfam_sequence) + ["current_timestamp()"])
            col_names = ", ".join(insert_columns)

            query = f"INSERT INTO {uc_table} ({col_names}) VALUES ({placeholders})"
            cursor.execute(query, pfam_sequence)

            conn.commit()
            cursor.close()
            conn.close()

            st.success(f"Processed sequence inserted into UC table {uc_table} for accession: {pfam_sequence[0]}")

        except Exception as e:
            st.error(f"Error inserting processed sequence for accession {pfam_sequence[0]} into UC table {uc_table}: {e}")
  
    


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


