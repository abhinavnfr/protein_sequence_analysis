from Bio import Entrez, SeqIO
from Bio.Blast import NCBIWWW, NCBIXML
from Bio.SeqUtils.ProtParam import ProteinAnalysis
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


# Filter only new sequences from input
def filter_new_sequences(accessions: list) -> list:
    if len(accessions) == 0:
        return []
    
    uc_table = "workspace.raw.protein"
    
    try:
        conn = dbh.get_databricks_connection()
        cursor = conn.cursor()
        cursor.execute(f"SELECT id FROM {uc_table}")
        existing_ids = set(row[0] for row in cursor.fetchall())
        new_accessions = [acc for acc in accessions if acc not in existing_ids]
        st.success(f"New accessions to add into UC table {uc_table}: {len(new_accessions)}")
        cursor.close()
        conn.close()
        return new_accessions
    
    except Exception as e:
        st.error(f"Error filtering new sequences: {str(e)}")
        return []


# Add new accession IDs into UC table raw.protein
def add_new_accession_uc_table(accessions: list):
    if len(accessions) == 0:
        return
    
    uc_table = "workspace.raw.protein"
    
    try:
        with st.spinner(f"Adding new accession IDs into UC table {uc_table}", show_time=True):
            conn = dbh.get_databricks_connection()
            cursor = conn.cursor()
            # Fetch column names for the table
            cursor.execute(f"DESCRIBE TABLE {uc_table}")
            columns_info = cursor.fetchall()
            table_columns = [row[0] for row in columns_info if row != ""]

            for acc in accessions:
                query = f"INSERT INTO {uc_table} (id, record_create_ts, record_update_ts) VALUES (?, current_timestamp(), current_timestamp())"
                cursor.execute(query, (acc,))

            conn.commit()
            cursor.close()
            conn.close()
            st.success(f"Successfully added new accessions {accessions} into UC table {uc_table}")

    except Exception as e:
        st.error(f"Error adding new accessions {accessions} into UC table {uc_table}: {e}")


# Fetch FASTA sequence from accession number
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


# Add fasta sequence into UC table raw.protein
def add_fasta_uc_table():
    uc_table = "workspace.raw.protein"

    try:
        with st.spinner(f"Adding FASTA sequences for new accessions into UC table {uc_table}", show_time=True):
            conn = dbh.get_databricks_connection()
            cursor = conn.cursor()
            cursor.execute(f"SELECT id FROM {uc_table} WHERE fasta_sequence IS NULL")
            accessions = set(row[0] for row in cursor.fetchall())
            accession_lst = [id for id in accessions]

            if len(accession_lst) == 0:
                return
            
            for acc in accession_lst:
                fasta_sequence = fetch_fasta_sequence(acc, acc)
                update_query = f"""UPDATE {uc_table} 
                                        SET fasta_sequence = '{fasta_sequence}', 
                                        record_update_ts = current_timestamp() 
                                    WHERE id = '{acc}'
                                """
                cursor.execute(update_query)
            
            conn.commit()
            cursor.close()
            conn.close()
            st.success(f"Successfully added FASTA sequences for new accessions into UC table {uc_table}")
    
    except Exception as e:
        st.error(f"Failed to add FASTA sequences for new accessions into UC table {uc_table}: {str(e)}")


# Get sequences to BLAST
def get_seq_to_blast():
    uc_table = "workspace.raw.protein"
    
    try:
        with st.spinner(f"Identifying sequences to BLAST from UC table {uc_table}"):
            conn = dbh.get_databricks_connection()
            cursor = conn.cursor()
            cursor.execute(f"SELECT id, fasta_sequence FROM {uc_table} WHERE id NOT IN (SELECT DISTINCT blast_of_id FROM {uc_table} WHERE blast_of_id IS NOT NULL)")
            rows = cursor.fetchall()
            fasta_sequences = [(row[0], row[1]) for row in rows]
            st.success(f"Successfully identified sequences to BLAST from UC table {uc_table}")
            return fasta_sequences
    
    except Exception as e:
        st.error(f"Failed to identify sequences to BLAST from UC table {uc_table}: {str(e)}")


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

                # Get the top n unique percentage identities
                top_hits = []
                unique_identities = set()

                for hit in hits:
                    if len(top_hits) >= num_hits: # This decides how many top hits you want
                        break
                    if hit["percent_identity"] not in unique_identities:
                        top_hits.append(hit)
                        unique_identities.add(hit["percent_identity"])

            blasted_sequence = []
            for i, hit in enumerate(top_hits, start=1):
                temp_list = [hit['accession'], fetch_fasta_sequence(hit['accession'], accession), accession, i, f"{hit['percent_identity']:.2f}%"]
                blasted_sequence.append(temp_list)
            st.success(f"Successfully BLASTed sequence for accession: {accession}")
            return blasted_sequence

    except Exception as e:
        st.error(f"Failed to BLAST sequence for {accession}: {str(e)}")


# add BLAST sequences to UC table raw.protein
def add_blast_uc_table(accession: str, blasted_sequence: list) -> None:
    uc_table = "workspace.raw.protein"
    with st.spinner(f"Adding BLAST sequences into UC table {uc_table} for accession: {accession}", show_time=True):
        try:
            conn = dbh.get_databricks_connection()
            cursor = conn.cursor()

            # Fetch column names for the table
            cursor.execute(f"DESCRIBE TABLE {uc_table}")
            columns_info = cursor.fetchall()
            table_columns = [row[0] for row in columns_info if row != ""]
 
            for seq in blasted_sequence:
                # Trim to only number of provided values
                insert_columns = table_columns[2:len(seq)+2] + ["record_create_ts", "record_update_ts"]

                # Prepare parameter placeholders (use ? for Databricks SQL)
                placeholders = ", ".join(["?"] * len(seq) + ["current_timestamp()", "current_timestamp()"])
                col_names = ", ".join(insert_columns)

                insert_query = f"INSERT INTO {uc_table} ({col_names}) VALUES ({placeholders})"
                cursor.execute(insert_query, seq)

            conn.commit()
            cursor.close()
            conn.close()

            st.success(f"BLAST sequences added into UC table {uc_table} for accession: {accession}")

        except Exception as e:
            st.error(f"Error adding BLAST sequences for accession {accession} into UC table {uc_table}: {e}")
  

# perform Effector P of protein sequence
def predict_effectorp():
    uc_table = "workspace.raw.protein"
    with st.spinner(f"Predicting EffectorP for sequences", show_time=True):
        try:
            conn = dbh.get_databricks_connection()
            cursor = conn.cursor()

            cursor.execute(f"SELECT fasta_sequence FROM {uc_table} WHERE prediction IS NULL")
            sequences = set(row[0] for row in cursor.fetchall())
            blasted_sequence = [seq for seq in sequences]

            if len(blasted_sequence) == 0:
                st.success(f"EffectorP predicted for all sequences")
                return

            # Submit sequence to the form endpoint. The HTML shows the textarea name is 'seq'
            submit_url = st.secrets["effectorp"]["submit_url"]
            
            for seq in blasted_sequence:
                data = {
                    'seq': seq,
                    'submit': 'Run EffectorP'
                }
                session = requests.Session()
                response = session.post(submit_url, data=data)
                if response.status_code != 200:
                    raise Exception(f"Form submission failed for sequence: {seq}")
                
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
                update_sql = f"""UPDATE {uc_table} 
                                    SET cytoplasmic_effector = '{results[1][1]}', 
                                        apoplastic_effector = '{results[1][2]}', 
                                        non_effector = '{results[1][3]}',
                                        prediction = '{results[1][4]}',
                                        record_update_ts = current_timestamp()
                                    WHERE fasta_sequence = '{seq}'
                            """
                cursor.execute(update_sql)
            
            conn.commit()
            cursor.close()
            conn.close()
            
            st.success(f"EffectorP predicted for all sequences")

        except Exception as e:
            st.error(f"Error predicting EffectorP for sequences: {e}")


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


# perform InterproScan PFAM Domain
def pfam_domain_search():
    uc_table = "workspace.raw.protein"

    try:
        with st.spinner(f"Performing InterproScan PFAM Domain Search for sequences", show_time=True):
            conn = dbh.get_databricks_connection()
            cursor = conn.cursor()

            cursor.execute(f"SELECT fasta_sequence FROM {uc_table} WHERE pfam_domain_acc_1 IS NULL")
            sequences = set(row[0] for row in cursor.fetchall())
            blasted_sequence = [seq for seq in sequences]

            if len(blasted_sequence) == 0:
                st.success(f"Completed InterproScan PFAM Domain Search for all sequences")
                return

            for seq in blasted_sequence:
                # Step 1: Submit sequence to InterPro
                job_id = submit_to_interpro(seq)
                # st.write(f"Submitted Job ID for accession: {accession}: {job_id}") # optional print statement
    
                # Step 2: Check job status
                while True:
                    status = check_status(job_id)
                    # st.write(f"Job status: {status}") # optional print statement
                    if status == "FINISHED":
                        break
                    elif status == "ERROR":
                        raise Exception(f"Error occurred during InterProScan job for accession: {accession}")
                    time.sleep(30)  # Wait for 30 seconds before checking again
    
                # Step 3: Retrieve and process results
                results = retrieve_results(job_id)
                domain_num = 1
                for line in results.strip().split("\n"):
                    if line.startswith("#"):  # Ignore comment lines
                        continue
                    parts = line.split("\t")
                    database = parts[3]
                    domain_acc = parts[4]  # Domain accession number (e.g., PF00931)
                    domain_name = parts[5]  # Domain name (e.g., NB-ARC)

                    # Filter only PFAM domains
                    if database == "Pfam":
                        update_query = f"""UPDATE {uc_table}
                                            SET pfam_domain_acc_{domain_num} = '{domain_acc}',
                                                pfam_domain_name_{domain_num} = '{domain_name}',
                                                record_update_ts = current_timestamp()
                                            WHERE fasta_sequence = '{seq}'
                                        """
                        cursor.execute(update_query)
                        domain_num += 1
            
            conn.commit()
            cursor.close()
            conn.close()

            st.success(f"Completed InterproScan PFAM Domain Search for all sequences")
    
    except Exception as e:
        st.error(f"Error performing InterproScan PFAM Domain Search for sequences: {e}")


# calculate molecular weights of protein sequences
def calculate_molecular_weight_kda():
    uc_table = "workspace.raw.protein"
    st.spinner(f"Calculating molecular weights for sequences", show_time=True)
    try:
        conn = dbh.get_databricks_connection()
        cursor = conn.cursor()

        cursor.execute(f"SELECT fasta_sequence FROM {uc_table} WHERE molecular_weight_kda IS NULL")
        sequences = set(row[0] for row in cursor.fetchall())
        blasted_sequence = [seq for seq in sequences]

        if len(blasted_sequence) == 0:
            st.success(f"Calculated Molecular Weights, Isoelectric Points and Lengths for all sequences")

        for seq in blasted_sequence:
            trimmed_seq = ''.join([line.strip() for line in seq.splitlines() if not line.startswith('>')])
            analysis = ProteinAnalysis(trimmed_seq)
            mw_kda = round(analysis.molecular_weight() / 1000, 2)  # Convert Da to kDa
            pi = round(analysis.isoelectric_point(), 2)
            aa_length = len(trimmed_seq)

            update_sql = f"""UPDATE {uc_table} 
                                    SET molecular_weight_kda = {mw_kda}, 
                                        isoelectric_point_pi = {pi}, 
                                        sequence_length = {aa_length},
                                        record_update_ts = current_timestamp()
                                    WHERE fasta_sequence = '{seq}'
                            """
            cursor.execute(update_sql)

        conn.commit()
        cursor.close()
        conn.close()

        st.success(f"Calculated Molecular Weights, Isoelectric Points and Lengths for all identified sequences")
    
    except Exception as e:
        st.error(f"Error calculating Molecular Weights, Isoelectric Points and Lengths for identified sequences: {e}")


