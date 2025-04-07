import requests
import time
import pandas as pd
from Bio import SeqIO
import streamlit as st


def submit_to_interpro(fasta_file):
    """
    Submit protein sequences to the InterProScan REST API.

    Parameters:
        fasta_file (str): Path to the FASTA file containing sequences.

    Returns:
        str: Job ID for the submitted job.
    """
    url = "https://www.ebi.ac.uk/Tools/services/rest/iprscan5/run/"
    with open(fasta_file, "rb") as fasta:
        # Include both email and sequence file in the POST request
        files = {"sequence": fasta}
        data = {
            "email": "abhinavrana18july@gmail.com",  # Replace with a valid email address
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


def split_fasta(input_file, chunk_size=750):
    output_prefix = input_file.replace(".fasta", "_chunk")
    chunk_files = []

    with open(input_file, "r") as f:
        sequences = []
        current_seq = ""
        header = ""

        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if header and current_seq:
                    sequences.append((header, current_seq))
                header = line
                current_seq = ""
            else:
                current_seq += line

        if header and current_seq:
            sequences.append((header, current_seq))

        for i in range(0, len(sequences), chunk_size):
            chunk = sequences[i:i + chunk_size]
            output_file = f"{output_prefix}_part{i//chunk_size + 1}.fasta"
            chunk_files.append(output_file)
            with open(output_file, "w") as out_f:
                for h, seq in chunk:
                    out_f.write(f"{h}\n{seq}\n")

    return chunk_files


def generate_pfam_dataframe(input_fasta_file):
    chunk_files = split_fasta(input_fasta_file)

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

# df.shape[0]

# # Step 4: Save to Excel file or display
# file_name = "output_615.xlsx"

# df.to_excel(file_name, index=False)
# print(f"DataFrame stored in {file_name}")
