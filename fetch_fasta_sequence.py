from Bio import Entrez, SeqIO
import streamlit as st
from io import BytesIO
import databricks_handler as dbh


def fetch_fasta_sequence(accession: str) -> str:
    """
    Retrieves a FASTA sequence from the NCBI database.
    Parameters:
        accession (str): Accession number of the sequence.
    Returns:
        str: FASTA sequence.
    """
    Entrez.email = dbh.get_entrez_email()

    # Search for the accession in the NCBI database
    try:
        handle = Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text")
        record = handle.read()
        handle.close()
        return record
    except Exception as e:
        return f"Failed to retrieve {accession}: {str(e)}"


def update_uc_table_accession(accessions: list) -> None:
    """
    Updates UC table workspace.raw.accession with accession numbers and their corresponding FASTA sequences.
    Parameters:
        accessions (list): List of accession numbers.
    Returns:
        None
    """
    uc_table = "workspace.raw.accession"
    try:
        conn = dbh.get_databricks_connection()
        cursor = conn.cursor()
        
        # Fetch existing accession IDs
        cursor.execute(f"SELECT id FROM {uc_table}")
        existing_ids = set(row[0] for row in cursor.fetchall())
        
        new_accessions = [acc for acc in accessions if acc not in existing_ids]
        
        # Insert only new ones
        new_accessions_count = len(new_accessions)
        progress_bar = st.progress(0)  # Initialize the progress bar
        percentage_text = st.empty()  # Placeholder for percentage text
        if new_accessions:
            for acc in new_accessions:
                sequence = fetch_fasta_sequence(acc)
                cursor.execute(f"INSERT INTO {uc_table} (id, sequence, record_create_ts) VALUES (?)", (acc, sequence, CURRENT_TIMESTAMP()))
                progress_percentage = (i + 1) / new_accessions_count  # Calculate percentage
                progress_bar.progress(progress_percentage)  # Update the progress bar
                percentage_text.text(f"Progress: {int(progress_percentage * 100)}%")  # Update percentage text
            conn.commit()
        
        cursor.close()
        conn.close()
        st.success(f"New accession numbers added to UC table {uc_table}: {new_accessions_count}")
        st.success(f"Input accession numbers already existing in UC table {uc_table}: {len(accessions)-new_accessions_count}")
        st.success(f"Total accession numbers now in UC table {uc_table}: {len(existing_ids)+len(new_accessions)}")
    
    except Exception as e:
        st.error(f"Error updating UC table {uc_table}: {str(e)}")


def generate_fasta_file(accessions: list):
    """
    Generates a FASTA file from a list of accession numbers provided.
    Parameters:
        accessions (list): List of accession numbers.
    Returns:
        BytesIO: In-memory FASTA file containing the retrieved sequences.
    """
    uc_table = "workspace.raw.accession"
    results = {}
    total_accessions = len(accessions)

    try:
        conn = dbh.get_databricks_connection()
        cursor = conn.cursor()

        # for acc accessions:
        #     results[acc] = 

        # Save the retrieved sequences into a BytesIO object
        fasta_io = BytesIO()
        for accession, sequence in results.items():
            fasta_io.write(f"{sequence}\n".encode('utf-8'))

        fasta_io.seek(0)  # Reset the file pointer to the beginning of the BytesIO object
        return fasta_io, total_accessions, results
    
    except Exception as e:
        st.error(f"Error generating FASTA file: {str(e)}")
