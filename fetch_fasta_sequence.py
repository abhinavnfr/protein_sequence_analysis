from Bio import Entrez
import streamlit as st
from io import BytesIO
import databricks_handler as dbh
from datetime import datetime


def fetch_fasta_sequence(accession: str) -> str:
    """
    Retrieves a FASTA sequence from the NCBI database.
    """
    Entrez.email = dbh.get_entrez_email()
    try:
        with Entrez.efetch(db="protein", id=accession, rettype="fasta", retmode="text") as handle:
            record = handle.read()
        return record
    except Exception as e:
        return f"Failed to retrieve {accession}: {str(e)}"


def update_uc_table_accession(accessions: list) -> None:
    """
    Updates UC table workspace.raw.accession with accession numbers and their corresponding FASTA sequences.
    """
    uc_table = "workspace.raw.accession"
    try:
        conn = dbh.get_databricks_connection()
        cursor = conn.cursor()

        cursor.execute(f"SELECT id FROM {uc_table}")
        existing_ids = set(row[0] for row in cursor.fetchall())

        new_accessions = [acc for acc in accessions if acc not in existing_ids]
        new_accessions_count = len(new_accessions)

        progress_bar = st.progress(0)
        percentage_text = st.empty()

        if new_accessions:
            for i, acc in enumerate(new_accessions):
                sequence = fetch_fasta_sequence(acc)
                cursor.execute(
                    f"""INSERT INTO {uc_table} (id, sequence, record_create_ts) 
                        VALUES (?, ?, CURRENT_TIMESTAMP)""",
                    (acc, sequence)
                )
                progress_percentage = (i + 1) / new_accessions_count
                progress_bar.progress(progress_percentage)
                percentage_text.text(f"Updating sequences in database: {int(progress_percentage * 100)}%")

            conn.commit()

        cursor.close()
        conn.close()

        st.success(f"Already existing sequences in UC table {uc_table}: {len(accessions) - new_accessions_count}")
        st.success(f"New sequences added to UC table {uc_table}: {new_accessions_count}")
        st.success(f"Total sequences now in {uc_table}: {len(existing_ids) + new_accessions_count}")

    except Exception as e:
        st.error(f"Error updating UC table {uc_table}: {str(e)}")


def generate_fasta_file(accessions: list):
    """
    Generates a FASTA file from accession numbers stored in Databricks.
    """
    uc_table = "workspace.raw.accession"
    results = {}
    total_accessions = len(accessions)
    try:
        conn = dbh.get_databricks_connection()
        cursor = conn.cursor()
        # Retrieve sequences from UC table
        for acc in accessions:
            cursor.execute(f"SELECT sequence FROM {uc_table} WHERE id = ?", (acc,))
            row = cursor.fetchone()
            if row and row[0]:
                results[acc] = row
        cursor.close()
        conn.close()

        fasta_io = BytesIO()
        for sequence in results.values():
            fasta_io.write(f"{sequence.strip()}\n".encode('utf-8'))

        fasta_io.seek(0)
        return fasta_io, total_accessions, results

    except Exception as e:
        st.error(f"Error generating FASTA file: {str(e)}")
        return None, 0, {}
