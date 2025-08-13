import streamlit as st
from databricks import sql
from Bio import Entrez


# Set databricks connection from secrets
def get_databricks_connection():
    return sql.connect(
        server_hostname=st.secrets["databricks"]["server_hostname"],
        http_path=st.secrets["databricks"]["http_path"],
        access_token=st.secrets["databricks"]["access_token"]
    )

# Set Entrez email from secrets
def get_entrez_email():
    return st.secrets["entrez"]["email"]


# Update UC table workspace.raw.accession
def update_uc_table_accession(accession_numbers):
    uc_table = "workspace.raw.accession"
    conn = get_databricks_connection()
    cursor = conn.cursor()
    
    # Fetch existing accession IDs
    cursor.execute(f"SELECT id FROM {uc_table}")
    existing_ids = set(row[0] for row in cursor.fetchall())
    
    new_accessions = [acc for acc in accession_numbers if acc not in existing_ids]
    
    # Insert only new ones
    if new_accessions:
        for acc in new_accessions:
            cursor.execute(f"INSERT INTO {uc_table} (id) VALUES (?)", (acc,))
        conn.commit()
    
    cursor.close()
    conn.close()
    return uc_table, existing_ids, new_accessions