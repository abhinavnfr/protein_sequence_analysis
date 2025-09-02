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


# def get_curated_blast_sequence():
#     uc_table = "workspace.curated.blast_sequence"
#     try:
#         conn = dbh.get_databricks_connection()
#         cursor = conn.cursor()

#         cursor.execute(f"SELECT fasta_sequence FROM {uc_table} WHERE id IN ({})")
#         existing_ids = set(row[0] for row in cursor.fetchall())
