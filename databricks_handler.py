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




