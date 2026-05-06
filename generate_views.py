import databricks_handler as dbh
import pandas as pd
import streamlit as st


# Generate view
def generate_view(curated_view, accessions: list):
    try:
        conn = dbh.get_databricks_connection()
        cursor = conn.cursor()

        placeholders = ""
        for acc in accessions:
            placeholders += f"'{acc}',"
        placeholders = placeholders[:-1]

        if curated_view == "workspace.curated.blast_sequence":
            sql_query = f"""SELECT DISTINCT * FROM {curated_view}
                            WHERE TRIM(accession_number) IN ({placeholders})
                        """
        else:
            sql_query = f"""SELECT DISTINCT * FROM {curated_view}
                            WHERE TRIM(accession_number) IN ({placeholders}) OR TRIM(blast_of_id) IN ({placeholders})
                        """
        
        cursor.execute(sql_query)
        rows = cursor.fetchall()
        columns = [desc[0] for desc in cursor.description] # Extract column names from cursor description
        df = pd.DataFrame(rows, columns=columns) # Return result as pandas DataFrame
        
        cursor.close()
        conn.close()
        
        return df
        
    except Exception as e:
        st.write(f"Failed to generate view {curated_view.split('.')[-1]}: {e}")
        return pd.DataFrame()  # Return empty DataFrame on error
