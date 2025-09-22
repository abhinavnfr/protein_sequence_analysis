import databricks_handler as dbh
import pandas as pd


# Generate view EffectorP
def generate_view_effectorp(accessions: list):
    curated_view = "workspace.curated.effectorp"
    try:
        conn = dbh.get_databricks_connection()
        cursor = conn.cursor()
        placeholders = ','.join(["?"] * len(accessions))
        sql_query = f"""SELECT * FROM {curated_view}
                        WHERE accession_number IN ({placeholders})
                    """
        cursor.execute(sql_query, accessions)
        rows = cursor.fetchall()
        columns = [desc[0] for desc in cursor.description] # Extract column names from cursor description
        df = pd.DataFrame(rows, columns=columns) # Return result as pandas DataFrame
        
        cursor.close()
        conn.close()
        
        return df
        
    except Exception as e:
        print(f"Failed to generate view EffectorP: {e}")
        return pd.DataFrame()  # Return empty DataFrame on error

