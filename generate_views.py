import databricks_handler as dbh
import pandas as pd


# Generate view
def generate_view(curated_view, accessions: list):
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
        print(f"Failed to generate view {curated_view.split('.')[-1]}: {e}")
        return pd.DataFrame()  # Return empty DataFrame on error
