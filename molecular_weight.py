from Bio.SeqUtils.ProtParam import ProteinAnalysis
from pyspark.sql.functions import udf
from pyspark.sql.types import DoubleType, IntegerType


# Define UDFs for molecular weight, isoelectric point, and length calculation
def calc_mw_kda(sequence):
    return ProteinAnalysis(sequence).molecular_weight() / 1000

def calc_pi(sequence):
    return ProteinAnalysis(sequence).isoelectric_point()

def calc_length(sequence):
    return len(sequence)

calc_mw_kda_udf = udf(calc_mw_kda, DoubleType())
calc_pi_udf = udf(calc_pi, DoubleType())
calc_length_udf = udf(calc_length, IntegerType())

# Register UDFs
spark.udf.register("calc_mw_kda_udf", calc_mw_kda_udf)
spark.udf.register("calc_pi_udf", calc_pi_udf)
spark.udf.register("calc_length_udf", calc_length_udf)

# Perform update with a SQL query
uc_table = "workspace.raw.protein"

update_sql = f"""
UPDATE {uc_table}
SET molecular_weight_kda = ROUND(calc_mw_kda_udf(regexp_replace(fasta_sequence, '>(.|\\r|\\n)*?\\n', '')), 2),
    isoelectric_point_pi = ROUND(calc_pi_udf(regexp_replace(fasta_sequence, '>(.|\\r|\\n)*?\\n', '')), 2),
    sequence_length = calc_length_udf(regexp_replace(fasta_sequence, '>(.|\\r|\\n)*?\\n', ''))
WHERE molecular_weight_kda IS NULL
AND blast_of_id IS NOT NULL 
"""

spark.sql(update_sql)
