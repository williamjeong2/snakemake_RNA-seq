import subprocess
import sys
def install(package):
    subprocess.check_call([sys.executable, "-m", "pip", "install", package])

try:
    import pybiomart, pandas, openpyxl
except:
    install("pybiomart, pandas, openpyxl")
    import pybiomart, pandas

import pandas as pd
from pybiomart import Server

counts = pd.read_table(snakemake.input[0],
                      sep = "\t")
counts['rowSum'] = counts.iloc[:, 6:].sum(axis=1)
counts = counts[counts['rowSum'] != 0].drop(columns="rowSum")

server = Server(host='http://www.ensembl.org')

dataset = (server.marts['ENSEMBL_MART_ENSEMBL']
          .datasets['mmusculus_gene_ensembl'])
bmIDs = dataset.query(attributes=['ensembl_gene_id', 'external_gene_name'])
counts = pd.merge(counts, bmIDs, left_on=['Geneid'], right_on=['Gene stable ID'], how="left")
counts.drop(columns=['Gene stable ID'], inplace=True)
col1 = counts.columns[-1:].to_list()
col2 = counts.columns[1:-1].to_list()
col3 = [counts.columns[0]]
new_col = col3 + col1 + col2
counts = counts[new_col]
counts.drop_duplicates(['Geneid'], keep='first')
counts
counts.to_excel(snakemake.output[0],
               sheet_name="Sheet1",
               index=False)

cmd = "rm -f " + snakemake.input[0]
os.system(cmd)