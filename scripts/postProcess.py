import subprocess
import sys
import os
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
                      sep = "\t",
                      skiprows = 1)
counts['rowSum'] = counts.iloc[:, 6:].sum(axis=1)
counts = counts[counts['rowSum'] != 0].drop(columns="rowSum")
counts.drop(counts.columns[1:6], axis=1, inplace=True)
server = Server(host='http://www.ensembl.org')

if snakemake.params[0].find("HOMO") >= 0 or snakemake.params[0].find("HUMAN") >= 0:
    dataset = (server.marts['ENSEMBL_MART_ENSEMBL'].datasets['hsapiens_gene_ensembl'])
if snakemake.params[0].find("MUS") >= 0 or snakemake.params[0].find("MOUSE") >= 0:
    dataset = (server.marts['ENSEMBL_MART_ENSEMBL'].datasets['mmusculus_gene_ensembl'])

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
counts.to_excel(snakemake.output[0].split(".")[0]+".xlsx",
               sheet_name="Sheet1",
               index=False)
counts.to_csv(snakemake.output[0],
        sep="\t",
        index=False)

cmd = "rm -f " + snakemake.input[0]
os.system(cmd)
