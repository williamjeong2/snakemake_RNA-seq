import sys
import pandas as pd

counts = pd.read_table(snakemake.input[0],
                      sep = "\t",
                      skiprows = 1)
counts['rowSum'] = counts.iloc[:, 6:].sum(axis=1)
counts = counts[counts['rowSum'] != 0].drop(columns="rowSum")

counts.to_csv(snakemake.output[0],
              sep = "\t",
              index = False)