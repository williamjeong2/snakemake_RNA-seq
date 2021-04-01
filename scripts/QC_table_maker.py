import os
import pandas as pd

def search(dirname, extension = '.txt'):
    global file_list
    file_list = []
    filenames = os.listdir(dirname)
    for filename in filenames:
        full_filename = os.path.join(dirname, filename)
        ext = os.path.splitext(full_filename)[-1]
        if ext == extension:
            file_list.append(full_filename)
    return file_list
            
file_list = search(snakemake.input[0])
file_list.sort()

res = pd.DataFrame(columns=['Sample Name', 'before_total reads', 'before_total bases',
                                             'before_Q20 bases', 'before_Q30 bases',
                                             'after_total reads', 'after_total bases',
                                             'after_Q20 bases', 'after_Q30 bases'])
index_list = [f.split('.')[1].replace('/','') for f in file_list]
res_index = i = 0

for file_name in file_list:
    df = pd.read_csv(file_name, sep = "\t", skip_blank_lines = True, header = None)
    df = df[df.iloc[:,0].str.startswith(('Read1', 'Read2', 'total', 'Q20', 'Q30'))]
    df.columns = ['key']
    df[["key", "value"]] = df['key'].str.split(':', expand=True)
    df.reset_index(drop = True, inplace = True)
    df.drop(axis = 0, index = [0, 5, 10, 15])
    df = df.reindex([1,2,3,4,11,12,13,14,6,7,8,9,16,17,18,19])
    df.reset_index(drop = True, inplace = True)
    
    value1 = df.value[0:8]
    value2 = df.value[8:]
    
    res.loc[res_index] = [index_list[i] + "_1"] + list(value1)
    res_index += 1
    res.loc[res_index] = [index_list[i] + "_2"] + list(value2)
    res_index += 1
    i += 1

res.to_csv(snakemake.output[0], sep="\t", index = False)

