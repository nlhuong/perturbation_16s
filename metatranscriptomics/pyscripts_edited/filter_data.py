import pandas as pd
import feather
import os.path

minSum = 5
chunksize = 10 ** 5

#refseq =  pd.DataFrame()
infile_names = ["abund_refseq_len_ratio.csv", "abund_seed_len_ratio.csv",
                "abund_uniref_len_ratio.csv"]
outfile_names = [infile.replace(".csv", "_fltr.csv") for infile in infile_names]

for i in range(len(infile_names)):
    tab  = pd.DataFrame()
    if os.path.isfile(outfile_names[i].replace(".csv", ".feather")) :
        continue
    for chunk in pd.read_csv(infile_names[i], chunksize=chunksize, low_memory= False):
        rowsums_chunk = chunk.drop(["GeneID", "Organism", "Function ", "Length"], axis = 1)
        rowsums_chunk = rowsums_chunk.sum(axis=1)
        chunk_fltr = chunk.loc[rowsums_chunk > minSum]
        tab = pd.concat([tab, chunk_fltr], axis = 0)
        tab.reset_index()
    tab.to_csv(outfile_names[i])
    feather.write_dataframe(tab, outfile_names[i].replace(".csv", ".feather"))

