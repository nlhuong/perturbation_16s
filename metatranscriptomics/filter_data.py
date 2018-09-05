import pandas as pd

minSum = 5
chunksize = 10 ** 5

refseq =  pd.DataFrame()
for chunk in pd.read_csv("abund_refseq_len_ratio.csv", chunksize=chunksize, low_memory= False):
    print(chunk.head())
    rowsums_chunk = chunk.drop(["GeneID", "Organism", "Function ", "Length"], axis = 1)
    rowsums_chunk = rowsums_chunk.sum(axis=1)
    chunk_fltr = chunk.loc[rowsums_chunk > minSum]
    refseq = pd.concat([refseq, chunk_fltr], axis = 0)
    refseq.reset_index(inplace=True)
refseq.to_csv("abund_refseq_len_ratio_fltr.csv")

seed =  pd.DataFrame()
for chunk in pd.read_csv("abund_seed_len_ratio.csv", chunksize=chunksize, low_memory= False):
    rowsums_chunk = chunk.drop(["GeneID", "Organism", "Function", "Length"], axis = 1)
    rowsums_chunk = rowsums_chunk.sum(axis=1)
    chunk_fltr = chunk.loc[rowsums_chunk > minSum]
    refseq = pd.concat([refseq, chunk_fltr], axis = 0)
    refseq.reset_index(inplace=True)

seed.to_csv("abund_seed_len_ratio_fltr.csv")

