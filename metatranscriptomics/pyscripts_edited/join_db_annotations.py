
import argparse, sys, os
import pandas as pd


parser = argparse.ArgumentParser(description='Add db annotations to count matrix.') 
parser.add_argument('infile', metavar='I', type=str, nargs=1,
                    help='path to input count matrix file')
parser.add_argument('dbfile', metavar='DB', type=str, nargs=1,
                    help='path to database annotations mapping file (.tsv)')
parser.add_argument('-outfile', dest='outfile', type=str, default=None,
                    help='path to output file. If not specified, based ' +
                         'on the input file')

args = parser.parse_args()
print(args)

if args.outfile is None:
    outdir, outfile = os.path.split(args.infile[0])
    annotfile = os.path.join(outdir, "Annotated_" + outfile)
    genelenfile = os.path.join(outdir, "gene_len_" + outfile.strip("abund_")) 

cnt_mat = pd.read_csv(args.infile[0])
print(cnt_mat.head())

db_map = pd.DataFrame()
for chunk in pd.read_csv(args.dbfile[0], sep='\t', chunksize = 100000, low_memory = False):
    chunk_fltr = chunk[chunk['GeneID'].isin(cnt_mat['GeneID'])]
    db_map = pd.concat([db_map, chunk_fltr])
print(db_map.head())
db_map.to_csv(genelenfile, index=False)



#res = pd.merge(cnt_mat, db_map, on='GeneID', how='left')
#print(res.head())
#
#res.to_csv(annotfile, index=False)



