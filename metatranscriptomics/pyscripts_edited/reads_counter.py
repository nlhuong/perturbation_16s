#!/usr/bin/env Python
##########################################################################
#
# Copyright (C) 2018
#
# This program is free software; you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by the
# Free Software Foundation;
#
# This program is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
# See the GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
##########################################################################
#
# reads_counter.py
# Created 03/22/2018
# Lan Huong Nguyen, nlhuong90@gmail.com
#
# This program parses through the standard BLAST (m8) output files from 
# DIAMOND or other aligners to compute the number of reads and the total
# length of bp aligned to each reference sequence (gene).
#
##########################################################################

# imports
import argparse, os
import numpy as np
import pandas as pd


def sample_reads(infile, dbfile, outfile):
    colnames = ['qseqid', 'GeneID', 'pident', 'Length_Aligned',
                'mismatch', 'gapopen', 'qstart', 'qend', 'sstart', 'send', 
                'evalue', 'bitscore']

    infile_kwargs = {'sep': '\t', 'chunksize': 100000, 'low_memory': False, \
                     'header': None, 'names': colnames}

    gene_abund = pd.DataFrame()
    for chunk in pd.read_csv(infile, **infile_kwargs):
        chunk = chunk[['GeneID', 'Length_Aligned']]
        chunk['Raw'] = 1
        gene_abund = pd.concat([gene_abund, chunk], axis = 0)
        gene_abund = gene_abund.groupby('GeneID')[['Length_Aligned', 'Raw']].sum()
    
    gene_abund.index.name = 'GeneID'
    gene_abund.reset_index(inplace=True)
    
    db_map = pd.DataFrame()
    for chunk in pd.read_csv(dbfile, sep='\t', chunksize = 100000, 
                             low_memory = False):
        chunk_fltr = chunk[chunk['GeneID'].isin(gene_abund['GeneID'])]
        db_map = pd.concat([db_map, chunk_fltr], axis = 0)

    gene_abund = pd.merge(gene_abund, db_map[['GeneID', 'Length']],
                          on='GeneID', how='left')
    
    sum_bp_aligned = np.sum(gene_abund['Length_Aligned'])
    gene_abund['Length_Ratio'] = \
        gene_abund['Length_Aligned'] / gene_abund['Length']
    gene_abund['Log2_Length_Ratio'] = \
        np.log2(gene_abund['Length_Aligned']) - np.log2(gene_abund['Length'])
    gene_abund['Log2_Normed_Length_Aligned'] = \
        np.log2(gene_abund['Length_Aligned']) - np.log2(sum_bp_aligned)
    gene_abund['Log2_Normed_Length_Ratio'] = \
        gene_abund['Log2_Normed_Length_Aligned'] - np.log2(gene_abund['Length'])
    
    gene_abund = pd.merge(gene_abund, db_map.drop(['Length'], axis = 1),
                          on='GeneID', how='left')
    gene_abund.to_csv(outfile, index=False)
    



###############################################################################
## Run program
###############################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog='reads_counter',
        description='Returns read counts and total number if bp aligned ' + 
                    'to each reference sequence using standard BLAST output',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
        
    parser.add_argument('infile', metavar='infile', type=str, nargs=1,
                        help='input file - BLAST format output file from an aligner')
    parser.add_argument('dbfile', metavar='dbfile', type=str, nargs=1,
                        help='reference sequence lengths and info tsv file')
    parser.add_argument('-outfile', dest='outfile', type=str, default=None,
                        help='output file. If not specified, the same as '+
                             'infile with "_abund.csv" suffix')
    args = parser.parse_args()

    if args.outfile is None:
        outdir = '.'
        outfile = args.infile[0].split(".")[0]
        outfile = os.path.join(outdir, outfile + "_abund.csv")
    else:
        outfile = args.outfile
   
    if (not '.tsv' in args.dbfile[0]):
        print("ERROR: dbfile must be a .tsv file with info on gene lengths")
        exit(1) 
    
    print("Saving to: " + outfile)
    sample_reads(args.infile[0], args.dbfile[0], outfile)

