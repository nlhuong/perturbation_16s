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
# abundance_table.py
# Created 04/12/2018
# Lan Huong Nguyen, nlhuong90@gmail.com
#
# This program merges abundance  matrices from multiple batches into
# a single table.
#
##########################################################################

# imports
import argparse, os, glob
import numpy as np
import pandas as pd
import multiprocessing as mp

def merge_batches(batch_files, outfile, min_feature_sum=None):
    gene_info_cols = ['Organism', 'Function', 'SEED1', 'SEED2', 'SEED3', 'SEED4']
    smp_df = pd.DataFrame()
    gene_df = pd.DataFrame()    
    for bfile in batch_files:
        df = pd.read_csv(bfile)
        cols = [col for col in df if col.startswith('M')] 
        abund = df.copy()[['GeneID'] + cols]
        annot = df.copy().drop(cols, axis = 1)
        gene_df = pd.concat([gene_df, annot]).drop_duplicates('GeneID').reset_index(drop=True)
        if smp_df.empty:
            smp_df = abund
        else:
            smp_df = pd.merge(smp_df, abund, on='GeneID', how='outer') 
    
    smp_df = pd.merge(smp_df, gene_df, on='GeneID', how='left')
    if min_feature_sum is None:
        smp_df.to_csv(outfile, index = False)
        return
    # Sample columns start with "M" for measurement:
    sample_col = [col for col in merged_df if col.startswith('M')] 
    abundance = smp_df[sample_col]
    feat_sum = abundance.sum(axis=1, skipna=True)
    smp_df = smp_df.loc[feat_sum > min_feature_sum]
    smp_df.to_csv(outfile, index = False)
    return
    
###############################################################################
## Run program
###############################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog='merge_batches',
        description='Merge abundance tables from different batches.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('indir', metavar='indir', type=str, nargs=1,
                        help='directory with all batch abundance tables')
    parser.add_argument('-outdir', dest='outdir', type=str, default='.',
                        help='path to output directory')
    parser.add_argument('-thresh', dest='thresh', type=float, default=None,
                        help='minimum feature sum for filtering')
    args = parser.parse_args()
    print(args)
   
    suffixes = [ db + meas \
                for db in ["refseq", "nr_contigs", "nr_unassembled", "seed"] \
                for meas in ["_raw.csv", "_bp_aligned.csv", "_len_ratio.csv"]]
 
    pool = mp.Pool(processes=len(suffixes))   
    for sff in suffixes:
        batch_files = glob.glob(args.indir[0] + "/*_" + sff )
        outfile = os.path.join(args.outdir, "abund_" + sff)
        pool.apply_async(merge_batches(batch_files, outfile, args.thresh))
    pool.close()
    pool.join()

