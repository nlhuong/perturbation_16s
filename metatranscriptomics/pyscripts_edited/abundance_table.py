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
# This program aggregates sample abundance measures into a abundance
# matrix.
#
##########################################################################

# imports
import argparse, os, glob
import numpy as np
import pandas as pd
import multiprocessing as mp 
from functools import reduce

def combine_samples(smp_dict, column, outfile):
    print('Aggregating columns of the samples in a list')    
    smp_lst = []
    for sname in list(smp_dict.keys()):
        res = smp_dict[sname]
        sdf = res[['GeneID', column]]
        sdf.columns = ['GeneID', sname]
        smp_lst.append(sdf)
    
    print('Merging into a data frame from list of samples')
    table_df = reduce(lambda x, y: pd.merge(x, y, on = 'GeneID', how='outer'), smp_lst)
    return table_df

def generate_table(sample_files, outfile):  
    outdir, outfile = os.path.split(outfile)
    outfile = outfile.split(".")[0]
    outfiles = [os.path.join(outdir, outfile + suffix) for suffix in \
        ["_raw.csv", "_bp_aligned.csv", "_len_ratio.csv"]]
    
    abund_meas = ['Raw', 'Length_Aligned', 'Length', 'Length_Ratio',
                  'Log2_Length_Ratio', 'Log2_Normed_Length_Aligned', 
                  "Log2_Normed_Length_Ratio"] 
    #gene_dict = {}
    smp_dict = {}    
    gene_df = pd.DataFrame()
    for smp_file in sample_files:
        indir, smp_name = os.path.split(smp_file)
        smp_name = smp_name.split(".")[0]
        df = pd.read_csv(smp_file)
        smp_dict[smp_name] = df.copy()[['GeneID'] + abund_meas]
        annot = df.copy().drop(abund_meas, axis = 1)
        gene_df = pd.concat([gene_df, annot]).drop_duplicates('GeneID').reset_index(drop=True)
        #annot = annot.set_index('GeneID') 
        #gene_dict.update(annot.T.to_dict())
    
    #gene_df = pd.DataFrame(gene_dict).T
    #gene_df.index.name = 'GeneID'
    #gene_df.reset_index(inplace=True)
    
    pool = mp.Pool(processes=len(outfiles))
    results = {}
    for i, column in enumerate(['Raw', "Length_Aligned", 'Length_Ratio']):
        print("Joining samples for abundance measure: " + column)
        results[outfiles[i]] = pool.apply_async(combine_samples, args=(smp_dict, column, outfiles[i]))
    pool.close()
    pool.join()
   
    for outfile, res in list(results.items()): 
        df = res.get()
        df = pd.merge(df, gene_df, on='GeneID', how='left')
        df.to_csv(outfile, index = False)
    return


###############################################################################
## Run program
###############################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog='abundance_table',
        description='Generate abundance table from sample abundance files.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('procdir', metavar='pdir', type=str, nargs=1,
                        help='directory with all processed files')
    parser.add_argument('-outfile', dest='outfile', type=str,
                        default='abundance.csv',
                        help='path to output file')
    parser.add_argument('-subdir', dest='subdir', type=str, default=None,
                        help='specific subdirectory of processed directory' +
                             ' to limit the summary to')
    #parser.add_argument('-ncores', dest='ncores', type=int, default=10,
    #                    help='number of cores for multithreading')
    args = parser.parse_args()
    print(args)

    for sdir in os.listdir(args.procdir[0]):
        count_dir = os.path.join(args.procdir[0], sdir, "counts")
        if not os.path.isdir(count_dir) or \
            (args.subdir is not None and not args.subdir in sdir):
            continue
        print("MESSAGE: Adding files in " + sdir)
        refseq_dir = os.path.join(count_dir, "dmnd_RefSeq") 
        nr_dir = os.path.join(count_dir, "dmnd_NR") 
        seed_dir = os.path.join(count_dir, "dmnd_SEED") 
        uniref_dir = os.path.join(count_dir, "dmnd_UniRef")
        refseq_files = glob.glob(refseq_dir + "/*") 
        nr_contig_files = glob.glob(nr_dir + "/*_nr_contigs_abund.csv") 
        nr_unassembled_files = glob.glob(nr_dir + "/*_nr_unassembled_abund.csv") 
        seed_files = glob.glob(seed_dir + "/*") 
        uniref_files = glob.glob(uniref_dir + "/*")
           
 
    outfile = args.outfile.split(".")[0]
    if (len(glob.glob(outfile + '_refseq_*')) == 0):
        print('Processing RefSeq aligned abundances...')
        generate_table(refseq_files, outfile + "_refseq") 
    if (len(glob.glob(outfile + '_nr_contigs_*')) == 0):
        print('Processing contigs NR aligned abundances...')
        generate_table(nr_contig_files, outfile + "_nr_contigs") 
    if (len(glob.glob(outfile + '_nr_unassembled_*')) == 0):
        print('Processing unassembled NR aligned abundances...')
        generate_table(nr_unassembled_files, outfile + "_nr_unassembled")
    if (len(glob.glob(outfile + '_seed_*')) == 0):
        print('Processing SEED aligned abundances...')
        generate_table(seed_files, outfile + "_seed")
    if (len(glob.glob(outfile + '_uniref_*')) == 0):
        print('Processing SEED aligned abundances...')
        generate_table(seed_files, outfile + "_uniref")


