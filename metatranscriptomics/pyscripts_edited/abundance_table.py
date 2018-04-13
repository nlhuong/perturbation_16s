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
import argparse, os
import numpy as np
import pandas as pd
import multiprocessing as mp 
from functools import reduce

def generate_table(sample_files, outfile, num_cores=10):
    pool = mp.Pool(processes=num_cores)
    results = {}
         
    for smp_file in sample_files:
        indir, smp_name = os.path.split(smp_file)
        results[smp_name] = pool.apply_async(pd.read_csv, args=(smp_file))
    
    pool.close()
    pool.join()
    results = {sname: res.get() for sname, res in results.items()}
    
    raw = []
    bp_aligned = []
    log_coverage = []
    
    for sname in list(results.keys()):
        res = results[sname]
        sraw = res[['GeneID', 'Raw']]
        sbp = res[['GeneID', "Length_Aligned"]]
        slogcov = res[['GeneID', 'Log2_Normed_Length_Aligned']]
        sraw.columns = ['GeneID', sname]
        sbp.columns = ['GeneID', sname]
        slogcov.columns = ['GeneID', sname]
        
        raw.append(sraw)
        bp_aligned.append(sbp)
        log_coverage.append(slogcov)
        
    raw_df = reduce(
        lambda left, right: pd.merge(left, right, on='GeneID', how='outer'), 
        raw)
    bp_aligned_df = reduce(
        lambda left, right: pd.merge(left, right, on='GeneID', how='outer'), 
        bp_aligned)
    log_coverage_df = reduce(
        lambda left, right: pd.merge(left, right, on='GeneID', how='outer'), 
        log_coverage)

    outdir, outfile = os.path.split(outfile)
    outfile = outfile.split(".")[0]
    raw_df.to_csv(os.path.join(outdir, outfile + "_raw.csv"), index = False)
    bp_aligned_df.to_csv(os.path.join(outdir, outfile + "_bp_aligned.csv"), index = False)
    log_coverage_df.to_csv(os.path.join(outdir, outfile + "_log2_cov.csv"), index = False)
    return


###############################################################################
## Run program
###############################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog='abundance_table',
        description='Generate anundance table from sample abundance files.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('procdir', metavar='pdir', type=str, nargs=1,
                        help='directory with all processed files')
    parser.add_argument('-outfile', dest='outfile', type=str,
                        default='abundance.csv',
                        help='path to output file')
    parser.add_argument('-subdir', dest='subdir', type=str, default=None,
                        help='specific subdirectory of processed directory' +
                             ' to limit the summary to')
    parser.add_argument('-ncores', dest='ncores', type=int, default=10,
                        help='number of cores for multithreading')
    args = parser.parse_args()
    print(args)

    refseq_files = []
    nr_files = []
    seed_files = []
    
    for sdir in os.listdir(processed_dir):
        if (args.subdir is not None and not args.subdir in sdir):
            continue
        refseq_dir = os.path.join(args.procdir[0], sdir, "counts", "dmnd_RefSeq") 
        nr_dir = os.path.join(args.procdir[0], sdir, "counts", "dmnd_NR") 
        seed_dir = os.path.join(args.procdir[0], sdir, "counts", "dmnd_SEED") 
        refseq_files.append(os.listdir(refseq_dir)) 
        nr_files.append(os.listdir(nr_dir)) 
        seed_files.append(os.listdir(seed_dir)) 

    generate_table(refseq_files, args.outfile, args.ncores)
    generate_table(nr_files, args.outfile, args.ncores)
    generate_table(seed_files, args.outfile, args.ncores)


