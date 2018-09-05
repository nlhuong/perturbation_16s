#! /usr/bin/env python

"""

Compute Workflow Summaries

These functions look at the intermediate output from the metatranscriptomic
workflow, and computes statistics describing the quality of the output at each
stage.

author: lanhuong@stanford.edu
date: 06/09/2018
"""

import argparse, glob, os, subprocess
import pandas as pd
import multiprocessing as mp
from collections import OrderedDict
from re import split

###############################################################################
## Small utility functions used later
###############################################################################

def file_len(fname):
    """
    https://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
    """
    if '.gz' in fname:
        bash_command = ["gunzip -c " + fname + " | wc -l"]
    else:
        bash_command = ['wc -l '+  fname]
    p = subprocess.Popen(
        bash_command,
        shell = True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    try:
        res = int(result.strip().split()[0])
    except:
        res = "NA"
    return res




###############################################################################
## Functions to calculate summary statistics
###############################################################################

def raw_input_stats(raw_dir, sub_dir, sample_id):
    """
    Raw input file read counts.

    Example:
        raw_dir = "/scratch/PI/sph/resilience/metagenomics/raw"
        sub_dir = "Relman_DNA_Pool3.2017310"
        sample_id = "M3499_CAAst_28_5d_GAACAGGCAT_L006"
        raw_input_stats(raw_dir, sub_dir, sample_id)
    """
    print("Counting raw reads for sample " + sample_id)
    sub_dir = os.path.join(raw_dir, sub_dir)
    read_files = [ f for f in os.listdir(sub_dir) \
                   if os.path.isfile(os.path.join(sub_dir, f )) and sample_id in f]
    fwd_reads = [f for f in read_files \
                 if '_R1_001.fastq' in f or '_1P.fq.gz' in f]
    rev_reads = [f for f in read_files \
                 if '_R2_001.fastq' in f or '_2P.fq.gz' in f]
    fwd_reads = os.path.join(sub_dir, fwd_reads[0])
    rev_reads = os.path.join(sub_dir, rev_reads[0])
    ## statistics about trimming / quality filtering
    try:
        stats = OrderedDict()
        stats['batch'] = sub_dir 
        stats["input_fwd"] = file_len(fwd_reads) / 4
        stats["input_rev"] = file_len(fwd_reads) / 4
    except:
        stats = "NA"
    return stats


def summary_stats(raw_dir, merged_abund, outfile=None, sub_dir=None, num_cores=1):
    """
    Summaries across all samples

    Example:
        raw_dir = "/scratch/PI/sph/resilience/metagenomics/raw"
        merged_abund = "/scratch/PI/sph/resilience/metagenomics/merged/count_reads.txt"
        stats = summary_stats(raw_dir, merged_abund)
    """
    num_cores = min(num_cores, mp.cpu_count())
    stats = dict()
    for sdir in os.listdir(raw_dir):
        sdir_path = os.path.join(raw_dir, sdir)
        if sdir == "plate_8.201681" or not os.path.isdir(sdir_path):  
            print("MESSAGE: Skipping " + sdir)
            continue
        print("STARTING STATS SUMMARY FOR " + sdir + " subdirectory...")
        print(sdir_path)
        sample_ids = glob.glob(os.path.join(sdir_path, "*_R1_001.fastq"))
        sample_ids = [os.path.basename(x.replace("_R1_001.fastq", "")) 
                      for x in sample_ids]
        print(sample_ids)
        pool = mp.Pool(processes=num_cores)
        results = {}
        for sid in sample_ids:
            results[sid] = pool.apply_async(raw_input_stats, 
                args=(raw_dir, sdir, sid))
        pool.close()
        pool.join()
        results = {sid: res.get() for sid, res in results.items()}
        print(results[sample_ids[0]])
        stats.update(results)

    column_names = stats[sample_ids[0]].keys()
    stats = pd.DataFrame(stats).T
    stats = stats[column_names]
    stats.index.name = "Meas_ID"
    stats.reset_index(inplace=True)
    stats['Meas_ID'] = [measid.split("_")[0] for measid in stats['Meas_ID']]
    merged = pd.read_csv(merged_abund, sep = "\t")    
    merged_smp_count = merged.drop(['species_id'], axis=1)
    merged_smp_count = merged_smp_count.sum().to_frame()
    merged_smp_count.columns = ['annotated']
    merged_smp_count.index.name = "Meas_ID"
    merged_smp_count.reset_index(inplace=True)

    stats = pd.merge(stats, merged_smp_count, on = "Meas_ID", how = "left")
    if outfile is not None:
        stats.to_csv(outfile, index=False)
    return stats


###############################################################################
## Run program
###############################################################################

if __name__ == "__main__":

    parser = argparse.ArgumentParser(
        prog='summary_stats',
        description='Summarize read coverage througout metatranscriptomics pipeline.',
        formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    parser.add_argument('rawdir', metavar='rdir', type=str, nargs=1,
                        help='directory with raw read fastq files')
    parser.add_argument('mergefile', metavar='mergefile', type=str, nargs=1,
                        help='midas count table file path')
    parser.add_argument('-outfile', dest='outfile', type=str,
                        default='metag_midas_stats.csv',
                        help='path to output file')
    parser.add_argument('-ncores', dest='ncores', type=int, default=10,
                        help='number of cores for multithreading')
    args = parser.parse_args()
    print(args)

    stats = summary_stats(args.rawdir[0], args.mergefile[0],
                          args.outfile, args.ncores)


