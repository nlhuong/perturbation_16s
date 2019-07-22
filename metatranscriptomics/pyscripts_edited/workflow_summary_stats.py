#! /usr/bin/env python

"""

Compute Workflow Summaries

These functions look at the intermediate output from the metatranscriptomic
workflow, and computes statistics describing the quality of the output at each
stage.

author: krissankaran@stanford.edu, lanhuong@stanford.edu
date: 03/15/2018
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

def fasta_count(fname):
    """
    Only used for counting reads in .fasta files.
    """
    p = subprocess.Popen(
        ["grep -c '>' " +  fname],
        shell=True,
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

def kaiju_count(fname):
    """
    The first column of kaiju output (files _tax_class.tsv)
    are U/C letters for unclassified/classified respectively.
    We count the number of classified reads.
    """
    p = subprocess.Popen(
        ["cut -f1 " +  fname + " | grep -c 'C'"],
        shell=True,
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

def tail(fname, n):
    """
    https://stackoverflow.com/questions/136168/get-last-n-lines-of-a-file-with-python-similar-to-tail
    """
    p = subprocess.Popen(
        ["tail", "-n", str(n), fname],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return result

def contig_reads(fname):
    num_reads_in_contigs = 0
    contig_file = open (fname, "r")
    for line in contig_file:
        splitline = line.split("\t")
        num_reads_in_contigs += int(splitline[1])
    return(num_reads_in_contigs)

def bam_mapped(fname):
    """
    Number of mapped reads in a BAM file

    samtools must be loaded!

    https://www.biostars.org/p/138116/
    """
    p = subprocess.Popen(
        ["samtools view -F 0x40" + fname + "| cut -f1 | sort | uniq | wc -l"],
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    try:
       res = int(float(result))
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
        raw_dir = "/scratch/PI/sph/resilience/metatranscriptomics/raw"
        sub_dir = "Relman_RNAseq_16"
        sample_id = "M3303_DBUsw_2r"
        raw_input_stats(raw_dir, sub_dir, sample_id)
    """
    stats = OrderedDict()
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
    stats["input_fwd"] = file_len(fwd_reads) / 4
    stats["input_rev"] = file_len(fwd_reads) / 4
    return stats


def read_filter_stats(processed_dir, sub_dir, sample_id):
    """
    Statistics about Read Filtering

    Example:
        processed_dir = "/scratch/PI/sph/resilience/metatranscriptomics/processed"
        sub_dir = "Relman_RNAseq_16"
        sample_id = "M3303_DBUsw_2r"
        read_filter_stats(processed_dir, sub_dir, sample_id)
    """
    stats = OrderedDict()
    sub_dir = os.path.join(processed_dir, sub_dir)

    def make_path(dir, suffix):
        return os.path.join(sub_dir, dir, sample_id + suffix)

    ## statistics about trimming / quality filtering
    trim_file = make_path("trimmed", "_trim.fq")
    qual_file = make_path("trimmed", "_qual.fq")
    stats["trimmed"] = file_len(trim_file) / 4
    stats["qual_fltr"] = file_len(qual_file) / 4

    ## statistics about vector and host read removal
    unique = make_path("unique", "_unique.fq")
    vector_bwa = make_path("unique", "_univec_bwa.fq")
    vector_blat = make_path("unique", "_univec_blat.fq")
    host_bwa = make_path("unique", "_human_bwa.fq")
    host_blat = make_path("unique", "_human_blat.fq")
    stats["unique"] = file_len(unique) / 4
    stats["vector_bwa"] = file_len(vector_bwa) / 4
    stats["vector_blat"] = file_len(vector_blat) / 4
    stats["human_bwa"] = file_len(vector_bwa) / 4
    stats["human_blat"] = file_len(vector_bwa) / 4

    ## statistics about rRNA filtering
    mRNA_unique = make_path("unique", "_unique_mRNA.fq")
    mRNA = make_path("main", "_mRNA.fq")
    stats["mRNA_unique"] = file_len(mRNA_unique) / 4
    stats["mRNA"] = file_len(mRNA) / 4
    return stats

def assembly_stats(processed_dir, sub_dir, sample_id):
    """
    Assess the quality of the assembly

    Example:
        processed_dir = "/scratch/PI/sph/resilience/metatranscriptomics/processed"
        subj_dir = "Relman_RNAseq_16"
        sample_id = "M3303_DBUsw_2r"
        assembly_stats(processed_dir, sub_dir, sample_id)
    """
    stats = OrderedDict()
    sub_dir = os.path.join(processed_dir, sub_dir)

    def make_path(dir, suffix):
        return os.path.join(sub_dir, dir, sample_id + suffix)
    
    contig_map_file = make_path("assembled", "_contigs_map.tsv")
    unassembled_file = make_path("assembled", "_unassembled.fq")
    contig_file = make_path("assembled", "_contigs.fasta")
    stats["contigs"] = fasta_count(contig_file)
    stats["contig_reads"] = contig_reads(contig_map_file)
    stats["unassembled_reads"] = file_len(unassembled_file) / 4
    return stats

def annotation_stats(processed_dir, sub_dir, sample_id):
    """
     Assess the quality of read annotation

    Example:
        processed_dir = "/scratch/PI/sph/resilience/metatranscriptomics/processed"
        sub_dir = "Relman_RNAseq_16"
        sample_id = "M3303_DBUsw_2r"
        annotation_stats(processed_dir, sub_dir, sample_id)
    """
    stats = OrderedDict()
    sub_dir = os.path.join(processed_dir, sub_dir)

    def make_path(dir, suffix):
        return os.path.join(sub_dir, dir, sample_id + suffix)

    kaiju_tax_file = make_path("taxonomy", "_tax_class.tsv")
    kaiju_genus_file = make_path("taxonomy", "_genus_class.tsv")

    bwa_mcds_contigs_ann = make_path("assembled", "_contigs_bwa_aligned.fq")
    bwa_mcds_unassembled_ann = make_path("assembled", "_unassembled_bwa_aligned.fq")

    dmnd_refseq = make_path("diamond", "_refseq.dmdout")
    dmnd_seed = make_path("diamond", "_seed.dmdout")
    dmnd_nr_contigs = make_path("diamond", "_nr_contigs.dmdout")
    dmnd_nr_unassembled = make_path("diamond", "_nr_unassembled.dmdout")

    stats["kaiju_tax"] = kaiju_count(kaiju_tax_file)
    stats["kaiju_genus"] = kaiju_count(kaiju_genus_file)

    stats["bwa_mcds_contigs"] = file_len(bwa_mcds_contigs_ann) / 4
    stats["bwa_mcds_unassembled"] = file_len(bwa_mcds_unassembled_ann) / 4

    stats["dmnd_refseq"] = file_len(dmnd_refseq)
    stats["dmnd_seed"] = file_len(dmnd_seed)
    stats["dmnd_nr_contigs"] = file_len(dmnd_nr_contigs)
    stats["dmnd_nr_unassembled"] = file_len(dmnd_nr_unassembled)
    return stats


def process_sample(raw_dir, processed_dir, sub_dir, sample_id):
    print("Processing sample " + sample_id)
    try:
        smp_stats = OrderedDict()
        smp_stats['batch'] = sub_dir
        smp_stats.update(raw_input_stats(raw_dir, sub_dir, sample_id))
        smp_stats.update(read_filter_stats(processed_dir, sub_dir, sample_id))
        smp_stats.update(assembly_stats(processed_dir, sub_dir, sample_id))
        smp_stats.update(annotation_stats(processed_dir, sub_dir, sample_id))
    except:
        smp_stats = "NA"
    return smp_stats


def summary_stats(raw_dir, processed_dir, outfile=None, sub_dir=None, num_cores=1):
    """
    Summaries across all samples

    Example:
        raw_dir = "/scratch/PI/sph/resilience/metatranscriptomics/raw"
        processed_dir = "/scratch/PI/sph/resilience/metatranscriptomics/processed"
        stats = summary_stats(raw_dir, processed_dir)
    """
    num_cores = min(num_cores, mp.cpu_count())
    stats = dict()
    for sdir in os.listdir(processed_dir):
        main_proc_files = os.path.join(processed_dir, sdir, "main")
        if not os.path.isdir(main_proc_files) or \
            (sub_dir is not None and not sub_dir in main_proc_files):
            print("MESSAGE: Skipping " + sdir + " subdirectory.")
            continue
        print("STARTING STATS SUMMARY FOR " + sdir + " subdirectory...")
        sample_ids = glob.glob(os.path.join(main_proc_files,"*_mRNA.fq"))
        sample_ids = [os.path.basename(x.replace("_mRNA.fq", "")) for x in sample_ids]
        print(sample_ids)
        pool = mp.Pool(processes=num_cores)
        results = {}
        for sid in sample_ids:
            results[sid] = pool.apply_async(process_sample, args=(raw_dir, processed_dir, sdir, sid))
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
    parser.add_argument('procdir', metavar='pdir', type=str, nargs=1,
                        help='directory with all processed files')
    parser.add_argument('-outfile', dest='outfile', type=str,
                        default='metat_pipeline_stats.csv',
                        help='path to output file')
    parser.add_argument('-subdir', dest='subdir', type=str, default=None,
                        help='specific subdirectory of processed directory' +
                             ' to limit the summary to')
    parser.add_argument('-ncores', dest='ncores', type=int, default=10,
                        help='number of cores for multithreading')
    args = parser.parse_args()
    print(args)

    stats = summary_stats(args.rawdir[0], args.procdir[0],
                          args.outfile, args.subdir, args.ncores)
