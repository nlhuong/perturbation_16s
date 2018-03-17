#! /usr/bin/env python

"""

Compute Workflow Summaries

These functions look at the intermediate output from the metatranscriptomic
workflow, and computes statistics describing the quality of the output at each
stage.

author: krissankaran@stanford.edu
date: 03/15/2018
"""

import pandas as pd
import glob
import os
import subprocess
from collections import OrderedDict
from re import split

###############################################################################
## Small utility functions used later
###############################################################################

def file_len(fname):
    """
    https://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
    """
    p = subprocess.Popen(
        ['wc', '-l', fname],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return int(result.strip().split()[0])


def tail(f, n):
    """
    https://stackoverflow.com/questions/136168/get-last-n-lines-of-a-file-with-python-similar-to-tail
    """
    p = subprocess.Popen(
        ["tail", "-n", str(n), f],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return result


def bam_mapped(f):
    """
    Number of mapped reads in a BAM file

    samtools must be loaded!

    https://www.biostars.org/p/138116/
    """
    p = subprocess.Popen(
        ["samtools view -F 0x40" + f + "| cut -f1 | sort | uniq | wc -l"],
        shell=True,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE
    )
    result, err = p.communicate()
    if p.returncode != 0:
        raise IOError(err)
    return float(result)


###############################################################################
## Functions to calculate summary statistics
###############################################################################

def read_filter_stats(output_dir, subject_id, sample_id):
    """
    Statistics about Read Filtering

    Example:
        output_dir = "/scratch/PI/sph/resilience/metatranscriptomics/processed"
        subject_id = "DBUr_Sub"
        sample_id = "M3303_DBUsw_2r_TrM31"
        read_filter_stats(output_dir, subject_id, sample_id)
    """
    stats = OrderedDict()
    subject_dir = os.path.join(output_dir, subject_id)

    def make_path(dir, suffix):
        return os.path.join(subject_dir, dir, sample_id + suffix)

    ## statistics about trimming / quality filtering
    trim_file = make_path("trimmed", "_trim.fq")
    qual_file = make_path("main", "_qual.fq")
    stats["n_trim"] = file_len(trim_file) / 4.0
    stats["n_qual"] = file_len(qual_file) / 4.0

    ## statistics about vector and host read removal
    unique = make_path("main", "_unique.fq")
    vector_bwa = make_path("aligned", "_univec_bwa.fq")
    vector_blat = make_path("aligned", "_univec_blat.fq")
    host_bwa = make_path("aligned", "_human_bwa.fq")
    host_blat = make_path("aligned", "_human_blat.fq")

    stats["unique"] = file_len(unique) / 4.0
    stats["vector_bwa"] = file_len(vector_bwa) / 4.0
    stats["vector_blat"] = file_len(vector_blat) / 4.0
    stats["human_bwa"] = file_len(vector_bwa) / 4.0
    stats["human_blat"] = file_len(vector_bwa) / 4.0

    ## statistics about rRNA filtering
    mRNA = make_path("main", "_unique_mRNA.fq")
    stats["mRNA"] = file_len(mRNA) / 4.0
    return stats


def annotation_stats(output_dir, subject_id, sample_id):
    """
     Assess the quality of read annotation

    Example:
        output_dir = "/scratch/PI/sph/resilience/metatranscriptomics/processed"
        subject_id = "DBUr_Sub"
        sample_id = "M3303_DBUsw_2r_TrM31"
        annotation_stats(output_dir, subject_id, sample_id)
    """
    stats = OrderedDict()
    subject_dir = os.path.join(output_dir, subject_id)

    def make_path(dir, suffix):
        return os.path.join(subject_dir, dir, sample_id + suffix)

    genus_file = make_path("taxonomy", "_genus_class_summary.txt")
    genus_summary = split("\n|\t", tail(genus_file, 3))

    stats["genus_unassigned_reads"] = float(genus_summary[1])
    stats["genus_unclassified_reads"] = float(genus_summary[5])

    dmnd_refseq = make_path("diamond", "_refseq.dmdout")
    dmnd_seed = make_path("diamond", "_seed.dmdout")
    dmnd_nr_contigs = make_path("diamond", "_nr_contigs.dmdout")
    dmnd_nr_unassembled = make_path("diamond", "_nr_unassembled.dmdout")

    bwa_mcds_contigs_ann = make_path("assembled", "_contigs_annotation_bwa.bam")
    bwa_mcds_unassembled_ann = make_path("assembled", "_unassembled_annotation_bwa.bam")

    stats["dmnd_refseq"] = file_len(dmnd_refseq)
    stats["dmnd_seed"] = file_len(dmnd_seed)
    stats["dmnd_nr_contigs"] = file_len(dmnd_nr_contigs)
    stats["dmnd_nr_unassembled"] = file_len(dmnd_nr_unassembled)
    stats["bwa_mcds_contigs_ann"] = bam_mapped(bwa_mcds_contigs_ann)
    stats["bwa_mcds_unassembled_ann"] = bam_mapped(bwa_mcds_unassembled_ann)

    return stats


def assembly_stats(output_dir, subject_id, sample_id):
    """
    Assess the quality of the assembly

    Example:
        output_dir = "/scratch/PI/sph/resilience/metatranscriptomics/processed"
        subject_id = "DBUr_Sub"
        sample_id = "M3303_DBUsw_2r_TrM31"
        assembly_stats(output_dir, subject_id, sample_id)
    """
    stats = OrderedDict()
    subject_dir = os.path.join(output_dir, subject_id)

    def make_path(dir, suffix):
        return os.path.join(subject_dir, dir, sample_id + suffix)

    original = make_path("main", "_mRNA.fq")
    unassembled = make_path("assembled", "_unassembled.fq")
    contigs = make_path("assembled", "_contigs.fq")
    stats["original"] = file_len(original) / 4.0
    stats["unassembled"] = file_len(unassembled) / 4.0

    return stats


def summary_stats(output_dir):
    """
    Summaries across all samples

    Example:
        output_dir = "/scratch/PI/sph/resilience/metatranscriptomics/processed"
        stats = summary_stats(output_dir)
    """
    stats = dict()
    for subject_id in os.listdir(output_dir):
        sample_ids = glob.glob(os.path.join(output_dir, subject_id, "main", "*_unique.fq"))
        sample_ids = [os.path.basename(x.replace("_unique.fq", "")) for x in sample_ids]
        for sid in sample_ids:
            print("Processing sample " + sid)
            try:
                stats[sid] = OrderedDict()
                stats[sid].update(read_filter_stats(output_dir, subject_id, sid))
                stats[sid].update(assembly_stats(output_dir, subject_id, sid))
                stats[sid].update(annotation_stats(output_dir, subject_id, sid))

            except:
                stats[sid] = "NA"

    stats = pd.DataFrame(stats).T
    stats.index.name = "Meas_ID"
    stats.reset_index(inplace=True)
    return stats
