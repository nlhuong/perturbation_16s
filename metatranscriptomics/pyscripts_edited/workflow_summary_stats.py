#! /usr/bin/env python

"""

Compute Workflow Summaries

These functions look at the intermediate output from the metatranscriptomic
workflow, and computes statistics describing the quality of the output at each
stage.

author: krissankaran@stanford.edu
date: 03/15/2018
"""

import subprocess
import os
from collections import OrderedDict
from re import split

def file_len(fname):
    """
    https://stackoverflow.com/questions/845058/how-to-get-line-count-cheaply-in-python
    """
    p = subprocess.Popen(['wc', '-l', fname], stdout=subprocess.PIPE,
                                              stderr=subprocess.PIPE)
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
    return result


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
    stats["trim_to_qual"] = 1 - stats["n_qual"] / stats["n_trim"]

    ## statistics about vector and host read removal
    unique = make_path("main", "_unique.fq")
    vector_bwa = make_path("aligned", "_univec_bwa.fq")
    vector_blat = make_path("aligned", "_univec_blat.fq")
    host_bwa = make_path("aligned", "_human_bwa.fq")
    host_blat = make_path("aligned", "_human_blat.fq")

    stats["unique"] = file_len(unique) / 4.0
    stats["vector_bwa"] = file_len(vector_bwa) / 4.0
    stats["vector_blat"] = file_len(vector_blat) / 4.0
    stats["unique_to_vector"] = 1 - stats["vector_blat"] / stats["unique"]
    stats["human_bwa"] = file_len(vector_bwa) / 4.0
    stats["human_blat"] = file_len(vector_bwa) / 4.0
    stats["vector_to_human"] = 1 - stats["human_blat"] / stats["vector_blat"]

    ## statistics about rRNA filtering
    mRNA = make_path("main", "_unique_mRNA.fq")
    stats["mRNA"] = file_len(mRNA) / 4.0
    stats["rRNA_removal"] = 1 - stats["mRNA"] / stats["human_blat"]
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

    stats["genus_unassigned_perc"] = float(genus_summary[0])
    stats["genus_unassigned_reads"] = float(genus_summary[1])
    stats["genus_unclassified_perc"] = float(genus_summary[4])
    stats["genus_unclassified_reads"] = float(genus_summary[5])

    contigs = make_path("assembled", "_contigs.fasta")
    unassembled = make_path("assembled", "_unassembled.fq")
    contigs_unmapped = make_path("assembled", "_contigs_unmapped.fq")
    unassembled_unmapped = make_path("assembled", "_unassembled_unmapped.fq")

    stats["contigs"] = file_len(contigs) / 4.0
    stats["unassembled"] = file_len(unassembled) / 4.0
    stats["contigs_unmapped"] = file_len(contigs_unmapped) / 4.0
    stats["contigs_unmapped_perc"] = stats["contigs_unmapped"] / stats["contigs"]
    stats["unassembled_unmapped"] = file_len(unassembled_unmapped) / 4.0
    stats["unassembled_unmapped_perc"] = stats["unassembled_unmapped"] / stats["unassembled"]

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
    stats["unassembled_prop"] = stats["unassembled"] / stats["original"]

    return stats
