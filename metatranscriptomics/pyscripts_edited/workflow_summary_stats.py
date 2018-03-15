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


def read_filter_stats(output_dir, subject_id, sample_id):
    """
    Statistics about Read Filtering

    Example:
        output_dir = "/scratch/PI/sph/resilience/metatranscriptomics/processed"
        subject_id = "DBUr_Sub"
        sample_id = "M3303_DBUsw_2r_TrM31"
        read_filter_stats(output_dir, subject_id, sample_id)
    """
    stats = dict()

    ## statistics about trimming / quality filtering
    trim_file = os.path.join(output_dir, subject_id, "trimmed", sample_id + "_trim.fq")
    qual_file = os.path.join(output_dir, subject_id, "main", sample_id + "_qual.fq")
    stats["n_trim"] = file_len(trim_file)
    stats["n_qual"] = file_len(qual_file)
    stats["trim_to_qual"] = 1 - float(stats["n_qual"]) / stats["n_trim"]

    return stats
