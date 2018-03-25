#! /usr/bin/env python

"""
Wrapper for workflow_summary_stats.py.

This program summarizes the coverage statistics for 
the metatranscriptomics data processing pipeline.
The program outputs a table with rows corresponding
to individual samples, and columns corresponding to
the number of reads retained in each processing
steps from the raw input to the final output.

author: nlhuong90@gmail.com
date: 03/19/2018

Usage:

-R        raw data directory            
-P        processed data directory
-O        output file path
-N        number of threads for parallel processing. Default is 1.
"""

import sys
from workflow_summary_stats import *

# String searching function:
def string_find(usage_term):
    for idx, elem in enumerate(sys.argv):
        next_elem = sys.argv[(idx + 1) % len(sys.argv)]
	if elem == usage_term:
            return next_elem


# checking for argument specification
if "-R" not in sys.argv:
    sys.exit("ERROR: raw data directory must be specified with '-R' flag.")

if "-P" not in sys.argv:
    sys.exit("ERROR: processed data directory must be specified with '-P' flag.")

raw_dir = string_find("-R")
processed_dir = string_find("-P") 

if "-N" in sys.argv:
    ncores = int(string_find("-N"))

if "-O" in sys.argv:
    output = string_find("-O")
else:
    output = "count_matrix.csv"
    print("WARNING: output saved in " + output + " file.")

stats = summary_stats(raw_dir, processed_dir, output, ncores)





