#!/usr/bin/env Python
##########################################################################
#
# Copyright (C) 2015-2016 Sam Westreich
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
# generate_count_matrix.py
# Created 3/12/2018, this version created 3/12/2018
# Lan Huong Nguyen, nlhuong90@gmail.com, github.com/nlhuong
#
# This program aggregates function/gene/organism counts for individual
# samples into a matrix of counts with each column corresponding to 
# a different sample. Inputs are .tsv files generated with
# DIAMOND_analysis_counter.py, a modified version of a script from 
# SAMSA2 framework 
#
# Usage:
#
# -D            input directory         specifies the directory containing 
#                                       input .tsv files with sample counts
#                                       included in the second column and 
#                                       function/gene/organism name in the 
#                                       third column 
# -O            output                  output file path name            
# -S            suffix                  specifies suffix of the subset of
#                                       files of interest  
##########################################################################

import glob, sys, os
import pandas as pd

# String searching function:
def string_find(usage_term):
        for idx, elem in enumerate(sys.argv):
                next_elem = sys.argv[(idx + 1) % len(sys.argv)]
                if elem == usage_term:
                         return next_elem

# If not be specified input directory is the current directory
if "-D" in sys.argv:
        input_dir = string_find("-D")
else:
	input_dir = os.getcwd()
	print("WARNING: using " + input_dir + " as input directory.")
# If not specified use all files in input directory with suffix .dm_function.tsv
if "-S" in sys.argv:
        suffix = string_find("-S")
else:
	suffix = ".dm_function.tsv"
        print("WARNING: aggregating all *" + suffix + " files in the current directory.")

# If not specified output file with name count_matrix.csv in the input directory
if "-O" in sys.argv:
        output = string_find("-O")
else:
	output = os.path.join(input_dir, "count_matrix.csv")
        print("WARNING: output saved in " + output + " file.")

filenames = glob.glob( os.path.join(input_dir, '*' + suffix) )
sample_lst = []
for idx, filename in enumerate(filenames):
    sample = filename.split("/")[-1].strip(suffix)
    sample_lst.append(sample)
    if ".hierarchy" in filename:
        colnames = ['freq', sample] + ['SEED_' + str(i) for i in range(5)]
    else:
        colnames = ['freq', sample, 'feature_name']
    sample_df = pd.read_table(filename, header=None, names=colnames)
    col_intersect = list(set(colnames) - set([sample, 'freq'])) 
    sample_df = sample_df[col_intersect + [sample]]
    #print(str(idx) + " sample: " +  sample)
    if idx == 0:
        count_matrix = sample_df
        continue
    count_matrix = pd.merge(count_matrix, sample_df, on = col_intersect, how = 'outer')

count_matrix = count_matrix.fillna(0)
count_matrix[sample_lst] = count_matrix[sample_lst].astype(int)
count_matrix.to_csv(output)
