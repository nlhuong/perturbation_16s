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
# diamond_read_counter.py
# Created 03/22/2018
# Lan Huong Nguyen, nlhuong90@gmail.com
# Based on DIAMOND_analysis_counter.py by:
# Sam Westreich, stwestreich@ucdavis.edu, github.com/transcript
#
# This program parses through the results file from a DIAMOND annotation run
# (in BLAST m8 format) to get the results into something more compressed
# and readable.
#
##########################################################################

# imports
import argparse, operator, os,  sys, time, re

parser = argparse.ArgumentParser(description='Aggregate and count reads aligned with DIAMOND.')
parser.add_argument('infile', metavar='I', type=str, nargs=1,
                    help='path to input DIAMOND alignment file in m8 format')
parser.add_argument('dbfile', metavar='DB', type=str, nargs=1,
                    help='path to reference database as .fasta file')
parser.add_argument('-outfile', dest='outfile', type=str, default=None, 
                    help='path to output file. If not specified, the same as ' +
			  'the infile with a replaced suffix')
parser.add_argument('-contig-map', dest='contig_map', type=str,
                    help='contig mapping (.tsv) file if input is contig alignment.' +
                          'The first column is contig id, second is the number of corresponding' +
                          'reads. Counts are then adjusted for number of read in a contig')
parser.add_argument('--fun', dest='fun', action="store_true", default=False, 
                    help='return function names. Default True.')
parser.add_argument('--org', dest='org', action="store_true", 
                    default=False, help='return organism names. Default False')
parser.add_argument('--condense', action="store_true", default=False, 
                    help='condense results to identical function/organism names')


def condense_hit_db(hits_dict, db_dict):
    condensed_hits = {}
    for entry_id in hits_dict.keys():
        try:
	        entry_name = db_dict[entry_id]
	        if entry_name in condensed_hits.keys():
	            condensed_hits[entry_name] += hits_dict[entry_id]
	        else:
		        condensed_hits[entry_name] = hits_dict[entry_id]
        except KeyError:
	        print("KeyError:\t" + entry_id)
	        continue
    return condensed_hits

def save_to_file(outfile_name, hits_dict, db_dict=None, idx = 0, condense = False):
    if condense:
	hits_dict = condense_hit_db(hits_dict, db_dict)
    total_reads = sum(hits_dict.values())
    outfile = open(outfile_name, "w")
    # writing the output
    error_counter = 0
    for k, v in sorted(hits_dict.items(), key=lambda (k,v): -v):
	try:
            q = v * 100 / float(total_reads)
            if condense or db_dict is None:
	        outfile.write (str(q) + "\t" + str(v) + "\t" + k + "\n")
	    else:
		outfile.write (str(q) + "\t" + str(v) + "\t" + k + "\t" + 
                               db_dict[k][idx] + "\n")
	except KeyError:
		outfile.write (str(q) + "\t" + str(v) + "\tWARNING: Key not found for " + k + "\n")
		error_counter += 1
		continue
    print("\nAnnotations saved to file: '" + outfile_name + "'.")
    print("Number of errors: " + str(error_counter))
    outfile.close()
    return 


args = parser.parse_args()
print(args)

if args.contig_map is not None:
    contig_file = open(args.contig_map, "r")
    contig_map = {}
    for line in contig_file:
        splitline = line.split("\t")
        contig_map[splitline[0]] = int(splitline[1])


infile = open(args.infile[0], "r")
# setting up databases
RefSeq_hit_count_db = {}
unique_seq_db = {}
line_counter = 0

# reading through the infile - the DIAMOND results m8 format
print("\nNow reading through the m8 results infile.")
t0 = time.clock()
for line in infile:
    line_counter += 1
    splitline = line.split("\t")
    if line_counter % 1000000 == 0:
        t1 = time.clock()
	print(str(line_counter)[:-6] + "M lines processed so far in " + 
              str(t1-t0) + " seconds.")

    unique_seq_db[splitline[0]] = 1
    if args.contig_map is not None:
	increment = contig_map[splitline[0]]
    else:
        increment = 1
    try:
        RefSeq_hit_count_db[splitline[1]] += increment 
    except KeyError:
        RefSeq_hit_count_db[splitline[1]] = increment
	continue

t2 = time.clock()
print("\nAnalysis of " + args.infile[0] + " complete.")
print("Number of total lines: " + str(line_counter))
print("Number of unique sequences: " + str(len(unique_seq_db)))
print("Time elapsed: " + str(t2-t0) + " seconds.")

infile.close()
if args.outfile is None:
    outdir = '.'
    outfile = args.infile[0].split(".")[0]
else:
    outdir, outfile = os.path.split(args.outfile)
    outfile = outfile.split(".")[0]

save_to_file(os.path.join(outdir, outfile + "_genes.tsv"), 
             RefSeq_hit_count_db)

if not any([args.fun, args.org]):
    sys.exit(0)

# time to search for these in the reference database
print("\nStarting database analysis now.")
# optional outfile of specific organism results
db_gene_info = {}

db = open(args.dbfile[0], "r")    
t0 = time.clock()
db_line_counter = 0
for line in db:
    if line.startswith(">"):
	    db_line_counter += 1
    else:
	    continue
	# id, organism and function names [https://stackoverflow.com/questions/6109882/regex-match-all-characters-between-two-strings]
    db_id = re.search("(?<=>)[^ ]+", line)
    db_id = db_id.group()
    db_entry = re.search("(?<= )(.*)(?= \[)", line)
    if db_entry is None: # some lines missing a space
	    db_entry = re.search("(?<= )(.*)(?=\[)", line)
    if db_entry is None: # The NR database has entries without "[", "]"
	    db_entry = line.strip('>' + db_id + ' ')
    else:
	    db_entry = db_entry.group()

    db_org = re.search("(?<=\[)(.*)(?=\])", line)
    if db_org is None: # In case no labels inside "[Genus Species]"
	    db_org = ''
    else:
        db_org = db_org.group()

    db_gene_info[db_id] = [db_org, db_entry]

	# line counter to show progress
    if db_line_counter % 1000000 == 0:     # each million
	    t1 = time.clock()
	    print(str(db_line_counter)[:-6] + "M lines processed so far in " + \
		      str(t1-t0) + " seconds.")

t2 = time.clock()
print("\nSuccess!")
print("Time elapsed: " + str(t2-t0) + " seconds.")
print("Number of lines: " + str(db_line_counter))

db.close()


# creating the outfiles
if args.fun:
    out = os.path.join(outdir, outfile + "_function.tsv")
    save_to_file(out, RefSeq_hit_count_db, db_gene_info, 1, args.condense)
if args.org:
    out = os.path.join(outdir, outfile + "_organism.tsv")
    save_to_file(out, RefSeq_hit_count_db, db_gene_info, 0, args.condense)


