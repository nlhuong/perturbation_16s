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
# DIAMOND_analysis_counter.py
# Created 8/16/2016, this version created 1/10/2017
# Sam Westreich, stwestreich@ucdavis.edu, github.com/transcript
#
# This program parses through the results file from a DIAMOND annotation run
# (in BLAST m8 format) to get the results into something more compressed
# and readable.
#
# Usage:
#
# -I		infile			specifies the infile (a DIAMOND results file
#								in m8 format)
# -D		database		specifies a reference database to search against
#								for results
# -O		organism		returns organism results
# -F		function		returns functional results
# -SO		specific org	creates a separate outfile for results that hit
#							a specific organism
#
##########################################################################

# imports
import operator, sys, time, gzip, re

# String searching function:
def string_find(usage_term):
	for idx, elem in enumerate(sys.argv):
		next_elem = sys.argv[(idx + 1) % len(sys.argv)]
		if elem == usage_term:
			 return next_elem


# checking for an option (organism or function) to be specified
if "-O" not in sys.argv:
	if "-F" not in sys.argv:
		sys.exit("WARNING: need to specify either organism results (with -O flag in command) or functional results (with -F flag in command).")

# loading starting file
if "-I" in sys.argv:
	infile_name = string_find("-I")
else:
	sys.exit ("WARNING: infile must be specified using '-I' flag.")

# checking to make sure database is specified
if "-D" in sys.argv:
	db_name = string_find("-D")
else:
	sys.exit( "No database file indicated; skipping database search step.")


infile = open (infile_name, "r")

# setting up databases
RefSeq_hit_count_db = {}
unique_seq_db = {}
line_counter = 0

# reading through the infile - the DIAMOND results m8 format
print "\nNow reading through the m8 results infile."
t0 = time.clock()
for line in infile:
	line_counter += 1
	splitline = line.split("\t")
	if line_counter % 1000000 == 0:
		t1 = time.clock()
		print str(line_counter)[:-6] + "M lines processed so far in " + str(t1-t0) + " seconds."

	unique_seq_db[splitline[0]] = 1

	try:
		RefSeq_hit_count_db[splitline[1]] += 1
	except KeyError:
		RefSeq_hit_count_db[splitline[1]] = 1
		continue

t2 = time.clock()
print "\nAnalysis of " + infile_name + " complete."
print "Number of total lines: " + str(line_counter)
print "Number of unique sequences: " + str(len(unique_seq_db))
print "Time elapsed: " + str(t2-t0) + " seconds."

infile.close()

# time to search for these in the reference database
print "\nStarting database analysis now."

# optional outfile of specific organism results
if "-SO" in sys.argv:
	target_org = string_find("-SO")
	db_SO_dictionary = {}

# building a dictionary of the reference database
if "-F" in sys.argv:
	db_func_dictionary = {}

if "-O" in sys.argv:
	db_org_dictionary = {}

db_line_counter = 0
db_error_counter = 0

t0 = time.clock()
db = open (db_name, "r")
for line in db:
	if line.startswith(">"):
		db_line_counter += 1
	else:
		continue
	# id, organism and function names [https://stackoverflow.com/questions/6109882/regex-match-all-characters-between-two-strings]
        db_id = re.search("(?<=>)[^ ]+", line)
	db_id = db_id.group()
        
	db_entry = re.search("(?<= )(.*)(?= \[)", line)
	if db_entry is None:
		db_entry = re.search("(?<= )(.*)(?=\[)", line) # some lines missing a space
        if db_entry is None:
		# The NR database has entries without "[", "]"
            	db_entry = line.strip('>' + db_id + ' ')
        else:
		db_entry = db_entry.group()
        
	db_org = re.search("(?<=\[)(.*)(?=\])", line)
	if db_org is None:
		# In case no labels inside "[Genus Species]"
		db_org = ''
	else:
		db_org = db_org.group()
        
	if "-F" in sys.argv:
		db_func_dictionary[db_id] = db_entry
        if "-O" in sys.argv:
		db_org_dictionary[db_id] = db_org
        if "-SO" in sys.argv:
		if target_org in db_org:
			db_SO_dictionary[db_id] = db_entry
	# line counter to show progress
        if db_line_counter % 1000000 == 0:							# each million
		t1 = time.clock()
		print str(db_line_counter)[:-6] + "M lines processed so far in " + str(t1-t0) + " seconds."

t2 = time.clock()
print "\nSuccess!"
print "Time elapsed: " + str(t2-t0) + " seconds."
print "Number of lines: " + str(db_line_counter)
print "Number of errors: " + str(db_error_counter)

# condensing down the identical matches
condensed_RefSeq_hit_db = {}
t0 = time.clock()
for entry in RefSeq_hit_count_db.keys():
	try:
		if "-O" in sys.argv:
			org = db_org_dictionary[entry]
		if "-F" in sys.argv:
			org = db_func_dictionary[entry]
		if org in condensed_RefSeq_hit_db.keys():
			condensed_RefSeq_hit_db[org] += RefSeq_hit_count_db[entry]
		else:
			condensed_RefSeq_hit_db[org] = RefSeq_hit_count_db[entry]
	except KeyError:
		print "KeyError:\t" + entry
		continue

if "SO" in sys.argv:
	condensed_RefSeq_SO_hit_db = {}

	for entry in RefSeq_hit_count_db.keys():
		if entry in db_SO_dictionary.values():
			org = db_SO_dictionary[entry]
			if org in condensed_RefSeq_SO_hit_db.keys():
				condensed_RefSeq_SO_hit_db[org] += RefSeq_hit_count_db[entry]
			else:
				condensed_RefSeq_SO_hit_db[org] = RefSeq_hit_count_db[entry]

t1 = time.clock()
# dictionary output and summary
print "\nDictionary database assembled."
print "Time elapsed: " + str(t1-t0) + " seconds."
print "Number of errors: " + str(db_error_counter)

if "-O" in sys.argv:
	print "\nTop ten organism matches:"
if "-F" in sys.argv:
	print "\nTop ten function matches:"
for k, v in sorted(condensed_RefSeq_hit_db.items(), key=lambda (k,v): -v)[:10]:
	try:
		print (str(v) + "\t" + k )
	except KeyError:
		print (str(v) + "\tWARNING: Key not found for " + k)
		continue

# creating the outfiles
if "-O" in sys.argv:
	outfile_name = infile_name[:-4] + "_organism.tsv"
if "-F" in sys.argv:
	outfile_name = infile_name[:-4] + "_function.tsv"
if "=SO" in sys.argv:
	target_org_outfile = open(infile_name[:-4] + "_" + target_org + ".tsv", "w")

outfile = open (outfile_name, "w")

# writing the output
error_counter = 0
for k, v in sorted(condensed_RefSeq_hit_db.items(), key=lambda (k,v): -v):
	try:
		q = v * 100 / float(line_counter)
		outfile.write (str(q) + "\t" + str(v) + "\t" + k + "\n")
	except KeyError:
		outfile.write (str(q) + "\t" + str(v) + "\tWARNING: Key not found for " + k + "\n")
		error_counter += 1
		continue

# writing the output if optional specific organism flag is active
if "-SO" in sys.argv:
	for k, v in sorted(condensed_RefSeq_SO_hit_db.items(), key=lambda (k,v): -v):
		try:
			q = v * 100 / float(line_counter)
			target_org_outfile.write (str(q) + "\t" + str(v) + "\t" + k + "\n")
		except KeyError:
			target_org_outfile.write (str(q) + "\t" + str(v) + "\tWARNING: Key not found for " + k + "\n")
			error_counter += 1
			continue

print "\nAnnotations saved to file: '" + outfile_name + "'."
print "Number of errors: " + str(error_counter)

db.close()
outfile.close()

