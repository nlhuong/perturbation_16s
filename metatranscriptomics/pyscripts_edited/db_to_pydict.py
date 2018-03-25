# imports
import argparse, time, pickle, re
import pandas as pd

parser = argparse.ArgumentParser(description='Aggregate and count reads aligned with DIAMOND.')
parser.add_argument('dbfile', metavar='DB', type=str, nargs=1,
                    help='path to reference database. Either ' +
                          'a .fasta file or a pickle file containing python dictionary ' +
                          'with aligned function/organism ids and descriptions ' +
                          'as keys and values.')
parser.add_argument('outfile', metavar='outfile', type=str, nargs=1,
                    help='path to output data frame with gene data')
parser.add_argument('--prot-seq', dest='protein_seq', 
                     action="store_true", default=False,
                    help='Database contains canonical protein sequences not' + 
                         'nucleic acids. The gene lengths are multiplied by 3.')

args = parser.parse_args()

db = open(args.dbfile[0], "r")
outfile = open(args.outfile[0], "w")
outfile.write("GeneID" + "\t" + "Length" + "\t" "Organism" + "\t" + "Function \n")

t0 = time.clock()
db_line_counter = 0
gene_sequence = ''
all_genes_seq = ''
for line in db:
    if line.startswith(">"):
        db_line_counter += 1
    else:
        all_genes_seq += line.rstrip()
        gene_sequence += line.rstrip()
 	continue

    gene_len = len(gene_sequence)
    if args.protein_seq or not bool(re.match('^[ATCG]+$', gene_sequence)):
        gene_len *= 3 
    if db_line_counter > 1: # this is the previous gene
	outfile.write(db_id + "\t" + str(gene_len) + "\t" + 
                      db_org + "\t" + db_entry  + "\n")
	gene_sequence = ''
    # id, organism and function names [https://stackoverflow.com/questions/6109882/regex-match-all-characters-between-two-strings]
    db_id = re.search("(?<=>)[^ ]+", line)
    db_id = db_id.group()
    line = line.split(']', 1)[0] + ']' # limit to the first match
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
    
    # line counter to show progress
    if db_line_counter % 1000000 == 0:     # each million
        t1 = time.clock()
        print(str(db_line_counter)[:-6] + "M lines processed so far in " + \
              str(t1-t0) + " seconds.")
 
gene_len = len(gene_sequence)
if args.protein_seq or not bool(re.match('^[ATCG]+$', gene_sequence)):
    gene_len *= 3 
outfile.write(db_id + "\t" + str(gene_len)+ "\t" + db_org + "\t" + db_entry  + "\n")

t2 = time.clock()
print("\nSuccess!")
print("Time elapsed: " + str(t2-t0) + " seconds.")
print("Number of genes: " + str(db_line_counter))

db.close()
outfile.close()
