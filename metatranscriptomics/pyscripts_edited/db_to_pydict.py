# imports
import argparse, time, pickle, re
import pandas as pd

def parse_header(line, db_type):
    if (db_type == "seed" or db_type == "uniref"):
        db_id = line.split()[0].strip('>')
        db_entry = line.strip(line.split()[0])
    elif (db_type == "refseq" or db_type == "nr"):
        # id, organism and function names [https://stackoverflow.com/questions/6109882/regex-match-all-characters-between-two-strings]
        db_id = re.search("(?<=>)[^ ]+", line)
        db_id = db_id.group()
        db_entry = line.split(']', 1)[0] + ']' # limit to the first match
    else:
        raise ValueError("db_type argument supplied is not supported")
    return (db_id, db_entry)
    

def write_to_file (db_id, db_entry, db_type, gene_sequence, outfile):
    gene_len = len(gene_sequence)
    if not bool(re.match('^[ATCG]+$', gene_sequence)):
        gene_len *= 3 
    if(db_type == "seed"):
        db_entry = db_entry.rsplit('\t', 1)[0]
        outfile.write(db_id + "\t" + str(gene_len) + db_entry + "\n")
    elif (db_type == "uniref"):
        db_entry = db_entry.lstrip()
        db_entry = db_entry.replace('\t', '')
        for x in [" n=", " Tax=", " TaxID=", " RepID="]:
            db_entry = db_entry.replace(x, "\t")
        outfile.write(db_id + "\t" + str(gene_len) + "\t" + db_entry)
            
    elif db_type == "refseq" or db_type == "nr":
        db_org = re.search("(?<=\[)(.*)(?=\])", db_entry)
        if db_org is None:
            db_org = '' # no genus/species assignment
        else:
            db_org = db_org.group()
        db_function = re.search("(?<= )(.*)(?= \[)", db_entry)
        if db_function is None:
            db_function = re.search("(?<= )(.*)(?=\[)", db_entry) # some lines missing a space
        if db_function is None:
            # The NR database has entries without "[", "]"
            db_function = db_entry.strip(db_entry.split()[0])
        else:
            db_function = db_function.group()
        outfile.write(db_id + "\t" + str(gene_len) + "\t" + 
                      db_org + "\t" + db_function  + "\n")
    else:
        raise ValueError("db_type argument supplied is not supported")
    return 


parser = argparse.ArgumentParser(description='Aggregate and count reads aligned with DIAMOND.')
parser.add_argument('dbfile', metavar='dbfile', type=str, nargs=1,
                    help='path to reference database. Either ' +
                          'a .fasta file or a pickle file containing python dictionary ' +
                          'with aligned function/organism ids and descriptions ' +
                          'as keys and values.')
parser.add_argument('outfile', metavar='outfile', type=str, nargs=1,
                    help='path to output data frame with gene data')
parser.add_argument('dbtype', metavar='dbtype', type=str, nargs=1,
                    help='database type; one of: refseq, seed, nr, uniref')

# parser.add_argument('--prot-seq', dest='protein_seq', 
#                     action="store_true", default=False,
#                    help='Database contains canonical protein sequences not' + 
#                         'nucleic acids. The gene lengths are multiplied by 3.')
# parser.add_argument('--seed', dest='seed',
#                    action="store_true", default=False,
#                    help='Database is the hierarchical seed db')

args = parser.parse_args()
dbtype = args.dbtype[0]
print("Reading database of type: " + dbtype)
print("Located at: " + args.dbfile[0])
print("Saving results at: " + args.outfile[0])

db = open(args.dbfile[0], "r")
outfile = open(args.outfile[0], "w")

if (dbtype == "seed"):
    outfile.write("GeneID" + "\t" + "Length" + "\t" "SEED1" + "\t" + \
                  "SEED2" + "\t" "SEED3" + "\t" "SEED4" + "\n")
if (dbtype == "uniref"):
    outfile.write("GeneID" + "\t" + "Length" + "\t" "ClusterName" + "\t" + \
                  "No_Members" + "\t" "Taxon" + "\t" "TaxID" + "\t" + "RepID" + "\n")
elif (dbtype == "refseq" or dbtype == "nr"):
    outfile.write("GeneID" + "\t" + "Length" + "\t" "Organism" + "\t" + "Function \n")
else:
    raise ValueError("Unsuported input dbtype")

t0 = time.clock()
gene_counter = 0
gene_sequence = ''

for line in db:
    if not line.startswith(">"):
        gene_sequence += line.rstrip()
        continue
    gene_counter += 1 
    if len(gene_sequence) > 0: # this is the previous gene
        write_to_file(dbid, dbentry, dbtype, gene_sequence, outfile) 
    dbid, dbentry = parse_header(line, dbtype) 
    gene_sequence = ''     
    # line counter to show progress
    if gene_counter % 1000000 == 0:     # each million
        t1 = time.clock()
        print(str(gene_counter)[:-6] + "M genes processed so far in " + \
              str(t1-t0) + " seconds.")

write_to_file(dbid, dbentry, dbtype, gene_sequence, outfile) 
outfile.close()
db.close()

if(dbtype == "seed"):
    out = pd.read_csv(args.outfile[0], sep = "\t")
    out = out.drop_duplicates('GeneID')
    out.to_csv(args.outfile[0], sep="\t", index = False)

t2 = time.clock()
print("\nSuccess!")
print("Time elapsed: " + str(t2-t0) + " seconds.")
print("Number of genes: " + str(gene_counter))

