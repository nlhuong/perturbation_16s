#!/usr/bin/env python

import sys
import os
import os.path
import shutil
import subprocess
from Bio import SeqIO
from Bio.SeqRecord import SeqRecord

sequence_file = sys.argv[1]
sequences = list(SeqIO.parse(sequence_file, "fastq"))
Infernal_out = sys.argv[2]
mRNA_file = sys.argv[3]
rRNA_file = sys.argv[4]

Infernal_rRNA_IDs = set()

mRNA_seqs = {}
rRNA_seqs = {}

with open(Infernal_out, "r") as infile_read:
    for line in infile_read:
        if not line.startswith("#") and len(line) > 10:
            Infernal_rRNA_IDs.add(line[:line.find(" ")].strip())

for sequence in sequences:
    if sequence.id in Infernal_rRNA_IDs:
        rRNA_seqs[sequence.id] = sequence
    else:
        mRNA_seqs[sequence.id] = sequence

with open(mRNA_file, "w") as out:
    SeqIO.write(list(mRNA_seqs.values()), out, "fastq")

with open(rRNA_file, "w") as out:
    SeqIO.write(list(rRNA_seqs.values()), out, "fastq")

print(str(len(rRNA_seqs)) + " reads were aligned to the rRNA database")
print(str(len(mRNA_seqs)) +  " reads were not aligned to the rRNA database")
