#!/usr/bin/env Python
##########################################################################
#
# Copyright (C) 2018
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
#
# This program parses through the results file from a DIAMOND annotation run
# (in BLAST m8 format) to get the results into something more compressed
# and readable.
#
##########################################################################

# imports
import argparse
import pandas as pd

parser = argparse.ArgumentParser(description='Aggregate and count reads aligned with DIAMOND.')
parser.add_argument('infile', metavar='I', type=str, nargs=1,
                    help='path to input DIAMOND alignment file in m8 format')
parser.add_argument('dbfile', metavar='DB', type=str, nargs=1,
                    help='path to reference database as .fasta file')
parser.add_argument('-outfile', dest='outfile', type=str, default=None,
                    help='path to output file. If not specified, the same as '+
			             'the infile with a replaced suffix')
args = parser.parse_args()
print(args)

if args.outfile is None:
    outdir = '.'
    outfile = args.infile[0].split(".")[0]
    outfile = os.path.join(outdir, outfile + "_gene_count.csv")
else:
    outfile = args.outfile[0]

colnames = ['qseqid', 'sseqid', 'pident', 'tot_length_aligned', 'mismatch', \
    'gapopen', 'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']

infile_kwargs = {'sep': '\t', 'chunksize': 100000, 'low_memory': False, \
    'header': None, 'names': colnames}

gene_stats = pd.DataFrame()
for chunk in pd.read_csv(args.infile[0], **infile_kwargs):
   chunk = chunk[['sseqid', 'tot_length_aligned']]
   chunk['raw'] = 1
   gene_stats = pd.concat([gene_stats, chunk], axis = 0)
   gene_stats = gene_stats.groupby('sseqid')[['tot_length_aligned', 'raw']].sum()

gene_stats.to_csv(outfile, index=False)
