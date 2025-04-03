#!/usr/bin/env python

""" Extract a 1Mbp sequence from chromsome 21 and export an adapted GFF file """

from Bio import SeqIO
import operator
from functools import reduce
from itertools import batched

fasta_chunk_size = 70

input_seq_file = "data/reference/Homo_sapiens/chr21.fasta"
input_annotation_file = "data/reference/Homo_sapiens/chr21.gff3"

output_seq_file = "data/reference/Homo_sapiens/chr21_small.fasta"
output_annotation_file = "data/reference/Homo_sapiens/chr21_small.gff3"

pos_begin   = 16_000_000
pos_end     = 17_000_000

record = list(SeqIO.parse(input_seq_file, "fasta"))[0]

with open(output_seq_file, "w") as outfile:
    outfile.write('>' + record.id + '\n')
    for chunk in batched(record.seq[pos_begin:(pos_end + 1)], fasta_chunk_size):
        outfile.write(reduce(operator.add, chunk))
        outfile.write('\n')

with open(output_annotation_file, "w") as outfile:
    with open(input_annotation_file) as infile:
        for line in infile:
            if line[0] == '#':
                outfile.write(line)
                continue

            fields = line.split('\t')
            newstart = int(fields[3]) - pos_begin

            if newstart < 0:
                continue
            
            fields[3] = str(newstart)
            fields[4] = str(int(fields[4]) - pos_begin)

            outfile.write('\t'.join(fields))
