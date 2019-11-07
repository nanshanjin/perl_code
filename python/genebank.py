#!/usr/bin/env python
# -*- coding: utf-8 -*-

from Bio import SeqIO
import sys

gb_file = sys.argv[1]
gb_outfile = sys.argv[2]
f_out = open(gb_outfile, 'w')
f_out.write('DEFINITION' + '\t' + 'ACCESSION' + '\t' + 'ORGANISM' + '\t' + 'AUTHORS' + '\t' + 'TITLE' + '\t' +
            'JOURNAL' + '\t' + 'organism' + '\t' + 'strain' + '\t' + 'product' + '\t' + 'ORIGIN' + '\n')
for gb_record in SeqIO.parse(open(gb_file, "r"), "genbank"):
    seq = ''
    for i in range(len(gb_record.seq)):
        if i % 60 == 0 and i == 0:
            seq += str(i + 1)
        if i % 60 == 0 and i > 0:
            seq += ' ' + str(i + 1)
        if i % 10 == 0:
            seq += ' '
        seq += gb_record.seq[i].lower()
    if 'strain' not in gb_record.features[0].qualifiers:
        strain = ''
    else:
        strain = gb_record.features[0].qualifiers['strain'][0]
    if len(gb_record.features)<3 or  'product' not in gb_record.features[2].qualifiers:
        product = ''
    else:
        product = gb_record.features[2].qualifiers['product'][0]
    out = [gb_record.description, gb_record.annotations['accessions'][0],
           ';'.join([str(k) for k in gb_record.annotations['taxonomy']]),
           gb_record.annotations['references'][0].authors, gb_record.annotations['references'][0].title,
           gb_record.annotations['references'][0].journal, gb_record.features[0].qualifiers['organism'][0],
           strain,
           product, seq]
    f_out.write('\t'.join([str(k) for k in out]) + '\n')

