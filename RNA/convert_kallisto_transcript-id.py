#! /usr/bin/env python

import sys
import re

gtf_file = '/gnome/genome_database/ensembl/Homo_sapiens.GRCh38.98.convert.gtf'
kallisto_output = sys.argv[1]
out_file = kallisto_output.rsplit('.',1)[0] + '.update.tsv'

db = {}

with open(gtf_file, 'r') as fin:
    for line in fin:
        try:
            gene_name = re.search(r'gene_name(.*?);',line).group(1)
            transcript_id = re.search(r'transcript_id(.*?);',line).group(1)
        except AttributeError:
            continue
        trans = transcript_id.strip().strip('"')
        gene = gene_name.strip().strip('"')
        db[trans] = gene
with open(kallisto_output, 'r') as fin, open(out_file, 'w') as fout:
    header = next(fin)
    fout.write('{}\tgene_name\n'.format(header.strip()))
    for line in fin:
        l = line.strip().split()
        trans = l[0].split('.')[0]
        out = '\t'.join(l[1:])
        try:
            gene = db[trans]
        except KeyError:
            gene = 'NA'
        fout.write('{}\t{}\t{}\n'.format(trans,out,gene))
