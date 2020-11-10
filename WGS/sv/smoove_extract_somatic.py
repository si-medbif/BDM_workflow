#! /usr/bin/env python

# Compares: deletions

import sys
import gzip

infile = sys.argv[1]
outfile = '{}.somatic.vcf'.format(sys.argv[1].rsplit('.',2)[0]) # Remove 'vcf.gz'

def geno_process(k,v):
    # Sample fields does not contain the keys, so FORMAT field must be provided
    db = {}
    lk = k.split(':')
    lv = v.split(':')
    for key, value in zip(lk, lv):
        db[key] = value
    return db

with gzip.open(infile, 'rt') as fin, open(outfile, 'w') as fout:
    for line in fin:
        if line.startswith('##'):
            fout.write(line)
            continue
        if line.startswith('#'):
            CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,TUMOR,NORMAL = line.strip().split()
            fout.write(line)
            continue
        CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,TGENO,NGENO = line.strip().split()
        ngeno = geno_process(FORMAT, NGENO)
        tgeno = geno_process(FORMAT, TGENO)
        if ngeno['GT'] == '0/0' and tgeno['GT'] != '0/0':
            fout.write(line)
