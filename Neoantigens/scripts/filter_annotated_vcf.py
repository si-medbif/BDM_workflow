#!/usr/bin/env python
# Filter a VCF file to only include variants on the 23 main chromosomes
# Also removes any variant not having either 'PASS' or 'germline_risk' as filter
# criteria.
import sys

allowed_chroms = ['chr'+str(i) for i in range(1,23)] + ['chrX']
allowed_filters = ['PASS','germline_risk']

with open(sys.argv[1], 'r') as fin, open(sys.argv[2], 'w') as fout:
    for line in fin:
        if line.startswith('#'):
            fout.write(line)
            continue
        l = line.strip().split()
        if l[0] in allowed_chroms and l[6] in allowed_filters:
            fout.write(line)
