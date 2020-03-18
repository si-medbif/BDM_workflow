#!/usr/bin/env python3

"""
    Collects output from single sample featureCounts analysis.
    Input:  Folder name where featureCounts output is stored
            Output files should be named as <SAMPLE>_readcounts.featurecounts.txt
    Output: One file with TPM values for each sample
"""

import sys
import glob

infolder = sys.argv[1]

outfile = infolder+'/all-samples_tpm.txt'
try:
    column = int(sys.argv[2])
except IndexError:
    column = 6

files = glob.glob(infolder+'/*readcounts.featurecounts.txt')
sum_rpk = {}
counts = []
genes = []
for index, infile in enumerate(files):
    if infile not in sum_rpk:
        sum_rpk[infile] = 0
        counts.append([])
    with open(infile, 'r') as fin:
        for line in fin:
            if line.startswith('#'):
                next(fin)
                continue
            l = line.strip().split()
            rpk = 1000*int(l[column])/int(l[5])
            name = l[0]
            sum_rpk[infile] += rpk
            counts[index].append(rpk)
            if index == 0:
                genes.append(name)
with open(outfile, 'w') as fout:
    fout.write('gene_name')
    for name in files:
        t = name.split('/')[-1].split('_')[0]
        fout.write('\t{}'.format(t))
    fout.write('\n')
    for index, name in enumerate(genes):
        fout.write(name)
        for fi_index, infile in enumerate(files):
            rpk = counts[fi_index][index]
            tpm = rpk / sum_rpk[infile] * 1000000
            fout.write('\t{:.3f}'.format(tpm))
        fout.write('\n')
