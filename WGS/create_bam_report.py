#! /usr/bin/env python

import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import zipfile
import re

if len(sys.argv) < 3:
    sys.stderr.write('Usage: create_bam_report.py <BAM-folder> <sample1> <sample2> ...\n')
    sys.exit(1)

bamfolder = sys.argv[1]
samples = sys.argv[2:]


# Read samtools stats reports
input_reads = []
uniquely_mapped_reads = []
mapped_to_multiple_loci = []
mapped_to_too_many_loci = []
ratio = [] # (uniquely_mapped_reads + mapped_to_multiple_loci) / input_reads
    
for sample in samples:
    infile = '{}/{}_sorted.bam.stats.txt'.format(bamfolder, sample)
    with open(infile, 'r') as fin:
        for line in fin:
            try:
                l = line.strip().split('\t')
            except ValueError:
                continue
            tag = l[0]
            key = l[1]
            value = l[2]
            if tag != 'SN':
                continue
            if key == 'raw total sequences:':
                input_reads.append(int(value))
            elif key == 'reads properly paired:':
                uniquely_mapped_reads.append(int(value))
            elif key == 'reads duplicated:':
                mapped_to_multiple_loci.append(int(value))
            elif key == 'reads MQ0:':
                mapped_to_too_many_loci.append(int(value))
df = pd.DataFrame({'Samples': samples,
                   'Raw total sequences':input_reads,
                   'Properly paired reads': uniquely_mapped_reads,
                   'Reads marked as duplicates': mapped_to_multiple_loci,
                   'Reads with MQ0': mapped_to_too_many_loci
                   })
df['Ratio'] = (df['Properly paired reads']+df['Reads with MQ0']) / df['Raw total sequences']
df['Unmapped reads'] = df['Raw total sequences'] - df['Properly paired reads'] - df['Reads with MQ0']
x_size = min(20, int(len(samples)*0.5))
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(x_size,7))
df[['Samples',
    'Properly paired reads',
    'Reads marked as duplicates',
    'Unmapped reads']].set_index('Samples').plot(kind='bar', stacked=True, ax = ax)
fig.savefig('figure3.png')
