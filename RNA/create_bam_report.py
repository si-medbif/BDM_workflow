#! /usr/bin/env python

import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import numpy as np
import matplotlib.pyplot as plt
import zipfile
import re

if len(sys.argv) != 4:
    sys.stderr.write('Usage: create_bam_report.py <BAM-folder> <sample1> <sample2> ... \n')
    sys.exit(1)

bamfolder = sys.argv[1]

samples = sys.argv[2:]

# Read STAR-aligner reports
input_reads = []
uniquely_mapped_reads = []
mapped_to_multiple_loci = []
mapped_to_too_many_loci = []
ratio = [] # (uniquely_mapped_reads + mapped_to_multiple_loci) / input_reads
    
for sample in samples:
    infile = '{}/{}/{}_Log.final.out'.format(bamfolder, sample, sample)
    with open(infile, 'r') as fin:
        for line in fin:
            try:
                key, value = line.strip().split('|')
            except ValueError:
                continue
            key = key.strip()
            value = value.strip()
            if key == 'Number of input reads':
                input_reads.append(int(value))
            elif key == 'Uniquely mapped reads number':
                uniquely_mapped_reads.append(int(value))
            elif key == 'Number of reads mapped to multiple loci':
                mapped_to_multiple_loci.append(int(value))
            elif key == 'Number of reads mapped to too many loci':
                mapped_to_too_many_loci.append(int(value))
df = pd.DataFrame({'Samples': samples,
                   'Input reads':input_reads,
                   'Uniquely mapped reads': uniquely_mapped_reads,
                   'Mapped to multiple loci': mapped_to_multiple_loci,
                   'Mapped to too many loci': mapped_to_too_many_loci
                   })
df['Ratio'] = (df['Uniquely mapped reads']+df['Mapped to multiple loci']) / df['Input reads']
df['Unmapped reads'] = df['Input reads'] - df['Uniquely mapped reads'] - df['Mapped to multiple loci']
x_size = min(20, int(len(samples)*0.5))
fig, ax = plt.subplots(nrows=1, ncols=1, figsize=(x_size,7))
df[['Samples',
    'Uniquely mapped reads',
    'Mapped to multiple loci',
    'Unmapped reads']].set_index('Samples').plot(kind='bar', stacked=True, ax = ax)
fig.savefig('figure3.png')
