#! /usr/bin/env python

import sys
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

if len(sys.argv) != 4:
    sys.stderr.write('Usage: create_bam_report.py <sample_file> <BAM-folder> <report-file>\n')
    sys.exit(1)

samples_file = sys.argv[1]
bamfolder = sys.argv[2]
outfile = sys.argv[3]

samples = []
with open(samples_file, 'r') as fin:
    for line in fin:
        l = line.strip()
        samples.append(l)

input_reads = []
uniquely_mapped_reads = []
mapped_to_multiple_loci = []
mapped_to_too_many_loci = []
ratio = [] # (uniquely_mapped_reads + mapped_to_multiple_loci) / input_reads
    
for sample in samples:
    infile = '{}/{}_Log.final.out'.format(bamfolder, sample)
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
fig.savefig(outfile)
