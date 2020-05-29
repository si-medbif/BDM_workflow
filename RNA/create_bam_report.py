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
    sys.stderr.write('Usage: create_bam_report.py <sample_file> <BAM-folder> <FASTQ-folder>\n')
    sys.exit(1)

samples_file = sys.argv[1]
bamfolder = sys.argv[2]
fastqfolder = sys.argv[3]

samples = []
with open(samples_file, 'r') as fin:
    for line in fin:
        l = line.strip()
        samples.append(l)

# Read FASTQC reports

def get_bq(lines):
    x = []
    y = []
    for line in lines:
        l = line.split('\\t')
        if len(l) < 2:
            continue
        try:
            x.append(l[0])
            y.append(float(l[1]))
        except IndexError:
            continue
    return x,y

def get_gc(lines):
    x = []
    y = []
    for line in lines:
        l = line.split('\\t')
        if len(l) < 2:
            continue
        try:
            x.append(l[0])
            y.append(float(l[1]))
        except IndexError:
            continue
    y = np.array(y) / sum(y)
    return x,y


fig1, ax1 = plt.subplots(nrows = 1, ncols = 1, figsize=(14,14))
fig2, ax2 = plt.subplots(nrows = 1, ncols = 1, figsize=(14,14))
col = {'1':'green','2':'red'}
for sample in samples:
    for direction in ['1','2']:
        fastq_zip = '{}/{}_{}_fastqc.zip'.format(fastqfolder, sample, direction)
        archive = zipfile.ZipFile(fastq_zip, 'r')
        fastqc_report1 = str(archive.read('{}_{}_fastqc/fastqc_data.txt'.format(sample, direction)))
        baseq1 = re.search('>>Per base sequence quality(.*?)>>END_MODULE', fastqc_report1).group(1)
        gc_content1 = re.search('>>Per sequence GC content(.*?)>>END_MODULE', fastqc_report1).group(1)
        bq1 = baseq1.split('\\n')
        gc1 = gc_content1.split('\\n')
        bq_x1,bq_y1 = get_bq(bq1[2:])
        gc_x1,gc_y1 = get_gc(gc1[2:])
        x = np.arange(1,len(bq_x1)+1)
        ax1.plot(x, bq_y1, color=col[direction], label=sample+'_'+direction)
        ax2.plot(gc_y1, color=col[direction], label=sample+'_'+direction)
ax1.set_title('FastQC: Mean Quality Scores')
ax1.set_xticks(np.arange(1,len(bq_x1)+1,5))
ax1.set_xticklabels(bq_x1[0::5], rotation=0)
ax1.set_xlabel('Position (bp)')
ax1.set_ylabel('Phred Score')
ax1.set_xlim([0,len(bq_x1)+1])
ax1.set_ylim([0,40])
ax2.set_title('FastQC: Per Sequence GC Content')
ax2.set_xlabel('% GC')
ax2.set_ylabel('Percentage')
fig1.savefig('figure1.png')
fig2.savefig('figure2.png')

# Read STAR-aligner reports
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
fig.savefig('figure3.png')
