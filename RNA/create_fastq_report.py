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
    sys.stderr.write('Usage: create_bam_report.py <FASTQ-folder> <sample1> <sample2> ... \n')
    sys.exit(1)

fastqfolder = sys.argv[1]

samples = sys.argv[2:]

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
