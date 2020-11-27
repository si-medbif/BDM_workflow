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
    sys.stderr.write('Usage: create_fastqc_report.py <fastqc-folder> <sample1> <sample2> ... \n')
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

rows = int(np.ceil(len(samples) / 4))

fout = open('Quality_report.txt', 'w')
fig1, ax1 = plt.subplots(nrows = 1, ncols = 1, figsize=(14,14))
fig2, ax2 = plt.subplots(nrows = 1, ncols = 1, figsize=(14,14))
fig3, ax3 = plt.subplots(nrows = rows, ncols = 4, figsize=(14, 7*rows), sharex=True, sharey=True)
fig4, ax4 = plt.subplots(nrows = rows, ncols = 4, figsize=(14, 7*rows), sharex=True, sharey=True)
col = {'1':'green','2':'red'}
fails = [[],[]] # Per base sequence quality, Per sequence GC content
warns = [[],[]] # Per base sequence quality, Per sequence GC content
for index, sample in enumerate(samples):
    nc, nr = index%4, index//4
    for direction in ['1','2']:
        fastq_zip = '{}/{}_R{}_fastqc.zip'.format(fastqfolder, sample, direction)
        archive = zipfile.ZipFile(fastq_zip, 'r')
        fastqc_report1 = str(archive.read('{}_R{}_fastqc/fastqc_data.txt'.format(sample, direction)))
        baseq1 = re.search('>>Per base sequence quality(.*?)>>END_MODULE', fastqc_report1).group(1)
        gc_content1 = re.search('>>Per sequence GC content(.*?)>>END_MODULE', fastqc_report1).group(1)
        bq1 = baseq1.split('\\n')
        gc1 = gc_content1.split('\\n')
        if 'warn' in bq1[0]:
            warns[0].append(sample+'_'+direction)
        if 'warn' in gc1[0]:
            warns[1].append(sample+'_'+direction)
        if 'fail' in bq1[0]:
            fails[0].append(sample+'_'+direction)
        if 'fail' in gc1[0]:
            fails[1].append(sample+'_'+direction)
        bq_x1,bq_y1 = get_bq(bq1[2:])
        gc_x1,gc_y1 = get_gc(gc1[2:])
        x = np.arange(1,len(bq_x1)+1)
        ax1.plot(x, bq_y1, color=col[direction], label=sample+'_'+direction)
        ax2.plot(gc_y1, color=col[direction], label=sample+'_'+direction)
        if rows > 1:
            ax3[nr, nc].plot(x, bq_y1, color=col[direction], label=sample+'_'+direction)
            ax3[nr, nc].set_xlabel(sample)
            ax4[nr, nc].plot(gc_y1, color=col[direction], label=sample+'_'+direction)
            ax4[nr, nc].set_xlabel(sample)
        else:
            ax3[nc].plot(x, bq_y1, color=col[direction], label=sample+'_'+direction)
            ax3[nc].set_xlabel(sample)
            ax4[nc].plot(gc_y1, color=col[direction], label=sample+'_'+direction)
            ax4[nc].set_xlabel(sample)

fout.write('\nPer base sequence quality:\n')
if len(warns[0]) > 0:
    fout.write('Samples with warning:\n\t{}\n'.format('\n\t'.join(warns[0])))
else:
    fout.write('Samples with warning:\tNone\n')
if len(fails[0]) > 0:
    fout.write('Samples failing:\n\t{}\n'.format('\n\t'.join(fails[0])))
else:
    fout.write('Samples failing:\tNone\n')
fout.write('\nPer sequence GC content:\n')
if len(warns[1]) > 0:
    fout.write('Samples with warning:\n\t{}\n'.format('\n\t'.join(warns[1])))
else:
    fout.write('Samples with warning:\tNone\n')
if len(fails[1]) > 0:
    fout.write('Samples failing:\n\t{}\n'.format('\n\t'.join(fails[1])))
else:
    fout.write('Samples failing:\tNone\n')
fout.close()
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
if rows > 1:
    ax3[0,0].set_ylim([0,40])
    ax4[0,0].set_ylim([0,0.05])
else:
    ax3[0].set_ylim([0,40])
    ax4[0].set_ylim([0,0.05])
fig1.savefig('figure1.png')
fig2.savefig('figure2.png')
fig3.savefig('figureS1.png')
fig4.savefig('figureS2.png')
