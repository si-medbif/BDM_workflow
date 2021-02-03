#! /usr/bin/env python

import sys
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import itertools
import re
import numpy as np

chroms = ['1','2','3','4','5','6','7','8','9','10','11','12','13','14','15','16','17','18','19','20','21','22','X']
chrom_sizes = [248956422,242193529,198295559,190214555,
        181538259,170805979,159345973,145138636,
        138394717,133797422,135086622,133275309,
        114364328,107043718,101991189,90338345,
        83257441,80373285,58617616,64444167,
        46709983,50818468,156040895]
bed_file = sys.argv[1]
ploidy = 2
colors = ['red','green']
x_locs = []
df = pd.read_table(bed_file, header=None, names=['chromosome','start','stop','value'], dtype={'chromosome':str})
df['pos'] = (df['start'] + df['stop']) // 2
df['size'] = df['stop'] - df['start']
maxval = df['value'].max()
offset = 0
fig, axes = plt.subplots(1,1,figsize=(21,7))
index = 0
for chrom, chrom_size in zip(chroms, chrom_sizes):
    c = 'chr{}'.format(chrom)
    df1 = df[df['chromosome']==c]
    axes.plot(df1['pos']+offset, df1['value'], '.', ms=2, color=colors[index])
    offset += chrom_size
    x_locs.append(offset - chrom_size//2)
    axes.vlines(offset, 0, maxval, color='black')
    axes.set_xlim([0,offset])
    #axes.set_ylim([0,maxval+0.1])
    #axes.set_title(sample1)
    index = 1 - index

plt.xticks(x_locs, chroms)
fig.savefig('{}.png'.format(bed_file), bbox_inches='tight')
plt.close()
