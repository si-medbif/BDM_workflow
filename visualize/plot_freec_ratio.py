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
ratio_files = sys.argv[1]
ploidy = 2
with open(ratio_files, 'r') as fin:
    for line in fin:
        if line.startswith('#'):
            continue
        ratio_file2, sample2 = None, None
        l = line.strip().split()
        if len(l) == 4:
            sample1, sample2, ratio_file1, ratio_file2 = line.strip().split()
        elif len(l) == 2:
            sample1, ratio_file1 = line.strip().split()
        else:
            sys.exit()
        colors = ['lightblue','darkblue']
        maxval = 6
        x_locs = []
        # First sample in the pair
        print(sample1)
        df = pd.read_table(ratio_file1, header=0, dtype={'Chromosome':str})
        df.loc[df['CopyNumber']>maxval, 'CopyNumber'] = maxval
        df.loc[df['MedianRatio']>maxval, 'MedianRatio'] = maxval
        df.loc[df['Ratio']>maxval, 'Ratio'] = maxval
        offset = 0
        if ratio_file2 is None:
            fig, axes = plt.subplots(1,1,figsize=(21,7))
            ax1 = axes
            ax2 = None
        else:
            fig, axes = plt.subplots(2,1,figsize=(21,14), sharex=True)
            ax1 = axes[0]
            ax2 = axes[1]
        index = 0
        for chrom, chrom_size in zip(chroms, chrom_sizes):
            df1 = df[df['Chromosome']==chrom]
            ax1.plot(df1['Start']+offset, df1['Ratio']*ploidy, '.', ms=1, color=colors[index])
            #ax1.plot(df1['Start']+offset, df1['CopyNumber'], '.', ms=1, color='red')
            df2 = df[(df['Chromosome']==chrom) & (df['CopyNumber']>ploidy)]
            ax1.plot(df2['Start']+offset, df2['Ratio']*ploidy, '.', ms=1, color='red')
            df2 = df[(df['Chromosome']==chrom) & (df['CopyNumber']<ploidy)]
            ax1.plot(df2['Start']+offset, df2['Ratio']*ploidy, '.', ms=1, color='green')
            ax1.plot(df1['Start']+offset, df1['MedianRatio']*ploidy, '.', ms=1, color='black')
            offset += chrom_size
            x_locs.append(offset - chrom_size//2)
            ax1.vlines(offset, 0, maxval+0.1, color='black')
            ax1.set_xlim([0,offset])
            ax1.set_ylim([0,maxval+0.1])
            ax1.set_title(sample1)
            index = 1 - index

        if ratio_file2 is None:
            plt.xticks(x_locs, chroms)
            fig.savefig('{}.cnv.png'.format(sample1), bbox_inches='tight')
            plt.close()
            continue
        # Second sample in the pair
        print(sample2)
        df = pd.read_table(ratio_file2, header=0, dtype={'Chromosome':str})
        df.loc[df['CopyNumber']>maxval, 'CopyNumber'] = maxval
        df.loc[df['MedianRatio']>maxval, 'MedianRatio'] = maxval
        df.loc[df['Ratio']>maxval, 'Ratio'] = maxval
        offset = 0
        index = 0
        for chrom, chrom_size in zip(chroms, chrom_sizes):
            df1 = df[df['Chromosome']==chrom]
            ax2.plot(df1['Start']+offset, df1['Ratio'], '.', ms=1, color=colors[index])
            ax2.plot(df1['Start']+offset, df1['CopyNumber'], '.', ms=1, color='red')
            ax2.plot(df1['Start']+offset, df1['MedianRatio'], '.', ms=1, color='black')
            offset += chrom_size
            ax2.vlines(offset, 0, maxval+0.1, color='black')
            ax2.set_xlim([0,offset])
            ax2.set_ylim([0,maxval+0.1])
            ax2.set_title(sample2)
            index = 1 - index
 
        plt.xticks(x_locs, chroms)
        fig.savefig('{}_{}.cnv.png'.format(sample1, sample2))
        plt.close()
