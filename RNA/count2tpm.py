#!/usr/bin/env python3

"""
    Collects output from single sample featureCounts analysis.
    Input:  Folder name where featureCounts output is stored
            Output files should be named as <SAMPLE>_readcounts.featurecounts.txt
    Output: One file with TPM values for each sample
            One file with raw counts for each sample
"""

import sys
import glob
import pandas as pd

infolder = sys.argv[1]

outfile1 = infolder+'/all-samples_tpm.txt'
outfile2 = infolder+'/all-samples_raw.txt'

files = glob.glob(infolder+'/*readcounts.featurecounts.txt')

for findex, infile in enumerate(files):
    sample = infile.split('/')[-1].rsplit('_',1)[0]
    indf = pd.read_csv(infile, header=0, sep='\t', comment='#')
    if findex == 0:
        rawdf = indf.iloc[:,[0, 5, 6]]
        columns = ['Geneid', 'Length', sample]
        continue
    df1 = indf.iloc[:,6]
    columns.append(sample)
    rawdf = pd.concat([rawdf, df1], axis=1)
rawdf.columns = columns

rawdf.to_csv(outfile2, sep='\t', index=False)
tpmdf = rawdf.iloc[:,[0, 1]]

for col_ind in range(2, len(files)+2):
    df1 = 1000 * rawdf.iloc[:,col_ind] / rawdf.iloc[:,1]
    df1 = 1000000 * df1 / sum(df1)
    tpmdf = pd.concat([tpmdf, df1], axis=1)

tpmdf.columns = columns
tpmdf.to_csv(outfile1, sep='\t', index=False)

