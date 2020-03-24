#!/usr/bin/env python3

"""
    Collects output from single sample featureCounts analysis.
    Input:  Either: 
              Folder name where featureCounts output is stored as <SAMPLE>_readcounts.featurecounts.txt
            Or:
              File with list of featureCounts output files. First line is folder to store output.
    Output: One file with TPM values for each sample
            One file with raw counts for each sample
"""

import sys
import glob
import pandas as pd

inparam = sys.argv[1]

try:
    with open(inparam, 'r') as fin:
        files = []
        for line in fin:
            l = line.strip().split()
            if len(l) == 2:
                if l[0] == 'PROJECT':
                    project_name = l[1]
                elif l[0] == 'OUT':
                    out_folder = l[1]
                elif l[0] == 'FILE':
                    files.append(l[1])
            else:
                pass
except IsADirectoryError:
    files = glob.glob(inparam+'/*readcounts.featurecounts.txt')
    out_folder = inparam
    project_name = 'project'

outfile1 = '{}/{}_{}-samples_tpm.tsv'.format(out_folder, project_name, len(files))
outfile2 = '{}/{}_{}-samples_raw.tsv'.format(out_folder, project_name, len(files))

for findex, infile in enumerate(files):
    infile.strip().rsplit('_',1)[0] + 'counts.tsv'
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
outdf = rawdf.drop('Length', axis=1)
outdf.to_csv(outfile2, sep='\t', index=False, float_format='%.3f')
tpmdf = rawdf.iloc[:,[0, 1]]

for col_ind in range(2, len(files)+2):
    df1 = 1000 * rawdf.iloc[:,col_ind] / rawdf.iloc[:,1]
    df1 = 1000000 * df1 / sum(df1)
    tpmdf = pd.concat([tpmdf, df1], axis=1)

tpmdf.columns = columns

outdf = tpmdf.drop('Length', axis=1)
outdf.to_csv(outfile1, sep='\t', index=False, float_format='%.3f)

