#!/usr/bin/env python3

"""
    Collects output from single sample featureCounts analysis.
    Input:  Either: 
              Folder name where featureCounts output is stored as <SAMPLE>_readcounts.featurecounts.txt
            Or:
              File with list of featureCounts output files. 
    Output: One file with TPM values for each sample
            One file with raw counts for each sample
"""

import sys
import glob
import pandas as pd
import argparse

# construct the argument parse and parse the arguments
ap = argparse.ArgumentParser()
group = ap.add_mutually_exclusive_group()
group.add_argument("-f", "--folder", help="path to the folder with count-files")
group.add_argument("-s", "--samples", help="file with list of the count-files")
group.add_argument("-m", "--multi", help="Multi-sample file")
ap.add_argument("-n", "--names", help="Map for gene translation [alias1 alias2 ... new]")
ap.add_argument("-o", "--outfile", help="Output file prefix", default="project")
args = vars(ap.parse_args())

if args["samples"] is not None:
    with open(args["samples"], 'r') as fin:
        files = []
        for line in fin:
            files.append(line.strip())
elif args['folder'] is not None:
    files = glob.glob(args["folder"]+'/*readcounts.featurecounts.txt')
    out_folder = args["folder"]
    project_name = 'project'
elif args['multi'] is not None:
    pass
else:
    sys.stderr.write('ERROR: missing or incorrect input.\n')
    sys.exit(1)

db = {}
if args["names"] is not None:
    with open(args["names"], 'r') as fin:
        for line in fin:
            l = line.strip().split()
            for el in l[0:-1]:
                db[el] = l[-1]

if args['multi'] is not None:
    indf = pd.read_csv(args['multi'], header=0, sep='\t', comment='#')
    rawdf = indf.drop(['Chr','Start','End','Strand'], axis=1)
    num_samples = len(rawdf.columns) - 2
else:
    num_samples = len(files)
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

outfile1 = '{}_{}-samples_tpm.tsv'.format(args["outfile"], num_samples)
outfile2 = '{}_{}-samples_raw.tsv'.format(args["outfile"], num_samples)
if len(db) > 0:
    rawdf['Geneid'] = rawdf['Geneid'].map(db)

# Write raw count data for all samples to file
outdf = rawdf.drop('Length', axis=1)
outdf.to_csv(outfile2, sep='\t', index=False, float_format='%.3f')

# Calculate TPM and write data for all samples to file
columns = rawdf.columns
tpmdf = rawdf.iloc[:,[0, 1]]
for col_ind in range(2, num_samples+2):
    df1 = 1000 * rawdf.iloc[:,col_ind] / rawdf.iloc[:,1]
    df1 = 1000000 * df1 / sum(df1)
    tpmdf = pd.concat([tpmdf, df1], axis=1)

tpmdf.columns = columns

outdf = tpmdf.drop('Length', axis=1)
outdf.to_csv(outfile1, sep='\t', index=False, na_rep = 'NaN', float_format='%.3f')

