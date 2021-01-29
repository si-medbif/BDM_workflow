#!/usr/bin/env python3

"""
    Collects output from single sample featureCounts analysis, combines output from all samples into one report per output type (TPM/FPKM/Raw).
    Input:  Either: 
              Folder name where featureCounts output is stored as <SAMPLE>_readcounts.featurecounts.txt
            Or:
              File with list of featureCounts output files. 
    Output: One file with TPM values
            One file with FPKM values
            One file with raw counts
"""

import sys
import glob
import pandas as pd
import argparse

# construct the argument parse and parse the arguments
ap = argparse.ArgumentParser()
group = ap.add_mutually_exclusive_group()
group.add_argument("-s", "--samples", help="file with list of the count-files")
group.add_argument("-f", "--files", help="path to count-file")
ap.add_argument("-n", "--names", help="Map for gene translation [alias1 alias2 ... new]", default="/gnome/genome_database/ensembl/ensembl_gene-transcript-name.tsv")
ap.add_argument("-o", "--outfile", help="Output file prefix", default="project")
args = vars(ap.parse_args())


files = []

if args["samples"] is not None:
    with open(args["samples"], 'r') as fin:
        for line in fin:
            files.append(line.strip())
elif args['files'] is not None:
    files.append(args['files'])
else:
    sys.stderr.write('ERROR: missing or incorrect input.\n')
    sys.exit(1)

db = {}
if args["names"] is not None:
    with open(args["names"], 'r') as fin:
        for line in fin:
            l = line.strip().split()
            for el in l:
                db[el] = l[-1]
counts = {}
gene_db = {}
for findex, infile in enumerate(files):
    with open(infile, 'r') as fin:
        header = next(fin)
        if header.startswith('#'):
            header = next(fin)
        samples = header.split()[6:]
        for sample in samples:
            counts[sample] = {}
        for line in fin:
            l = line.strip().split()
            Geneid, Chr, Start, End, Strand, Length, *count = l
            if Geneid not in gene_db:
                gene_db[Geneid] = int(Length)
            for sample, s in zip(samples, count):
                counts[sample][Geneid] = int(s)
print('Finished reading counts')

num_samples = len(counts)
outfile1 = '{}.tpm.tsv'.format(args["outfile"])
outfile2 = '{}.raw.tsv'.format(args["outfile"])
outfile3 = '{}.fpkm.tsv'.format(args["outfile"])
# Write raw count data for all samples to file
samples = counts.keys()
with open(outfile2, 'w') as fout:
    outsamples = []
    for sample in samples:
        outsamples.append(sample.split('/')[-1].rsplit('_',1)[0])
    fout.write('geneid\tgene\tlength\t{}\n'.format('\t'.join(outsamples)))
    for geneid in gene_db:
        gene = db.get(geneid, 'NAN')
        fout.write('{}\t{}\t{}'.format(geneid, gene, gene_db[geneid]))
        for sample in samples:
            try:
                fout.write('\t{}'.format(counts[sample][geneid]))
            except KeyError:
                fout.write('\tNA')
        fout.write('\n')
print('Finished writing raw counts')

# Calculate TPM and write data for all samples to file
rawdf = pd.read_csv(outfile2, header=0, sep='\t', comment='#')
columns = rawdf.columns
tpmdf = rawdf.iloc[:,[0, 1, 2]]
for col_ind in range(3, len(columns)):
    df1 = 1000 * rawdf.iloc[:,col_ind] / rawdf.iloc[:,2]
    df1 = 1000000 * df1 / df1.sum()
    tpmdf = pd.concat([tpmdf, df1], axis=1)

tpmdf.columns = columns

outdf = tpmdf.drop('length', axis=1)
outdf.to_csv(outfile1, sep='\t', index=False, na_rep = 'NaN', float_format='%.3f')

# Calculate FPKM and write data for all samples to file
rawdf = pd.read_csv(outfile2, header=0, sep='\t', comment='#')
columns = rawdf.columns
tpmdf = rawdf.iloc[:,[0, 1, 2]]
for col_ind in range(3, len(columns)):
    df1 = rawdf.iloc[:,col_ind]
    df1 = 1000000 * df1 / df1.sum()
    df1 = 1000 * df1 / rawdf.iloc[:,2]
    tpmdf = pd.concat([tpmdf, df1], axis=1)

tpmdf.columns = columns

outdf = tpmdf.drop('length', axis=1)
outdf.to_csv(outfile3, sep='\t', index=False, na_rep = 'NaN', float_format='%.3f')

