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
group.add_argument("-f", "--folder", help="path to the folder with count-files")
group.add_argument("-s", "--samples", help="file with list of the count-files")
ap.add_argument("-n", "--names", help="Map for gene translation [alias1 alias2 ... new]", default="/gnome/genome_database/ensembl/ensembl_gene-transcript-name.tsv")
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
            geneid = db.get(Geneid, 'NAN')
            if geneid == 'NAN':
                continue
            if geneid not in gene_db:
                gene_db[geneid] = int(Length)
            for sample, s in zip(samples, count):
                counts[sample][geneid] = int(s)
print('Finished reading counts')

num_samples = len(counts)
outfile1 = '{}.tpm.tsv'.format(args["outfile"])
outfile2 = '{}.raw.tsv'.format(args["outfile"])
outfile3 = '{}.fpkm.tsv'.format(args["outfile"])
# Write raw count data for all samples to file
samples = counts.keys()
with open(outfile2, 'w') as fout:
    fout.write('geneid\tlength\t{}\n'.format('\t'.join(samples)))
    for gene in gene_db:
        fout.write('{}\t{}'.format(gene, gene_db[gene]))
        for sample in samples:
            try:
                fout.write('\t{}'.format(counts[sample][gene]))
            except KeyError:
                fout.write('\tNA')
        fout.write('\n')
print('Finished writing raw counts')

# Calculate TPM and write data for all samples to file
rawdf = pd.read_csv(outfile2, header=0, sep='\t', comment='#')
columns = rawdf.columns
tpmdf = rawdf.iloc[:,[0, 1]]
for col_ind in range(2, len(columns)):
    df1 = 1000 * rawdf.iloc[:,col_ind] / rawdf.iloc[:,1]
    df1 = 1000000 * df1 / df1.sum()
    tpmdf = pd.concat([tpmdf, df1], axis=1)

tpmdf.columns = columns

outdf = tpmdf.drop('length', axis=1)
outdf.to_csv(outfile1, sep='\t', index=False, na_rep = 'NaN', float_format='%.3f')

# Calculate FPKM and write data for all samples to file
rawdf = pd.read_csv(outfile2, header=0, sep='\t', comment='#')
columns = rawdf.columns
tpmdf = rawdf.iloc[:,[0, 1]]
for col_ind in range(2, len(columns)):
    df1 = rawdf.iloc[:,col_ind]
    df1 = 1000000 * df1 / df1.sum()
    df1 = 1000 * df1 / rawdf.iloc[:,1]
    tpmdf = pd.concat([tpmdf, df1], axis=1)

tpmdf.columns = columns

outdf = tpmdf.drop('length', axis=1)
outdf.to_csv(outfile3, sep='\t', index=False, na_rep = 'NaN', float_format='%.3f')

