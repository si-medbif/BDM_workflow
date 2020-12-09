#! /usr/bin/env python

# Input is a file with two columns: sample-name VCF-file
# Output:
#   TMB

import sys
import gzip
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

project = sys.argv[1]
sample_files_list = sys.argv[2]
chroms = ['chr'+str(c) for c in range(1,23)]+['chrX','chrY']
infiles=[]
samplelist = []

syn = ['synonymous_variant']
nonsyn = ['missense_variant', 'inframe_insertion', 'inframe_deletion']
high = ['transcript_ablation','splice_acceptor_variant','splice_donor_variant','stop_gained','frameshift_variant','stop_lost','start_lost','transcript_amplification']


with open(sample_files_list, 'r') as fin:
    for line in fin:
        l = line.strip().split()
        infiles.append([l[0],l[1]])
        samplelist.append(l[0])

def info_process(s):
    db = {}
    for el in s.split(';'):
        try:
            key,value = el.split('=')
            db[key] = value
        except ValueError:
            db[el] = True
    return db

def geno_process(k,v):
    # Sample fields does not contain the keys, so FORMAT field must be provided
    db = {}
    lk = k.split(':')
    lv = v.split(':')
    for key, value in zip(lk, lv):
        db[key] = value
    return db

mutations = []
sampledb = {} # Counts of mutation category (syn, non-syn, high) for each sample 
sampledb_file = '{}-samples.tsv'.format(project)
tmb = {} # Counts of mutation category (syn, non-syn, high) for each gene
tmb_file = '{}-tmb.tsv'.format(project)
pos = {} # Counts of specific mutation for each gene
pos_file = '{}-position.tsv'.format(project)
for sample, infile in infiles:
    sampledb[sample] = {}
    if infile.endswith('gz'):
        fileopen = gzip.open
    else:
        fileopen = open
    with fileopen(infile, 'rt') as fin:
        for line in fin:
            if line.startswith('##'):
                continue
            if line.startswith('#'):
                CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,*SAMPLES = line.strip().split()
                continue
            CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT,*GENOS = line.strip().split()
            if CHROM not in chroms:
                continue
            if FILTER != 'PASS':
                continue
            info = info_process(INFO)
            annotations = info['CSQ'].split(',')
            an = annotations[0].split('|')
            # Collect mutation consequences for each sample
            mut = an[1]
            if mut in syn:
                mut = 'synonymous'
            elif mut in nonsyn:
                mut = 'non-synonymous'
            elif mut in high:
                mut = 'high'
            else:
                continue
            if mut not in mutations:
                mutations.append(mut)
            sampledb[sample][mut] = sampledb[sample].get(mut, 0) + 1
            # Collect mutation consequences for each gene
            gene = an[3]
            if gene not in tmb:
                tmb[gene] = [0,0,0] # non-syn, syn, high
            if mut == 'synonymous':
                tmb[gene][0] += 1
            elif mut == 'non-synonymous':
                tmb[gene][1] += 1
            elif mut == 'high':
                tmb[gene][2] += 1
            # Count samples per (MODERATE/HIGH) mutation position within each gene
            if mut not in ['non-synonymous','high']:
                continue
            if gene not in pos:
                pos[gene] = {}
            if int(POS) not in pos[gene]:
                pos[gene][int(POS)] = []
            pos[gene][int(POS)].append(sample)

# Write sample report
ns = []
s = []
t = []
h = []
samples = []
with open(sampledb_file,'w') as fout:
    fout.write('sample\t{}\ttotal\n'.format('\t'.join(mutations)))
    for sample, infile in infiles:
        samples.append(sample)
        fout.write('{}'.format(sample))
        total = 0
        for mut in ['synonymous','non-synonymous','high']:
            an = sampledb[sample].get(mut, 0)
            if mut == 'synonymous':
                s.append(an)
            elif mut == 'non-synonymous':
                ns.append(an)
            elif mut == 'high':
                h.append(an)
            total += an
            fout.write('\t{}'.format(an))
        t.append(total)
        fout.write('\t{}\n'.format(total))
# Create sample plot
fig, ax = plt.subplots(nrows=3, ncols=1, figsize=(14,14), sharex=True)
N = len(samples)
first = s
second = ns
third = t
ind = np.arange(N)
width = 0.35
p1 = ax[2].bar(ind, first, width, color='red')
ax[2].set_title('Synonymous')
p2 = ax[1].bar(ind, second, width, color = 'blue')
ax[1].set_title('Non-Synonymous')
p3 = ax[0].bar(ind, third, width, color = 'green')
ax[0].set_title('Total')
t1 = plt.xticks(ind,samples, rotation=90)
#l = plt.legend((p1[0], p2[0], p3[0]), ('Synonymous', 'Non-Synonymous', 'Total'))
fig.savefig(sampledb_file+'.png')

# Write gene occurence report
with open(tmb_file,'w') as fout, open(pos_file,'w') as fout2:
    fout.write('gene\tsynonymous\tnon-synonymous\thigh\n')
    for gene in tmb:
        s, ns, h = tmb[gene]
        if s + ns + h >= 1:
            fout.write('{}\t{}\t{}\t{}\n'.format(gene, s, ns, h))
        if gene in pos:
            for p, s_list in pos[gene].items():
                fout2.write('{}\t{}'.format(gene, p))
                for s in s_list:
                    fout2.write('\t{}'.format(s))
                fout2.write('\n')
