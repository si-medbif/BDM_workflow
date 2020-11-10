#! /usr/bin/env python

# Compares: deletions

import sys
import gzip
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

sample_files_list = sys.argv[1]
chroms = ['chr'+str(c) for c in range(1,23)]+['chrX','chrY']
infiles=[]

syn = ['synonymous_variant']
nonsyn = ['missense_variant', 'inframe_insertion', 'inframe_deletion']
high = ['transcript_ablation','splice_acceptor_variant','splice_donor_variant','stop_gained','frameshift_variant','stop_lost','start_lost','transcript_amplification']

with open(sample_files_list, 'r') as fin:
    for line in fin:
        l = line.strip().split()
        infiles.append([l[1],l[0]])

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
db = {}
genes = {}
for sample, infile in infiles:
    db[sample] = {}
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
            for annotation in annotations:
                an = annotation.split('|')
                # Consequence
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
                # Gene
                gene = an[3]
                if gene not in genes:
                    genes[gene] = [0,0,0] # non-syn, syn, high
                break # Assuming the first annotation has the highest impact
            # Consequence
            db[sample][mut] = db[sample].get(mut, 0) + 1
            # Gene
            if mut == 'synonymous':
                genes[gene][0] += 1
            elif mut == 'non-synonymous':
                genes[gene][1] += 1
            elif mut == 'high':
                genes[gene][2] += 1
# Write TMB report
ns = []
s = []
t = []
h = []
samples = []
sys.stdout.write('sample\t{}\ttotal\n'.format('\t'.join(mutations)))
for sample, infile in infiles:
    samples.append(sample)
    sys.stdout.write('{}'.format(sample))
    total = 0
    for mut in ['synonymous','non-synonymous','high']:
        an = db[sample].get(mut, 0)
        if mut == 'synonymous':
            s.append(an)
        elif mut == 'non-synonymous':
            ns.append(an)
        elif mut == 'high':
            h.append(an)
        total += an
        sys.stdout.write('\t{}'.format(an))
    t.append(total)
    sys.stdout.write('\t{}\n'.format(total))
# Create TMB plot
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
fig.savefig('tmb.png')

# Write gene occurence report
for gene in genes:
    s, ns, h = genes[gene]
    if s + ns + h > 1:
        sys.stdout.write('{}\t{}\t{}\t{}\n'.format(gene, s, ns, h))
