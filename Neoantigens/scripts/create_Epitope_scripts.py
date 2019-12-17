#!/usr/bin/env python

import sys
sample = sys.argv[1]

# Collect data from HLA prediction files

def clean_hla(allele):
    a = allele.split(':')
    return ':'.join(a[0:2])

# Kourami parsing
kourami = []
infile = '{}_Kourami.txt'.format(sample)
with open(infile, 'r') as fin:
    for line in fin:
        allele, *rest = line.strip().split('\t')
        if not allele.startswith('D'):
            allele = 'HLA-' + allele
        kourami.append(clean_hla(allele))
# HLAHD
hlahd = []
infile = '{}_HLAHD.txt'.format(sample)
with open(infile, 'r') as fin:
    for line in fin:
        category, allele1, allele2 = line.strip().split('\t')
        if category in ['DRB1','DQA1','DQB1']:
            allele1 = allele1.strip('HLA-')
            allele2 = allele2.strip('HLA-')
        if category in ['A','B','C','DRB1','DQA1','DQB1']:
            hlahd.append(clean_hla(allele1))
            hlahd.append(clean_hla(allele2))
# PHLAT
phlat = []
infile = '{}_PHLAT.txt'.format(sample)
with open(infile, 'r') as fin:
    next(fin)
    for line in fin:
        category, allele1, allele2, t1,t2,t3 = line.strip().split('\t')
        if not allele1.startswith('D'):
            allele1 = 'HLA-' + allele1
        if not allele2.startswith('D'):
            allele2 = 'HLA-' + allele2
        phlat.append(clean_hla(allele1))
        phlat.append(clean_hla(allele2))
db = {}
for allele in kourami+hlahd+phlat:
    db[allele] = db.get(allele, 0) + 1
for allele, count in db.items():
    print(allele, count)

# Combines the HLA alleles from a file "{sample}.results" to the
# script template file "06_Pvac.sh.template" to create the final
# script file "06_pvac_{sample}.sh"

input_file = '06_Pvac.sh.template'
HLA_file = '{}.result'.format(sample)
script_file = '06_Pvac_{}.sh'.format(sample)

hla = []
with open(HLA_file, 'r') as fin:
    for line in fin:
        l = line.strip().split()
        try:
            allele = l[0]
        except IndexError:
            continue
        if allele[:2] in ['A*','B*','C*']:
            allele = 'HLA-' + allele
        t = ':'.join(allele.split(':')[:2])
        if t not in hla:
            hla.append(t)

with open(input_file, 'r') as fin, open(script_file, 'w') as fout:
    for line in fin:
        if line.startswith'#HLA=':
            fout.write('HLA={}\n'.format(','.join(hla)))
        elif line.startswith('#Abnm_ID='):
            fout.write('Abnm_ID={}\n'.format(sample))
        else:
            fout.write(line[1:])
