#! /usr/bin/env python

import sys
import gzip

file1 = sys.argv[1]
if file1.endswith('gz'):
    open1 = gzip.open
else:
    open1 = open
name1 = file1.split('/')[-1]
file2 = sys.argv[2]
if file2.endswith('gz'):
    open2 = gzip.open
else:
    open2 = open
name2 = file2.split('/')[-1]

db1 = {} # SNP
db2 = {} # INDEL

# results = 0:PASS-PASS, 1:PASS-FAIL, 2:PASS-MISSING, 3:FAIL-PASS, 4:FAIL-FAIL, 5:FAIL-MISSING, 6:MISSING-PASS, 7:MISSING-FAIL, 8:DISCORDANT 
results1 = [0, 0, 0, 0, 0, 0, 0, 0, 0]
results2 = [0, 0, 0, 0, 0, 0, 0, 0, 0]

with open1(file1, 'rt') as fin:
    for line in fin:
        if line.startswith('#'):
            continue
        CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT, *SAMPLES = line.strip().split()
        if len(REF+ALT) == 2:
            db1[CHROM,POS] = [REF,ALT,FILTER]
        else:
            db2[CHROM,POS] = [REF,ALT,FILTER]
with open2(file2, 'rt') as fin:
    for line in fin:
        if line.startswith('#'):
            continue
        CHROM,POS,ID,REF,ALT,QUAL,FILTER,INFO,FORMAT, *SAMPLES = line.strip().split()
        if len(REF+ALT) == 2:
            db = db1
            results = results1
        else:
            db = db2
            results = results2
        try:
            r, a, f = db[CHROM,POS]
            if r != REF or a != ALT:
                results[8] += 1
            elif f == 'PASS' and FILTER == 'PASS':
                results[0] += 1
            elif f == 'PASS' and FILTER != 'PASS':
                results[1] += 1
            elif f != 'PASS' and FILTER == 'PASS':
                results[3] += 1
            elif f != 'PASS' and FILTER != 'PASS':
                results[4] += 1
            else:
                continue
            db.pop((CHROM, POS))
        except KeyError:
            if FILTER == 'PASS':
                results[6] += 1
            elif FILTER != 'PASS':
                results[7] += 1
            else:
                pass
for key, value in db1.items():
    if value[2] == 'PASS':
        results1[2] += 1
    elif value[2] != 'PASS':
        results1[5] += 1
    else:
        pass
for key, value in db2.items():
    if value[2] == 'PASS':
        results2[2] += 1
    elif value[2] != 'PASS':
        results2[5] += 1
    else:
        pass


sys.stdout.write('{}\t{}\tSNP\t{}\n'.format(name1, name2, '\t'.join([str(i) for i in results1])))
sys.stdout.write('{}\t{}\tINDEL\t{}\n'.format(name1, name2, '\t'.join([str(i) for i in results2])))
