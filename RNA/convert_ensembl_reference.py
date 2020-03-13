#! /usr/bin/env python

import sys

fasta = '/gnome/genome_database/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.fa'
fasta_out = '/gnome/genome_database/ensembl/Homo_sapiens.GRCh38.dna.primary_assembly.convert.fa'
gtf = '/gnome/genome_database/ensembl/Homo_sapiens.GRCh38.98.gtf'
gtf_out = '/gnome/genome_database/ensembl/Homo_sapiens.GRCh38.98.convert.gtf'

chromosomes = {'1':'chr1','2':'chr2','3':'chr3','4':'chr4','5':'chr5','6':'chr6','7':'chr7',
        '8':'chr8','9':'chr9','10':'chr10','11':'chr11','12':'chr12','13':'chr13','14':'chr14',
        '15':'chr15','16':'chr16','17':'chr17','18':'chr18','19':'chr19','20':'chr20','21':'chr21',
        '22':'chr22','MT':'chrM','X':'chrX','Y':'chrY'}

chrom_order = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20',
        'chr21','chr22','chrX','chrY','chrM']

db = {}

with open(fasta, 'r') as fin, open(fasta_out, 'w') as fout:
    for line in fin:
        if line.startswith('>'):
            name = line.strip().split()[0][1:]
            try:
                newname = chromosomes[name]
            except KeyError:
                newname = name
                chrom_order.append(newname)
            db[newname] = []
        else:
            db[newname].append(line.strip())
    for chrom in chrom_order:
        fout.write('>{}\n'.format(chrom))
        for sequence in db[chrom]:
            fout.write('{}\n'.format(sequence))

db = {}
chrom_order = ['chr1','chr2','chr3','chr4','chr5','chr6','chr7','chr8','chr9','chr10',
        'chr11','chr12','chr13','chr14','chr15','chr16','chr17','chr18','chr19','chr20',
        'chr21','chr22','chrX','chrY','chrM']

with open(gtf, 'r') as fin, open(gtf_out, 'w') as fout:
    for line in fin:
        if line.startswith('#'):
            fout.write(line)
            continue
        name = line.strip().split()[0]
        try:
            newname = chromosomes[name]
        except KeyError:
            newname = name
        if newname not in db:
            db[newname] = []
        if newname not in chrom_order:
            chrom_order.append(newname)
        db[newname].append(line.strip())
    for chrom in chrom_order:
        for line in db[chrom]:
            l = line.split(None, 1)[1]
            fout.write('{}\t{}\n'.format(chrom, l))
