#!/bin/bash

REFERENCE=/gnome/genome_database/ensembl
OUTPUT=/gnome/chanitra_neoantigens/RNA2
FASTQ=/gnome/chanitra_neoantigens/fastq

docker run --rm -v ${REFERENCE}:/ref \
	-v ${OUTPUT}:/out \
	-v ${FASTQ}:/fq \
	zlskidmore/kallisto:latest kallisto quant \
	-i /ref/Homo_sapiens.GRCh38.cdna.all.fa.kallisto.indexed \
	-o /out/output_$1 \
	-b 100 \
	/fq/$1_1.fastq.gz /fq/$1_2.fastq.gz
