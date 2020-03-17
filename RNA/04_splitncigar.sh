#!/bin/bash

OUT=/gnome/tmp
REFERENCE=/gnome/genome_database/ensembl

#Split and trim reads and reassign mapping qualities 
docker run --rm -v ${OUT}:/out \
	-v ${REFERENCE}:/ref \
	broadinstitute/gatk:4.1.4.1 gatk --java-options "-Xmx8G" \
	SplitNCigarReads \
	-R /ref/Homo_sapiens.GRCh38.dna.primary_assembly.convert.fa \
	-I /out/$1_dedup.bam \
	-O /out/$1_splitncigar.bam 

