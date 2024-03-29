#!/bin/bash

OUT=$1
REFERENCE=/gnome/genome_database/ensembl
SAMPLE=$2

#Split and trim reads and reassign mapping qualities 
docker run --rm -v ${OUT}:/out \
	-v ${REFERENCE}:/ref \
	broadinstitute/gatk:4.1.9.0 gatk --java-options "-Xmx8G" \
	SplitNCigarReads \
	-R /ref/Homo_sapiens.GRCh38.dna.primary_assembly.convert.fa \
	-I /out/${SAMPLE}_dedup.bam \
	-O /out/${SAMPLE}_splitncigar.bam 

