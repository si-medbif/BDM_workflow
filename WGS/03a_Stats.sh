#!/bin/bash
# In version 1.9: -p will remove overlaps of pair-end reads from counting
BAMFOLDER=$1
SAMPLE=$2
samtools stats \
	-p \
	--threads 8 \
	${BAMFOLDER}/${SAMPLE}/BAM/${SAMPLE}_recal.bam \
        > ${BAMFOLDER}/${SAMPLE}/BAM/${SAMPLE}_recal.bam.stats.txt
