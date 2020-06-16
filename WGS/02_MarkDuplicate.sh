#!/bin/bash

BAMFOLDER=$1
SAMPLE=$2

docker run --rm -v ${BAMFOLDER}:/Output \
	broadinstitute/picard:latest \
	MarkDuplicates \
	CREATE_INDEX=true \
	I=/Output/BAM/${SAMPLE}_sorted.bam \
	O=/Output/BAM/${SAMPLE}_dedupped.bam \
	M=/Output/BAM/${SAMPLE}_dedup_output.metrics




