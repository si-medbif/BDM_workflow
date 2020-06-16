#!/bin/bash

# Usage: 02_MarkDuplicate.sh <output-folder> <sample_name>

# Input: BAM file, sorted by coordinate
# Output: BAM file, duplicate reads are tagged

BAMFOLDER=$1
SAMPLE=$2

docker run --rm -v ${BAMFOLDER}:/Output \
	broadinstitute/picard:latest \
	MarkDuplicates \
	CREATE_INDEX=true \
	I=/Output/${SAMPLE}_sorted.bam \
	O=/Output/${SAMPLE}_dedupped.bam \
	M=/Output/${SAMPLE}_dedup_output.metrics




