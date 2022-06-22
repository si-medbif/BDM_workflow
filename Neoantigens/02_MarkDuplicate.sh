#!/bin/bash

# Usage: 02_MarkDuplicate.sh <output-folder> <sample_name>

# Input: BAM file, sorted by coordinate
# Output: BAM file, duplicate reads are tagged

DIR_OUTPUT=$1
SAMPLE=$2

docker run --rm -v ${DIR_OUTPUT}:/Output \
	broadinstitute/picard:latest \
	MarkDuplicates \
	CREATE_INDEX=true \
	I=/Output/${SAMPLE}.sorted.bam \
	O=/Output/${SAMPLE}.dedupped.bam \
	M=/Output/${SAMPLE}.dedup_output.metrics




