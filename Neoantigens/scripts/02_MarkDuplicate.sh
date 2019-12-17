#!/bin/bash

# Usage: 02_MarkDuplicate.sh <output-folder> <sample_name>

# Input: BAM file, sorted by coordinate
# Output: BAM file, duplicate reads are tagged

dir_Output=$1

docker run --rm -v ${dir_Output}:/Output \
	broadinstitute/picard:latest \
	MarkDuplicates \
	CREATE_INDEX=true \
	I=/Output/$2/BAM/$2_sorted.bam \
	O=/Output/$2/BAM/$2_dedupped.bam \
	M=/Output/$2/BAM/$2_dedup_output.metrics




