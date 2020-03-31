#!/bin/bash

dir_Output=$1
SAMPLE=$2

docker run --rm -v ${dir_Output}:/Output \
	broadinstitute/picard:latest \
	MarkDuplicates \
	CREATE_INDEX=true \
	I=/Output/${SAMPLE}/BAM/${SAMPLE}_sorted.bam \
	O=/Output/${SAMPLE}/BAM/${SAMPLE}_dedupped.bam \
	M=/Output/${SAMPLE}/BAM/${SAMPLE}_dedup_output.metrics




