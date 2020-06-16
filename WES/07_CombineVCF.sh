#!/bin/bash

VCFFOLDER=$1

docker run --rm -v ${dir_Output}:/Output \
	broadinstitute/picard:latest \
	MergeVcfs \
	CREATE_INDEX=true \
	I=/Output/project.chr1.vcf.gz \
	I=/Output/project.chr2.vcf.gz \
	I=/Output/project.chr3.vcf.gz \
	I=/Output/project.chr4.vcf.gz \
	I=/Output/project.chr5.vcf.gz \
	I=/Output/project.chr6.vcf.gz \
	I=/Output/project.chr7.vcf.gz \
	I=/Output/project.chr8.vcf.gz \
	I=/Output/project.chr9.vcf.gz \
	I=/Output/project.chr10.vcf.gz \
	I=/Output/project.chr11.vcf.gz \
	I=/Output/project.chr12.vcf.gz \
	I=/Output/project.chr13.vcf.gz \
	I=/Output/project.chr14.vcf.gz \
	I=/Output/project.chr15.vcf.gz \
	I=/Output/project.chr16.vcf.gz \
	I=/Output/project.chr17.vcf.gz \
	I=/Output/project.chr18.vcf.gz \
	I=/Output/project.chr19.vcf.gz \
	I=/Output/project.chr20.vcf.gz \
	I=/Output/project.chr21.vcf.gz \
	I=/Output/project.chr22.vcf.gz \
	I=/Output/project.chrX.vcf.gz \
	I=/Output/project.chrY.vcf.gz \
	I=/Output/project.chrM.vcf.gz \
	O=/Output/project.raw.vcf.gz
