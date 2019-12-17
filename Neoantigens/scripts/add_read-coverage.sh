#!/bin/bash

# Extract read count data from BAM files, at the variant positions indicated in the VCF file
# The docker image mgibio/bam_readcount_helper-cwl:1.0.0 runs the command:
# $ bam-readcount -f <reference_fasta> -l <site_list> <bam_file> [-i] [-b 20]

dir_BAM=PATH_TO/BAM
dir_Reference=PATH_TO/hg38bundle
dir_Output=PATH_TO/results
mkdir -p ${dir_Output}/$1/coverage/{DNA,RNA}

docker run --rm -v ${dir_Reference}:/Reference -v ${dir_Output}:/Output \
	mgibio/bam_readcount_helper-cwl:1.0.0 \
	/usr/bin/python /usr/bin/bam_readcount_helper.py \
	/Output/$1/$1_m2_vep_filtered.vcf \
	$1 \
	/Reference/Homo_sapiens_assembly38.fasta \
	/Output/$1/BAM/$1_recal.bam \
	/Output/$1/coverage/DNA

docker run --rm -v ${dir_Reference}:/Reference -v ${dir_Output}:/Output \
	mgibio/bam_readcount_helper-cwl:1.0.0 \
	/usr/bin/python /usr/bin/bam_readcount_helper.py \
	/Output/$1/$1_m2_vep_filtered.vcf \
	$1 \
	/Reference/Homo_sapiens_assembly38.fasta \
	/Output/$1/BAM/$1Aligned.sortedByCoord.out.bam \
	/Output/$1/coverage/RNA

# Add gene expression data to the vcf files
vcf-expression-annotator \
	$1_m2_vep_filtered.vcf \
	$1_counts.txt \
	custom \
	gene \
	--id-column gene_id \
	--expression-column $1 \
	-s $1

# Sample EFDFC8
# Add tumour DNA read counts to the vcf files
vcf-readcount-annotator \
	$1_m2_vep_filtered.gx.vcf \
	${dir_Output}/coverage/DNA/$1_bam_readcount_snv.tsv \
	DNA \
	-s $1 \
	-t snv \
	-o $1_1.vcf
vcf-readcount-annotator \
	$1_1.vcf \
	${dir_Output}/coverage/DNA/$1_bam_readcount_indel.tsv \
	DNA \
	-s $1 \
	-t indel \
	-o $1_2.vcf

# Add RNA read counts to the vcf files
vcf-readcount-annotator \
	$1_2.vcf \
	${dir_Output}/coverage/RNA/$1_bam_readcount_snv.tsv \
	RNA \
	-s $1 \
	-t snv \
	-o $1-T_3.vcf
vcf-readcount-annotator \
	$1_3.vcf \
	${dir_Output}/coverage/RNA/$1_bam_readcount_indel.tsv \
	RNA \
	-s $1 \
	-t indel \
	-o $1_m2_vep_filtered.gx.cov.vcf
rm $1_?.vcf

