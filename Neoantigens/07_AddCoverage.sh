#!/bin/bash                                                                                                                                                              
set -e

# What pVACSec is looking for in the VCF file:
# Type 			: VCF sample 				: Format Fields
# Tumor RNA Gene Expr.	: 					: GX
# Tumor DNA Coverage 	: single-sample VCF or sample_name 	: AD, DP, and AF
# Tumor RNA Coverage 	: single-sample VCF or sample_name 	: RAD, RDP, and RAF
# Normal DNA Coverage 	: --normal-sample-name 			: AD, DP, and AF

DIR_VCF=$1
DIR_RNABAM=$2
DIR_TPM=$3
DIR_HG38=/gnome/genome_database/gatk_bundle/hg38bundle
RNASAMPLE=$4
TUMOR=$5

# Creates the *readcount_snv.tsv, and *readcount_indel.tsv from RNA BAM-files.
docker run --rm -v ${DIR_VCF}:/vcf \
        -v ${DIR_HG38}:/Hg38_dir \
	-v ${DIR_RNABAM}:/rnabam \
        mgibio/bam_readcount_helper-cwl:1.0.0 \
        /usr/bin/python /usr/bin/bam_readcount_helper.py \
        /vcf/${TUMOR}.vep.filtered.vcf \
        ${TUMOR} \
        /Hg38_dir/Homo_sapiens_assembly38.fasta \
        /rnabam/${RNASAMPLE}_Aligned.sortedByCoord.out.bam \
        /vcf/

# Add gene expression data to the vcf files
docker run --rm -v ${DIR_VCF}:/vcf \
	-v ${DIR_TPM}:/tpm \
        griffithlab/vatools \
        vcf-expression-annotator \
        /vcf/${TUMOR}.vep.filtered.vcf \
        /tpm/${RNASAMPLE}.readcounts.featurecounts.txt \
        custom \
        gene \
        --id-column gene \
        --expression-column ${RNASAMPLE} \
        -s ${TUMOR}

# Add RNA read counts to the vcf files
docker run --rm -v ${DIR_VCF}:/vcf \
        vatools:v1 \
        vcf-readcount-annotator \
        /vcf/${TUMOR}.vep.filtered.gx.vcf \
        /vcf/${TUMOR}_bam_readcount_snv.tsv \
        RNA \
        -s ${TUMOR} \
        -t snv \
        -o /vcf/${TUMOR}_1.vcf

docker run --rm -v ${DIR_VCF}:/vcf \
        vatools:v1 \
        vcf-readcount-annotator \
        /vcf/${TUMOR}_1.vcf \
        /vcf/${TUMOR}_bam_readcount_indel.tsv \
        RNA \
        -s ${TUMOR} \
        -t indel \
        -o /vcf/${TUMOR}.vep.filtered.gx.cov.vcf

