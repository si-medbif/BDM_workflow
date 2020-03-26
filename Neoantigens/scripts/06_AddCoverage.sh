#!/bin/bash                                                                                                                                                              
set -e

DIR_OUTPUT=$1
DIR_HG38=/gnome/genome_database/gatk_bundle/hg38bundle
SAMPLE=$2

# Creates the *readcount_snv.tsv, and *readcount_indel.tsv
docker run --rm -v ${DIR_OUTPUT}:/Output \
        -v ${DIR_HG38}:/Hg38_dir \
        mgibio/bam_readcount_helper-cwl:1.0.0 \
        /usr/bin/python /usr/bin/bam_readcount_helper.py \
        /Output/${SAMPLE}/VCF/${SAMPLE}_m2_vep_filtered.vcf \
        ${SAMPLE} \
        /Hg38_dir/Homo_sapiens_assembly38.fasta \
        /Output/RNA/${SAMPLE}_Aligned.sortedByCoord.out.bam \
        /Output/${SAMPLE}/VCF/

# Add gene expression data to the vcf files
docker run --rm -v ${DIR_OUTPUT}:/Output \
        griffithlab/vatools \
        vcf-expression-annotator \
        /Output/${SAMPLE}/VCF/${SAMPLE}_m2_vep_filtered.vcf \
        /Output/RNA/${SAMPLE}_readcounts.featurecounts.txt \
        custom \
        gene \
        --id-column gene \
        --expression-column ${SAMPLE} \
        -s ${SAMPLE}

# Add RNA read counts to the vcf files
docker run --rm -v ${DIR_OUTPUT}:/Output \
        vatools:v1 \
        vcf-readcount-annotator \
        /Output/${SAMPLE}/VCF/${SAMPLE}_m2_vep_filtered.gx.vcf \
        /Output/${SAMPLE}/VCF/${SAMPLE}_bam_readcount_snv.tsv \
        RNA \
        -s ${SAMPLE} \
        -t snv \
        -o /Output/${SAMPLE}/VCF/${SAMPLE}_1.vcf

docker run --rm -v ${DIR_OUTPUT}:/Output \
        vatools:v1 \
        vcf-readcount-annotator \
        /Output/${SAMPLE}/VCF/${SAMPLE}_1.vcf \
        /Output/${SAMPLE}/VCF/${SAMPLE}_bam_readcount_indel.tsv \
        RNA \
        -s ${SAMPLE} \
        -t indel \
        -o /Output/${SAMPLE}/VCF/${SAMPLE}_m2_vep_filtered.gx.cov.vcf
