#!/bin/bash

VCFFOLDER=$1
dir_Hg38=/gnome/genome_database/gatk_bundle/hg38bundle
SAMPLE=$2

docker run --rm -v ${VCFFOLDER}:/Output \
    -v ${dir_Hg38}:/Hg38_dir \
    broadinstitute/gatk:4.1.7.0 gatk --java-options "-Xmx16g" \
    HaplotypeCaller \
    -R /Hg38_dir/Homo_sapiens_assembly38.fasta \
    -I /Output/${SAMPLE}_recal.bam \
    -D /Hg38_dir/dbsnp_146.hg38.vcf.gz \
    -O /Output/${SAMPLE}.g.vcf \
    -ERC GVCF \
    --pcr-indel-model NONE \
    --annotation-group StandardHCAnnotation \
    --max-reads-per-alignment-start 0 \
    --native-pair-hmm-threads 1 \
    --max-mnp-distance 0
