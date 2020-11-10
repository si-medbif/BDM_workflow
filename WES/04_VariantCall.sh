#!/bin/bash

BAMFOLDER=$1
VCFFOLDER=$2
dir_Hg38=/gnome/genome_database/gatk_bundle/hg38bundle
SAMPLE=$3
BED=agilent_sureselect_V5_UTR_hg38_clean.bed

docker run --rm -v ${BAMFOLDER}:/Output \
	-v ${VCFFOLDER}:/vcf_folder \
	-v ${dir_Hg38}:/Hg38_dir \
    broadinstitute/gatk:4.1.9.0 gatk --java-options "-Xmx16g" HaplotypeCaller \
    -R /Hg38_dir/Homo_sapiens_assembly38.fasta \
    -I /Output/${SAMPLE}_recal.bam \
    -D /Hg38_dir/dbsnp_146.hg38.vcf.gz \
    -O /vcf_folder/${SAMPLE}.g.vcf \
    -L /Hg38_dir/${BED}
    -ERC GVCF \
    --pcr-indel-model NONE \
    --annotation-group StandardHCAnnotation \
    --max-reads-per-alignment-start 0 \
    --native-pair-hmm-threads 1 \
    --max-mnp-distance 0
