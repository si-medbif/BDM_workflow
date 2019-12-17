#!/bin/bash

# Usage: 04_SomaticVariantCall.sh <output-folder> <GATK-resource-folder> <Tumor_sample> <Normal_sample>

# Input: BAM file, base quality scores are recalibrated.
# Output: VCF file, somatic variants, filtered.


dir_Output=$1
dir_Hg38=/gnome/genome_database/gatk_bundle/hg38bundle
dir_GATK=$2
sampleT=$3
sampleN=$4

docker run --rm -v ${dir_Output}:/Output \
        -v ${dir_Hg38}:/Hg38_dir \
	-v ${dir_GATK}:/GATKtutorial \
	broadinstitute/gatk gatk --java-options "-Xmx8g" Mutect2 \
        -R /Hg38_dir/Homo_sapiens_assembly38.fasta \
        -I /Output/${sampleT}/BAM/${sampleT}_recal.bam \
        -tumor ${sampleT} \
        -I /Output/${sampleN}/BAM/${sampleN}_recal.bam \
        -normal ${sampleN} \
        -pon /GATKtutorial/1000g_pon.hg38.vcf.gz \
        --germline-resource /GATKtutorial/af-only-gnomad.hg38.vcf.gz \
        --af-of-alleles-not-in-resource 0.0000025 \
        -O /Output/${sampleT}/VCF/${sampleT}_m2.vcf.gz

docker run --rm -v ${dir_Output}:/out \
	broadinstitute/gatk gatk \
	--java-options "-Xmx8g" FilterMutectCalls \
	-V /out/${sampleT}/VCF/${sampleT}_m2.vcf.gz \
	-O /out/${sampleT}/VCF/${sampleT}_m2_filtered.vcf.gz

