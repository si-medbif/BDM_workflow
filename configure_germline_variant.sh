#!/bin/bash

PROJECT=$1

mkdir -p ${PROJECT}

cp DNA/alignment.sh ${PROJECT}/01_alignment.sh
cp DNA/mark_duplicate.sh ${PROJECT}/02_mark_duplicate.sh
cp DNA/recalibration.sh ${PROJECT}/03_recalibration.sh
cp DNA/combine_gvcfs.sh ${PROJECT}/04_combine_gvcfs.sh
cp DNA/genotype_gvcfs.sh ${PROJECT}/05_genotype_gvcfs.sh
cp DNA/combine_vcfs.sh ${PROJECT}/06_combine_vcfs.sh
cp DNA/variant_recalibrator.sh ${PROJECT}/07_variant_recalibrator.sh
cp DNA/annotate_vep.sh ${PROJECT}/08_annotate_vep.sh
cp reports/samtools_stats.sh ${PROJECT}/09_samtools_stats.sh
cp reports/create_dna-bam_report.py ${PROJECT}/10_create_dna-bam_report.py
