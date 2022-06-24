#!/bin/bash          

# Default filter values:
# --normal-cov 5
# --tdna-cov 10
# --trna-cov 10
# --normal-vaf 0.02
# --tdna-vaf 0.25
# --trna-vaf 0.25
# --expn-val 1
 

DIR_MHC=$1
SAMPLE=$2

docker run --rm -v ${DIR_MHC}:/mhc \
        griffithlab/pvactools:1.5.4 \
        pvacseq binding_filter \
        --binding-threshold 99999999 \
        --minimum-fold-change 0 \
        --top-score-metric median \
        /mhc/${SAMPLE}.all_epitopes.tsv \
        /mhc/${SAMPLE}.all_epitopes.binding.tsv

docker run --rm -v ${DIR_MHC}:/mhc \
        griffithlab/pvactools:1.5.4 \
        pvacseq coverage_filter \
        --normal-cov 0 \
        --tdna-cov 0 \
        --trna-cov 0 \
        --normal-vaf 1 \
        --tdna-vaf 0 \
        --trna-vaf 0 \
        --expn-val 0 \
        /mhc/${SAMPLE}.all_epitopes.binding.tsv \
        /mhc/${SAMPLE}.all_epitopes.binding.coverage.tsv

docker run --rm -v ${DIR_MHC}:/mhc \
        griffithlab/pvactools:1.5.4 \
        pvacseq transcript_support_level_filter \
        --maximum-transcript-support-level 1 \
        /mhc/${SAMPLE}.all_epitopes.binding.coverage.tsv \
        /mhc/${SAMPLE}.all_epitopes.binding.coverage.transcript.tsv

docker run --rm -v ${DIR_MHC}:/mhc \
        griffithlab/pvactools:1.5.4 \
        pvacseq top_score_filter \
        --top-score-metric median \
        /mhc/${SAMPLE}.all_epitopes.binding.coverage.transcript.tsv \
        /mhc/${SAMPLE}.man_filter.tsv
