#! /bin/bash

set -e
set -u
set -o pipefail

MANTA=/home/harald/manta-1.6.0.centos6_x86_64/bin
REFERENCE=/gnome/genome_database/gatk_bundle/hg38bundle
BAM=$1
OUT=$2
SAMPLE1=$3
SAMPLE2=$4

python2 ${MANTA}/configManta.py \
	--referenceFasta=${REFERENCE}/Homo_sapiens_assembly38.fasta \
	--runDir=${OUT}/manta \
	--bam=${BAM}/${SAMPLE1}_recal.bam \
	--bam=${BAM}/${SAMPLE2}_recal.bam \

${OUT}/manta/runWorkflow.py -j 32
