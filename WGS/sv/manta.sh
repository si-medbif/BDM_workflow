#! /bin/bash

set -e
set -u
set -o pipefail

MANTA=/home/harald/manta-1.6.0.centos6_x86_64/bin
REFERENCE=/gnome/genome_database/gatk_bundle/hg38bundle
TUMOR_BAM=$1
NORMAL_BAM=$2
TUMOR=$3
NORMAL=$4

python2 ${MANTA}/configManta.py \
	--normalBam=${NORMAL_BAM}/${NORMAL}/BAM/${NORMAL}_recal.bam \
	--tumorBam=${TUMOR_BAM}/${TUMOR}/BAM/${TUMOR}_recal.bam \
	--referenceFasta=${REFERENCE}/Homo_sapiens_assembly38.fasta \
	--runDir=${TUMOR_BAM}/${TUMOR}/manta/

${TUMOR_BAM}/${TUMOR}/manta/runWorkflow.py -j 8
