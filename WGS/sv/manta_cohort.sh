#! /bin/bash

set -e
set -u
set -o pipefail

MANTA=/home/harald/manta-1.6.0.centos6_x86_64/bin
REFERENCE=/gnome/genome_database/gatk_bundle/hg38bundle
BAM=/gnome/genomics_thailand/breast-20191107/Blood_WGS/neoantigens/BAM
OUT=/gnome/genomics_thailand/breast-20191107/Blood_WGS/neoantigens

python2 ${MANTA}/configManta.py \
	--referenceFasta=${REFERENCE}/Homo_sapiens_assembly38.fasta \
	--runDir=${OUT}/manta \
	--bam=${BAM}/BB-B0006-DNA_recal.bam \
	--bam=${BAM}/BB-B0012-DNA_recal.bam \
	--bam=${BAM}/BB-B0013-DNA_recal.bam \
	--bam=${BAM}/BB-B0030-DNA_recal.bam \
	--bam=${BAM}/BB-B0040-DNA_recal.bam \
	--bam=${BAM}/BB-B0066-DNA_recal.bam \
	--bam=${BAM}/BB-B0068-DNA_recal.bam \
	--bam=${BAM}/BB-B0069-DNA_recal.bam \
	--bam=${BAM}/BB-B0079-DNA_recal.bam \
	--bam=${BAM}/BB-B0082-DNA_recal.bam \
	--bam=${BAM}/BB-B0099-DNA_recal.bam \
	--bam=${BAM}/BB-B0102-DNA_recal.bam \
	--bam=${BAM}/BB-B0107-DNA_recal.bam \
	--bam=${BAM}/BB-B0110-DNA_recal.bam \
	--bam=${BAM}/BB-B0119-DNA_recal.bam \
	--bam=${BAM}/BB-B0126-DNA_recal.bam \
	--bam=${BAM}/BB-B0130-DNA_recal.bam \
	--bam=${BAM}/BB-B0136-DNA_recal.bam \
	--bam=${BAM}/BB-B0145-DNA_recal.bam \
	--bam=${BAM}/BB-B0158-DNA_recal.bam \
	--bam=${BAM}/BB-B0160-DNA_recal.bam \
	--bam=${BAM}/BB-B0168-DNA_recal.bam \
	--bam=${BAM}/BB-B0181-DNA_recal.bam \
	--bam=${BAM}/BB-B0184-DNA_recal.bam \
	--bam=${BAM}/BB-B0195-DNA_recal.bam \
	--bam=${BAM}/BB-B0225-DNA_recal.bam \
	--bam=${BAM}/BB-B0242-DNA_recal.bam \
	--bam=${BAM}/BB-B0272-DNA_recal.bam \
	--bam=${BAM}/BB-B0275-DNA_recal.bam \
	--bam=${BAM}/BB-B0277-DNA_recal.bam \
	--bam=${BAM}/BB-B0283-DNA_recal.bam \
	--bam=${BAM}/BB-B0287-DNA_recal.bam \
	--bam=${BAM}/BB-B0293-DNA_recal.bam \
	--bam=${BAM}/BB-B0305-DNA_recal.bam \
	--bam=${BAM}/BB-B0314-DNA_recal.bam \
	--bam=${BAM}/BB-B0367-DNA_recal.bam \
	--bam=${BAM}/BB-B0402-DNA_recal.bam \
	--bam=${BAM}/BB-B0410-DNA_recal.bam \
	--bam=${BAM}/BB-B0413-DNA_recal.bam \
	--bam=${BAM}/BB-B0443-DNA_recal.bam \
	--bam=${BAM}/BB-B0639-DNA_recal.bam \
	--bam=${BAM}/BB-B0860-DNA_recal.bam \
	--bam=${BAM}/BB-B0863-DNA_recal.bam \
	--bam=${BAM}/BB-B0892-DNA_recal.bam

${OUT}/manta/runWorkflow.py -j 32
