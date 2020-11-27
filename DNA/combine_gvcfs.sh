#!/bin/bash

# Add as many samples as are available
# Can be run pr. chromosome or for the whole genome
# cohort.sample_map:
#     sample1      /gvcf/sample1.vcf.gz
#     sample2      /gvcf/sample2.vcf.gz
#     sample3      /gvcf/sample3.vcf.gz

VCFFOLDER=$1
CHROM=$2
SAMPLES=$3

docker run --rm \
	-v ${VCFFOLDER}:/gvcf \
	-v $PWD:/here \
	broadinstitute/gatk:4.1.9.0 gatk --java-options "-Xmx64g" \
	GenomicsDBImport \
	--batch-size 50 \
	--sample-name-map /here/${SAMPLES} \
	--genomicsdb-workspace-path /gvcf/project_db.${CHROM} \
	--tmp-dir /gvcf/tmp \
	--reader-threads 5 \
	-L ${CHROM}

# Alternative way of running, with samples specified on the command line

#SAMPLE1=$3
#SAMPLE2=$4

#docker run --rm -v ${VCFFOLDER}:/Output \
#	broadinstitute/gatk:4.1.9.0 gatk --java-options "-Xmx64g" \
#	GenomicsDBImport \
#	-V /Output/${SAMPLE1}.g.vcf \
#	-V /Output/${SAMPLE2}.g.vcf \
#	--genomicsdb-workspace-path /Output/dengue_db.${CHROM} \
#	--tmp-dir=/Output/tmp \
#	-L ${CHROM}
