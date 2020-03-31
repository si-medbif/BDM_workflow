#!/bin/bash

# Add as many samples as are available
# Can be run pr. chromosome or for the whole genome
# cohort.sample_map:
#     sample1      /PATH/sample1.vcf.gz
#     sample2      /PATH/sample2.vcf.gz
#     sample3      /PATH/sample3.vcf.gz

dir_Output=$1
CHROM=$4

#SAMPLES=$2
#docker run --rm -v ${dir_Output}:/Output \
#    broadinstitute/gatk:4.1.5.0 gatk --java-options "-Xmx64g" GenomicsDBImport \
#    --batch-size 50 \
#    --sample-name-map ${SAMPLES} \
#    --genomicsdb-workspace-path /Output/project_db.${CHROM} \
#    --tmp-dir=/Output/tmp \
#    --reader-threads 5 \
#    -L ${CHROM}


SAMPLE1=$2
SAMPLE2=$3

docker run --rm -v ${dir_Output}:/Output \
	broadinstitute/gatk:4.1.5.0 gatk --java-options "-Xmx64g" GenomicsDBImport \
	-V /Output/${SAMPLE1}/${SAMPLE1}.g.vcf \
	-V /Output/${SAMPLE2}/${SAMPLE2}.g.vcf \
	--genomicsdb-workspace-path /Output/dengue_db.$1 \
	--tmp-dir=/Output/tmp \
	-L $1
