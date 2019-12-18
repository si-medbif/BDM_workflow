#!/bin/bash

# Preparation: Add HLA-alleles to the line "HLA="

# Usage: 06_EpitopePrediction.sh <output-folder> <sample-name> <HLA-alleles>
# HLA-alleles are provided as a comma-separated list of 4-digit alleles: 
#   HLA-A*24:02,HLA-B*54:01,HLA-C*03:02,DQA1*01:01,DQB1*05:03,DRB1*03:01
# There is no limit to the number of alleles that can be provided,
# although runtime will increase.

dir_Output=$1
sample=$2
HLA=$3

mkdir -p ${dir_Output}/${sample}/pvac

docker run --rm -v ${dir_Output}:/Output \
       -v ${dir_Output}/${sample}/pvac:/output_dir \
       griffithlab/pvactools \
	pvacseq run \
       /Output/${sample}/VCF/${sample}_m2_vep_filtered.vcf \
       -t 4 \
       --iedb-install-directory /opt/iedb \
	${sample} \
	${HLA} \
	MHCflurry MHCnuggetsI MHCnuggetsII NNalign NetMHC PickPocket SMM SMMPMBEC SMMalign \
	/output_dir/ \
       -e 8,9,10,11,12
