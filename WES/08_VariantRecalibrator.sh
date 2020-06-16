#!/bin/bash

dir_Output=$1
REFERENCE=/gnome/genome_database/gatk_bundle/hg38bundle

docker run --rm -v ${dir_Output}:/Output \
	-v ${REFERENCE}:/ref \
	broadinstitute/gatk:4.1.5.0 gatk --java-options "-Xmx64g" \
	VariantRecalibrator \
	-R /ref/Homo_sapiens_assembly38.fasta \
	-V /Output/project.raw.vcf.gz \
	--resource:hapmap,known=false,training=true,truth=true,prior=15.0 /ref/hapmap_3.3.hg38.vcf.gz \
	--resource:omni,known=false,training=true,truth=false,prior=12.0 /ref/1000G_omni2.5.hg38.vcf.gz \
	--resource:1000G,known=false,training=true,truth=false,prior=10.0 /ref/1000G_phase1.snps.high_confidence.hg38.vcf.gz \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /ref/dbsnp_146.hg38.vcf.gz \
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	-mode SNP \
	-O /Output/dengue.recal \
	--tranches-file /Output/project.tranches \
	--rscript-file /Output/project.plots.R

docker run --rm -v ${dir_Output}:/Output \
	-v ${REFERENCE}:/ref \
	broadinstitute/gatk:4.1.5.0 gatk --java-options "-Xmx64g" \
	ApplyVQSR \
	-R /ref/Homo_sapiens_assembly38.fasta \
	-V /Output/project.raw.vcf.gz \
	-O /Output/project.filtered.vcf.gz \
	--truth-sensitivity-filter-level 99.0 \
	--tranches-file /Output/project.tranches \
	--recal-file /Output/project.recal \
	-mode SNP

docker run --rm -v ${dir_Output}:/Output \
	-v ${REFERENCE}:/ref \
	broadinstitute/gatk:4.1.5.0 gatk --java-options "-Xmx64g" \
	VariantRecalibrator \
	-R /ref/Homo_sapiens_assembly38.fasta \
	-V /Output/project.raw.vcf.gz \
	--resource:mills,known=false,training=true,truth=true,prior=12.0 /ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz \
	--resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /ref/dbsnp_146.hg38.vcf.gz \
	-an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR \
	-mode INDEL \
	-O /Output/project.indel.recal \
	--tranches-file /Output/project.indel.tranches \
	--rscript-file /Output/project.indel.plots.R

docker run --rm -v ${dir_Output}:/Output \
	-v ${REFERENCE}:/ref \
	broadinstitute/gatk:4.1.5.0 gatk --java-options "-Xmx64g" \
	ApplyVQSR \
	-R /ref/Homo_sapiens_assembly38.fasta \
	-V /Output/project.raw.vcf.gz \
	-O /Output/project.filtered.indel.vcf.gz \
	--truth-sensitivity-filter-level 99.0 \
	--tranches-file /Output/project.indel.tranches \
	--recal-file /Output/project.indel.recal \
	-mode INDEL
	
