This is a pipeline to detect germline variants in WES data.

## Suggested workflow

1. Check read quality
	* ./00_Fastqc.sh FASTQFOLDER SAMPLE
	* ./create_fastqc_report.py FASTQFOLDER SAMPLE
2. Align reads to reference
	* ./01_Alignment.sh FASTQFOLDER BAMFOLDER SAMPLE
	* ./02_MarkDuplicate.sh BAMFOLDER SAMPLE
	* ./03_Recalibration.sh BAMFOLDER SAMPLE
	* ./collect_bam_stats.sh BAMFOLDER SAMPLE
	* ./create_bam_report.py BAMFOLDER SAMPLE
3. Produce quality report for reads and alignment
	* ./create_final_report.py
4. Call variants
	* ./04_VariantCall.sh BAMFOLDER VCFFOLDER SAMPLE
	* ./05_CombineGVCFs.sh VCFFOLDER CHROM SAMPLE1 SAMPLE2
	* ./06_GenotypeGVCFs.sh VCFFOLDER CHROM
	* ./07_CombineVCF.sh VCFFOLDER
	* ./08_VariantRecalibrator.sh VCFFOLDER
	* ./09_Annotation.sh VCFFOLDER
