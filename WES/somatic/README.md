#NOTE:
This pipeline uses several resources to help call somatic variants:

A region file showing where the genome has been sequenced.
	agilent_sureselect_V5_UTR_hg38_clean.bed
Panel-of-Normal (pon) with known sequencing artifacts.
	1000g_pon.hg38.vcf.gz
Germline SNPs and their allele frequency.
	af-only-gnomad.hg38.vcf.gz
	small_exac_common_3.hg38.vcf.gz
