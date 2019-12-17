#!/bin/bash
# Add gene expression data to the vcf files
vcf-expression-annotator EFDFC8-T_m2_vep_filtered.matched.vcf EFDFC8-T_expression2.txt custom gene --id-column gene_id --expression-column rna_count_EFDFC8 -s EFDFC8-T
vcf-expression-annotator C53C69-T_m2_vep_filtered.matched.vcf C53C69-T_expression2.txt custom gene --id-column gene_id --expression-column rna_count_C53C69 -s C53C69-T

# Sample EFDFC8
# Add tumour DNA read counts to the vcf files
vcf-readcount-annotator EFDFC8-T_m2_vep_filtered.matched.gx.vcf coverage/DNA/EFDFC8-T_bam_readcount_snv.tsv DNA -s EFDFC8-T -t snv -o EFDFC8-T_1.vcf
vcf-readcount-annotator EFDFC8-T_1.vcf coverage/DNA/EFDFC8-T_bam_readcount_indel.tsv DNA -s EFDFC8-T -t indel -o EFDFC8-T_2.vcf

# Add normal DNA read counts to the vcf files
vcf-readcount-annotator EFDFC8-T_2.vcf coverage/DNA/EFDFC8-N_bam_readcount_snv.tsv DNA -s EFDFC8-N -t snv -o EFDFC8-T_3.vcf
vcf-readcount-annotator EFDFC8-T_3.vcf coverage/DNA/EFDFC8-N_bam_readcount_indel.tsv DNA -s EFDFC8-N -t indel -o EFDFC8-T_4.vcf

# Add RNA read counts to the vcf files
vcf-readcount-annotator EFDFC8-T_4.vcf coverage/RNA/EFDFC8-T_bam_readcount_snv.tsv RNA -s EFDFC8-T -t snv -o EFDFC8-T_5.vcf
vcf-readcount-annotator EFDFC8-T_5.vcf coverage/RNA/EFDFC8-T_bam_readcount_indel.tsv RNA -s EFDFC8-T -t indel -o EFDFC8-T_m2_vep_filtered.matched.gx.cov.vcf
#rm EFDFC8-T_?.vcf

# Sample C53C69
# Add tumour DNA read counts to the vcf files
vcf-readcount-annotator C53C69-T_m2_vep_filtered.matched.gx.vcf coverage/DNA/C53C69-T_bam_readcount_snv.tsv DNA -s C53C69-T -t snv -o C53C69-T_1.vcf
vcf-readcount-annotator C53C69-T_1.vcf coverage/DNA/C53C69-T_bam_readcount_indel.tsv DNA -s C53C69-T -t indel -o C53C69-T_2.vcf

# Add normal DNA read counts to the vcf files
vcf-readcount-annotator C53C69-T_2.vcf coverage/DNA/C53C69-N_bam_readcount_snv.tsv DNA -s C53C69-N -t snv -o C53C69-T_3.vcf
vcf-readcount-annotator C53C69-T_3.vcf coverage/DNA/C53C69-N_bam_readcount_indel.tsv DNA -s C53C69-N -t indel -o C53C69-T_4.vcf

# Add RNA read counts to the vcf files
vcf-readcount-annotator C53C69-T_4.vcf coverage/RNA/C53C69-T_bam_readcount_snv.tsv RNA -s C53C69-T -t snv -o C53C69-T_5.vcf
vcf-readcount-annotator C53C69-T_5.vcf coverage/RNA/C53C69-T_bam_readcount_indel.tsv RNA -s C53C69-T -t indel -o C53C69-T_m2_vep_filtered.matched.gx.cov.vcf
#rm C53C69-T_?.vcf
