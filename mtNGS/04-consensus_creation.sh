#!/bin/bash

#Create consensus sequences by applying VCF variants
for PREFIX in $(cat /turtle/BSL/mon/sample.id); do 
  bcftools consensus -e 'INFO/DP<30' -f /penguin/resources/rCRS/rCRS.fasta \
-o /turtle/BSL/mon/consensus/$PREFIX.consensus.fasta \
/turtle/BSL/mon/variants/$PREFIX.vcf.gz; 
done
#Ref: https://samtools.github.io/bcftools/howtos/consensus-sequence.html
