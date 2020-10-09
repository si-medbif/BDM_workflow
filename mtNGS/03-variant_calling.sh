#!/bin/bash

## Call varints
for PREFIX in $(cat /turtle/BSL/mon/sample.id); do 
  bcftools mpileup -Ou -f /penguin/resources/rCRS/rCRS.fasta \
/turtle/BSL/mon/aligned/$PREFIX.readgroup.bam | \
bcftools call -mv -Oz -o /turtle/BSL/mon/variants/$PREFIX.vcf.gz \
tabix /turtle/BSL/mon/variants/$PREFIX.vcf.gz
done

## Normalize indels
for PREFIX in $(cat /turtle/BSL/mon/sample.id); do 
  bcftools norm -f /penguin/resources/rCRS/rCRS.fasta \
/turtle/BSL/mon/variants/$PREFIX.vcf.gz \
-Oz -o /turtle/BSL/mon/variants/$PREFIX.norm.vcf.gz
done

## Filter adjacent indels within 5 bp
for PREFIX in $(cat /turtle/BSL/mon/sample.id); do 
  bcftools filter --IndelGap 5 /turtle/BSL/mon/variants/$PREFIX.norm.vcf.gz \
-Oz -o /turtle/BSL/mon/variants/$PREFIX.filtered.vcf.gz
done

## Collect the total number of SNPs
for PREFIX in $(cat /turtle/BSL/mon/sample.id); do 
  bcftools stats /turtle/BSL/mon/variants/$PREFIX.filtered.vcf.gz | \
grep 'number of SNPs:' | \
awk '{print "Total Number of SNPs = " $6}' \
>> /turtle/BSL/mon/stats/$PREFIX.summary.metrics
done

## Collect the total number of indels
for PREFIX in $(cat /turtle/BSL/mon/sample.id); do 
  bcftools stats /turtle/BSL/mon/variants/$PREFIX.filtered.vcf.gz \
| grep 'number of indels:' \
| awk '{print "Total Number of Indels = " $6}' \
>> /turtle/BSL/mon/stats/$PREFIX.summary.metrics
done
