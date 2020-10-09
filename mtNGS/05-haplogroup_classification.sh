#!/bin/bash

#Classify haplogroup
for PREFIX in $(cat /turtle/BSL/mon/sample.id); do java -jar ~/haplogrep-2.1.20/haplogrep-2.1.20.jar --hits 5 --format vcf --in /turtle/BSL/mon/variants/$PREFIX.filtered.vcf.gz --out /turtle/BSL/mon/haplogroup/$PREFIX.haplogroup; done

#Change the sample name from NC_012920.1 to the appropriate ID
for PREFIX in $(cat /turtle/BSL/mon/sample.id); do sed -i "s/NC_012920\.1/$PREFIX/" /turtle/BSL/mon/consensus/$PREFIX.consensus.fasta; done
