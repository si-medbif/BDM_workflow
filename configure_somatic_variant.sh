#!/bin/bash

PROJECT=$1

mkdir -p ${PROJECT}

cp DNA/alignment.sh ${PROJECT}/01_alignment.sh
cp DNA/mark_duplicate.sh ${PROJECT}/02_mark_duplicate.sh
cp DNA/recalibration.sh ${PROJECT}/03_recalibration.sh
cp DNA/somatic_variant_call.sh ${PROJECT}/04_somatic_variant_call.sh
cp DNA/annotate_vep.sh ${PROJECT}/08_annotate_vep.sh 
