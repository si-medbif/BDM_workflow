#!/bin/bash

set -e

FOLDERNAME=/gnome/tmp

for file in ${FOLDERNAME}/*_1.fq.gz
do
	f=$(echo "${file##*/}");
	filename=$(echo $f| cut  -d'_' -f 1,2);
	./01_align.sh "${filename}"
done

