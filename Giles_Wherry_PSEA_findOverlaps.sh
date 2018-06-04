#!/bin/bash

for i in *_peaks.bed;
do
 filename="${i%_*.*}"
 echo "${filename}"
 bedtools intersect -a unionPeakList_all.bed -b $i -wb -wb -F 0.25 > ${filename}_overlap.txt
done
