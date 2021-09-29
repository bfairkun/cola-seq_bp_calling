#!/usr/bin/env bash
if [ "$1" == "-h" ]; then
  echo "Usage: `basename $0` <BedFileIn.gz> <MergeDistance> <BedFileOut.gz>"
  echo ""
  echo "Merge features in input bedfile within <MergeDistance> into the single best from score column"
  exit 0
fi

# set -xe #Debug mode

InputFile=$1
Distance=$2
OutputFile=$3

zcat $InputFile | awk -F'\t' -v OFS='\t' '{$5=sprintf("%.5f",$5); print $0}' | bedtools merge -i - -d $Distance -s -o count,max,first  -c 4,5,6 | bedtools intersect -a <(zcat $InputFile | awk -F'\t' -v OFS='\t' '{$5=sprintf("%.5f",$5); print $0}') -b -  -wao -s | awk -F'\t' -v OFS='\t' '$5==$11 {print $1,$2,$3,$4,$5,$6 }' | gzip - > $OutputFile
