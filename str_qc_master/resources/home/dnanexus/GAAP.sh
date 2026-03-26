#!/bin/bash

filepath=$1
GAAP_threshold=$2
#remove GAAP
bcftools query -f '[%AAP\t]\n' $filepath | sed 's/\\t$//' > extracted_AAP.tsv
bcftools query -l $filepath > samples_all.txt
awk '                                                    
{
    for (i=1; i<=NF; i++) {
        sum[i] += $i;
        count[i]++;
    }
}
END {
    for (i=1; i<=NF; i++) {
        avg[i] = sum[i] / count[i];
        print avg[i];
    }
}' extracted_AAP.tsv > averages.txt
paste samples_all.txt averages.txt > sample_averages.txt
awk -v threshold="$GAAP_threshold" '$2 < threshold {print $1}' sample_averages.txt > low_GAAP.txt

