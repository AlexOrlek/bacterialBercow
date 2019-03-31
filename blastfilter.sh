#!/bin/bash
set -e
set -u
set -o pipefail

#$1 is longerqueries boolean; $2 is formatcontigcol boolean $3 is blastoutputpath; $4 is sortedfilepath; $5 is pidthresh; $6 is coveragethresh

if [ $1 == 'True' ]; then
    if [ $2 == 'True' ]; then
	echo 'long format'
        cat ${3} | awk '{ print $0 "\t" ($3/100)*(($4 - $6)/$16)}' | awk -F $"\t" 'BEGIN {OFS = FS} {$14 = (($4 - $6)/$16); print $0}' | awk -v var1="$5" -v var2="$6" '$3 > var1 && $14 > var2 { print $0 }' | awk '{split($0,a,"|"); print a[1] "\t" $0}' | sort -t$'\t' -k1,1V -k18,18nr > ${4}
    else
	echo 'long dont format'
	cat ${3} | awk '{ print $0 "\t" ($3/100)*(($4 - $6)/$16)}' | awk -F $"\t" 'BEGIN {OFS = FS} {$14 = (($4 - $6)/$16); print $0}' | awk -v var1="$5" -v var2="$6" '$3 > var1 && $14 > var2 { print $0 }' | sort -t$'\t' -k1,1V -k17,17nr > ${4}
    fi
else
    if [ $2 == 'True' ]; then
	echo 'short format'
	cat ${1} | awk '{ print $0 "\t" ($3/100)*$14}' | awk -v var1="$5" -v var2="$6" '$3 > var1 && $14 > var2 { print $0 }' | awk '{split($0,a,"|"); print a[1] "\t" $0}' | sort -t$'\t' -k1,1V -k18,18nr > ${4}
    else
	echo 'short dont format'
	cat ${1} | awk '{ print $0 "\t" ($3/100)*$14}' | awk -v var1="$5" -v var2="$6" '$3 > var1 && $14 > var2 { print $0 }' | sort -t$'\t' -k1,1V -k17,17nr > ${4}
    fi
fi
