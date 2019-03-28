#!/bin/bash
set -e
set -u
set -o pipefail


batchsize=${1}
emailaddress=${2}
outdir=${3}

accessions=($(cut -f 1 "${outdir}/accessions_filtered.tsv"))

> ${outdir}/accessions_filtered.fa

econtact -email ${emailaddress} -tool plasmiddownload

len=${#accessions[@]}
chunklen=${batchsize}


for i in $(eval echo {0..$len..$chunklen})
do
    sum=$(( ($i + $chunklen) + 1 ))
    if [ $i -eq $len ]; then
        break
    elif [ $sum -eq $len ]; then
        echo $i
        chunklen=$(( $chunklen + 1 ))
        chunkedaccessions=${accessions[@]:$i:$chunklen} #slice accessions array
	chunkedaccessionsinput=$(echo $chunkedaccessions | sed 's/ /\n/g')  #converting array to data column to use as epost input
	echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | efetch -format fasta >> ${outdir}/accessions_filtered.fa
	break
    else
        echo $i
        chunkedaccessions=${accessions[@]:$i:$chunklen} #slice accessions array                                                                                                                       
        chunkedaccessionsinput=$(echo $chunkedaccessions | sed 's/ /\n/g')  #converting array to data column to use as epost input                                                                    
        echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | efetch -format fasta >> ${outdir}/accessions_filtered.fa
        sleep 1
    fi
done

echo 'finished sequence download'



###OLD CODE - retrieving full gb file
#cat ./${1}/${month}-${day}-${year}/accessions_filtered.tsv | cut -f1 | epost -db nuccore -format acc | efetch -format gbwithparts > ./${1}/${month}-${day}-${year}/accessions_filtered.gb

#OLD CODE - epost doesnt seem to work with array of accession ids
# accessions=($(cut -f1 ./${1}/${1}accessions_filtered_${month}-${day}-${year}.tsv))
# len=${#accessions[@]}
# echo $len
# testaccessions=(${accessions[@]:1:10})
# echo ${testaccessions[@]}
