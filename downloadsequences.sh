#!/bin/bash
set -e
set -u
set -o pipefail


batchsize=${1}
emailaddress=${2}
outdir=${3}

accessions=($(cut -f 1 "${outdir}/accessions_filtered.tsv" | sed '1d')) #first column with header removed

> ${outdir}/accessions_filtered.fa
echo -e 'Accession\tBioSample' > ${outdir}/accessions_filtered_biosamples.tsv
echo -e 'Accession\tFirst\tLast' > ${outdir}/accessions_filtered_metadata.tsv

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
	#echo "$chunkedaccessionsinput"	
	echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element AccessionVersion BioSample >> ${outdir}/accessions_filtered_biosamples.tsv
	echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | elink -target biosample | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element Accession First Last >> ${outdir}/accessions_filtered_metadata.tsv 
	echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | efetch -format fasta >> ${outdir}/accessions_filtered.fa
	break
    else
        echo $i
        chunkedaccessions=${accessions[@]:$i:$chunklen} #slice accessions array                                                           
        chunkedaccessionsinput=$(echo $chunkedaccessions | sed 's/ /\n/g')  #converting array to data column to use as epost input
	#echo "$chunkedaccessionsinput"
	echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element AccessionVersion BioSample >> ${outdir}/accessions_filtered_biosamples.tsv
	echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | elink -target biosample | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element Accession First Last >> ${outdir}/accessions_filtered_metadata.tsv 
	echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | efetch -format fasta >> ${outdir}/accessions_filtered.fa
        sleep 1
    fi
done

echo 'finished sequence download'


#OLD CODE - don't need both person's name and owner?
#	echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | elink -target biosample | efetch -format docsum | tee >(xtract -pattern DocumentSummary -element Accession First Last >> temp1) | xtract -pattern DocumentSummary -block Owner -element Name >>temp2
#OLD CODE - BEFORE USING -DEF ARGUMENT
# for i in $(eval echo {0..$len..$chunklen})
# do
#     sum=$(( ($i + $chunklen) + 1 ))
#     if [ $i -eq $len ]; then
#         break
#     elif [ $sum -eq $len ]; then
#         echo $i
#         chunklen=$(( $chunklen + 1 ))
#         chunkedaccessions=${accessions[@]:$i:$chunklen} #slice accessions array
# 	chunkedaccessionsinput=$(echo $chunkedaccessions | sed 's/ /\n/g')  #converting array to data column to use as epost input
# 	#echo "$chunkedaccessionsinput"
# 	echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | efetch -format docsum | xtract -pattern DocumentSummary -element AccessionVersion BioSample | awk -F "\t" '{ if (NF<2) { print $0"\t""not set" } else { print $0 }}' >> ${outdir}/accessions_filtered_biosamples.tsv
# 	echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | efetch -format fasta >> ${outdir}/accessions_filtered.fa
# 	break
#     else
#         echo $i
#         chunkedaccessions=${accessions[@]:$i:$chunklen} #slice accessions array                                                           
#         chunkedaccessionsinput=$(echo $chunkedaccessions | sed 's/ /\n/g')  #converting array to data column to use as epost input
# 	#echo "$chunkedaccessionsinput"
# 	echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | efetch -format docsum | xtract -pattern DocumentSummary -element AccessionVersion BioSample | awk -F "\t" '{ if (NF<2) { print $0"\t""not set" } else { print $0 }}' >> ${outdir}/accessions_filtered_biosamples.tsv
# 	echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | efetch -format fasta >> ${outdir}/accessions_filtered.fa
#         sleep 1
#     fi
# done

# echo 'finished sequence download'


###OLD CODE

#this fails
#echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | tee >(efetch -format docsum | xtract -pattern DocumentSummary -element AccessionVersion BioSample | awk -F "\t" '{ if (NF<2) { print $0"\t""not set" } else { print $0 }}' >> ${outdir}/accessions_filtered_biosamples.tsv) | efetch -format fasta >> ${outdir}/accessions_filtered.fa


###OLD CODE - retrieving full gb file
#cat ./${1}/${month}-${day}-${year}/accessions_filtered.tsv | cut -f1 | epost -db nuccore -format acc | efetch -format gbwithparts > ./${1}/${month}-${day}-${year}/accessions_filtered.gb

#OLD CODE - epost doesnt seem to work with array of accession ids
# accessions=($(cut -f1 ./${1}/${1}accessions_filtered_${month}-${day}-${year}.tsv))
# len=${#accessions[@]}
# echo $len
# testaccessions=(${accessions[@]:1:10})
# echo ${testaccessions[@]}
