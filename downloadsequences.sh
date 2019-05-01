#!/bin/bash
set -e
set -u
set -o pipefail


batchsize=${1}
emailaddress=${2}
outdir=${3}

accessions=($(cut -f1 "${outdir}/accessions_filtered.tsv" | sed '1d')) #first column with header removed

> ${outdir}/accessions_filtered.fa
echo -e 'Accession\tBioSample\tBioProject' > ${outdir}/accessions_filtered_dblinks.tsv

econtact -email ${emailaddress} -tool plasmiddownload


len=${#accessions[@]}
chunklen=${batchsize}

#download sequences as well as dblinks (biosample,bioproject)
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
	echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element AccessionVersion BioSample ProjectId >> ${outdir}/accessions_filtered_dblinks.tsv
	sleep 1
	#echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | elink -target biosample | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -first -VAR1 Accession -VAR2 First -VAR3 Last -block Owner -def "-" -sep " " -element "&VAR1" "&VAR2","&VAR3" -first Name >> ${outdir}/accessions_filtered_metadata.tsv
	#sleep 1
	echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | efetch -format fasta >> ${outdir}/accessions_filtered.fa
	break
    else
        echo $i
        chunkedaccessions=${accessions[@]:$i:$chunklen} #slice accessions array                                                           
        chunkedaccessionsinput=$(echo $chunkedaccessions | sed 's/ /\n/g')  #converting array to data column to use as epost input
	#echo "$chunkedaccessionsinput"
	echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element AccessionVersion BioSample ProjectId >> ${outdir}/accessions_filtered_dblinks.tsv
	sleep 1
	#echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | elink -target biosample | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -first -VAR1 Accession -VAR2 First -VAR3 Last -block Owner -def "-" -sep " " -element "&VAR1" "&VAR2","&VAR3" -first Name >> ${outdir}/accessions_filtered_metadata.tsv
	#sleep 1
	echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | efetch -format fasta >> ${outdir}/accessions_filtered.fa
        sleep 1
    fi
done


echo 'finished sequence download, and dblink retrieval'



#get biosample metadata

echo -e 'Accession\tName\tOwner' > ${outdir}/accessions_filtered_metadata.tsv

biosampleaccessions=($(cut -f2 "${outdir}/accessions_filtered_dblinks.tsv" | sed '1d' | awk '$1 != "-"' | sort | uniq))

len=${#biosampleaccessions[@]}

for i in $(eval echo {0..$len..$chunklen})
do
    sum=$(( ($i + $chunklen) + 1 ))
    if [ $i -eq $len ]; then
        break
    elif [ $sum -eq $len ]; then
        echo $i
        chunklen=$(( $chunklen + 1 ))
        chunkedaccessions=${biosampleaccessions[@]:$i:$chunklen} #slice accessions array
        chunkedaccessionsinput=$(echo $chunkedaccessions | sed 's/ /\n/g')  #converting array to data column to use as epost input                                                                          
        #echo "$chunkedaccessionsinput"                                                                                                                                                                     
        echo "$chunkedaccessionsinput" | epost -db biosample -format acc | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -first -VAR1 Accession -VAR2 First -VAR3 Last -block Owner -def "-" -sep " " -element "&VAR1" "&VAR2","&VAR3" -first Name >> ${outdir}/accessions_filtered_metadata.tsv
        break
    else
        echo $i
        chunkedaccessions=${biosampleaccessions[@]:$i:$chunklen} #slice accessions array
        chunkedaccessionsinput=$(echo $chunkedaccessions | sed 's/ /\n/g')  #converting array to data column to use as epost input                                                                          
        #echo "$chunkedaccessionsinput"                                                                                                                    
        echo "$chunkedaccessionsinput" | epost -db biosample -format acc | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -first -VAR1 Accession -VAR2 First -VAR3 Last -block Owner -def "-" -sep " " -element "&VAR1" "&VAR2","&VAR3" -first Name >> ${outdir}/accessions_filtered_metadata.tsv
        sleep 1
    fi
done

echo 'finished biosample metadata download'


#for refseq accessions get the cognate genbank accession (because I'll be using genbank bioproject ids)

refseqaccessions=()
for accession in ${accessions[@]}; do
    if [[ $accession == *_* ]]; then
	refseqaccessions+=($accession)
    fi
done

echo -e 'RefseqAccession\tCognateGenbankAccession' > ${outdir}/accessions_filtered_refseq_gb.tsv

if [ ${#refseqaccessions[@]} -gt 0 ]; then
    len=${#refseqaccessions[@]}
    chunklen=${batchsize}
    #get refseq/genbank cognates
    for i in $(eval echo {0..$len..$chunklen})
    do
	sum=$(( ($i + $chunklen) + 1 ))
	if [ $i -eq $len ]; then
            break
	elif [ $sum -eq $len ]; then
            echo $i
            chunklen=$(( $chunklen + 1 ))
            chunkedaccessions=${refseqaccessions[@]:$i:$chunklen} #slice accessions array
	    chunkedaccessionsinput=$(echo $chunkedaccessions | sed 's/ /\n/g')  #converting array to data column to use as epost input
	    echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element AccessionVersion AssemblyAcc >> ${outdir}/accessions_filtered_refseq_gb.tsv
	    break
	else
            echo $i
            chunkedaccessions=${refseqaccessions[@]:$i:$chunklen} #slice accessions array                                                           
            chunkedaccessionsinput=$(echo $chunkedaccessions | sed 's/ /\n/g')  #converting array to data column to use as epost input
	    echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element AccessionVersion AssemblyAcc >> ${outdir}/accessions_filtered_refseq_gb.tsv
	    sleep 1
	fi
    done
fi


echo 'finished retrieving cognate genbank accession ids for refseq accessions'


#retrieve any additional bioproject ids for cognate genbank accessions included in accessions_filtered_refseq_gb which are not included in accessions_filtered and hence not included in accessions_filtered_dblinks; this may occur for the following reasons:
#1 "refseq_genbank" accessions were retrived when querying NCBI but annotation data differs between corresponding refseq and genbank accessions such that the refseq accessions was included but the genbank was excluded e.g. A14660 is not annotated as "complete" whereas refseq cognate NC_025165.1 is. 
#2 only "refseq" accessions retrieved when querying NCBI or custom --accessions with missing cognate genbank accessions
#3 there can be cases where "refseq_genbank" query is used but a cognate genbank is not represented in accessions_filtered, due to removal from NCBI - in this case using such an accession as a query could lead to error (if there are no valid accessions in the query), hence "|| true" prevents for loop from breaking due to such an error  !but by ussing accession rather than accessions.version this shuldn't be a problem

#diff <(cut -f2 "${outdir}/accessions_filtered_refseq_gb.tsv" | sed '1d' | sort) <(cut -f1 "${outdir}/accessions_filtered.tsv" | sed '1d' | cut -f1 -d'.' | sort) | grep '<' | cut -f2 -d' ' > ${outdir}/accessions_filtered_missinggb.txt  #if there is no difference, accessions_filtered_missinggb.txt will be a blank file; missinggbaccessions will be empty array; len will be 0

#missinggbaccessions=($(cut -f1 "${outdir}/accessions_filtered_missinggb.txt"))

Array1=($(cut -f2 "${outdir}/accessions_filtered_refseq_gb.tsv" | sed '1d' | sort)) #cognate genbank accessions
Array2=($(cut -f1 "${outdir}/accessions_filtered.tsv" | sed '1d' | cut -f1 -d'.' | sort)) #filtered accessions
missinggbaccessions=(`echo ${Array1[@]} ${Array2[@]} ${Array2[@]} | tr ' ' '\n' | sort | uniq -u`)

len=${#missinggbaccessions[@]}

echo -e 'Accession\tBioProject' > ${outdir}/accessions_filtered_missinggb.tsv

if [[ ${#refseqaccessions[@]} -gt 0 && $len -gt 0 ]]; then  #if there are cognate genbank accessions and at least some of these are missing from accessions_filtered...
    #download additional bioprojects for missing genbank accessions
    for i in $(eval echo {0..$len..$chunklen})
    do
	sum=$(( ($i + $chunklen) + 1 ))
	if [ $i -eq $len ]; then
	    break
	elif [ $sum -eq $len ]; then
	    echo $i
	    chunklen=$(( $chunklen + 1 ))
	    chunkedaccessions=${missinggbaccessions[@]:$i:$chunklen} #slice accessions array
	    chunkedaccessionsinput=$(echo $chunkedaccessions | sed 's/ /\n/g')  #converting array to data column to use as epost input
	    #echo "$chunkedaccessionsinput"	
	    echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element AccessionVersion ProjectId >> ${outdir}/accessions_filtered_missinggb.tsv
	    break
	else
	    echo $i
	    chunkedaccessions=${missinggbaccessions[@]:$i:$chunklen} #slice accessions array                                                           
	    chunkedaccessionsinput=$(echo $chunkedaccessions | sed 's/ /\n/g')  #converting array to data column to use as epost input
	    #echo "$chunkedaccessionsinput"
	    echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element AccessionVersion ProjectId >> ${outdir}/accessions_filtered_missinggb.tsv
	    sleep 1
	fi
    done

    echo 'finished retrieving additional bioproject ids'
fi





#OLD CODE
        #echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | elink -target biosample | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element Accession First Last >> ${outdir}/accessions_filtered_metadata.tsv
	#echo "$chunkedaccessionsinput" | epost -db nuccore -format acc | elink -target biosample | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element Accession First Last >> ${outdir}/accessions_filtered_metadata.tsv

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
