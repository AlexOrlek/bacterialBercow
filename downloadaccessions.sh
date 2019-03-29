#!/bin/bash
set -e
set -u
set -o pipefail

datepresent=${1}
if [ "$datepresent" == "datepresent" ]; then
    taxonomyterm=${2}
    dateterm=${3}
    outdir=${4}
else
    taxonomyterm=${2}
    outdir=${3}
fi

date=$(date)
echo ${date} > ${outdir}/downloaddate.txt

#initial retrieval of accessions and titles, which will then be text-mined to decide which genbank files to retrieve
#no length filters applied at this stage

if [ "$datepresent" == "datepresent" ]; then
  esearch -db nuccore -query "${taxonomyterm} AND ${dateterm} AND biomol_genomic[PROP] AND plasmid[filter]" | efilter -source refseq | efetch -format docsum | xtract -pattern DocumentSummary -element AccessionVersion Topology Slen Title Completeness | awk -F "\t" '{ if (NF<5) { print $0"\t""not set" } else { print $0 }}' | sort -k3,3n | tee >(awk -F "\t" '{ if ($5 !~ /complete/) { print $0 }}' > ${1}/${downloaddate}/incompleteaccessions.tsv) | awk -F "\t" '{ if ($5 ~ /complete/) { print $0 }}' > ${outdir}/accessions.tsv
else
  esearch -db nuccore -query "${taxonomyterm} AND biomol_genomic[PROP] AND plasmid[filter]" | efilter -source refseq | efetch -format docsum | xtract -pattern DocumentSummary -element AccessionVersion Topology Slen Title Completeness | awk -F "\t" '{ if (NF<5) { print $0"\t""not set" } else { print $0 }}' | sort -k3,3n | tee >(awk -F "\t" '{ if ($5 !~ /complete/) { print $0 }}' > ${outdir}/incompleteaccessions.tsv) | awk -F "\t" '{ if ($5 ~ /complete/) { print $0 }}' > ${outdir}/accessions.tsv
fi

echo "finished downloading accessions"




#OLD CODE

#elif [ ${1} == 'bacteria' ]; then
#    esearch -db nuccore -query "bacteria[porgn:__txid2] AND biomol_genomic[PROP] AND plasmid[filter]" | efilter -source refseq | efetch -format docsum | xtract -pattern DocumentSummary -element AccessionVersion Topology Slen Title Completeness | awk '$5 ~ /complete/ { print $0 }' | sort -k3,3n > ${1}/${downloaddate}/accessions.tsv

#awk '$2 ~ /circular/ { print $0 }'  #may want to filter using length once sensible filters established + may want to filter by completeness


#OLD
# #initial retrieval of accessions and titles, which will then be text mined to decide which genbank files to retrieve (can also apply the inclusion criteria of complete sequence/plasmid/genome)
# if [ ${1} == 'ncbigproteobac' ]; then
#     esearch -db nuccore -query "g-proteobacteria[porgn] AND biomol_genomic[PROP] AND 1000:2000000[SLEN] AND plasmid[filter]" | efilter -source refseq | efetch -format docsum | xtract -pattern DocumentSummary -element AccessionVersion Topology Slen Title | awk '$2 ~ /circular/ { print $0 }' | sort -k3,3n > ${1}/${downloaddate}/accessions.tsv
# else
#     echo 'search not recognised'
# fi


#OLD - applying exclusion criteria in the initial search is uneccessary - just apply on title with python
# if [ ${1} == 'ncbigproteobac' ]; then
#     esearch -db nuccore -query "g-proteobacteria[porgn] AND biomol_genomic[PROP] AND 2000:2000000[SLEN] AND plasmid[filter] NOT complete cds[Title] NOT gene[Title] NOT genes[Title] NOT contig[Title] NOT scaffold[Title] NOT whole genome map[Title] NOT partial sequence[Title] NOT partial plasmid[Title] NOT locus[Title] NOT region[Title] NOT fragment[Title] NOT integron[Title] NOT transposon[Title] NOT insertion sequence[Title] NOT insertion element[Title] NOT phage[Title] NOT operon[Title]" | efilter -source refseq | efetch -format docsum | xtract -pattern DocumentSummary -element AccessionVersion Topology Slen Title | awl '$2 ~ /circular/ { print $0 }' > ${1}/${1}accessions_${downloaddate}.tsv
# else
#     echo 'search not recognised'
# fi


