#!/bin/bash
set -e
set -u
set -o pipefail


accessionsfile=${1}
outdir=${2}

echo -e 'Accession\tTopology\tLength\tTitle\tCompleteness' > ${outdir}/incompleteaccessions.tsv
echo -e 'Accession\tTopology\tLength\tTitle\tCompleteness' > ${outdir}/accessions.tsv

cat ${accessionsfile} | cut -f1 | epost -db nuccore -format acc | efetch -format docsum | xtract -pattern DocumentSummary -def "-" -element AccessionVersion Topology Slen Title Completeness | sort -k3,3n | tee >(awk -F "\t" '{ if ($5 !~ /complete/) { print $0 }}' >> ${outdir}/incompleteaccessions.tsv) | awk -F "\t" '{ if ($5 ~ /complete/) { print $0 }}' >> ${outdir}/accessions.tsv
