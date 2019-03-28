#!/bin/bash
set -e
set -u
set -o pipefail


accessionsfile=${1}
outdir=${2}

cat ${accessionsfile} | cut -f1 | epost -db nuccore -format acc | efetch -format docsum | xtract -pattern DocumentSummary -element AccessionVersion Topology Slen Title Completeness | awk -F "\t" '{ if (NF<5) { print $0"\t""not set" } else { print $0 }}' | sort -k3,3n | tee >(awk -F "\t" '{ if ($5 !~ /complete/) { print $0 }}' > ${outdir}/incompleteaccessions.tsv) | awk -F "\t" '{ if ($5 ~ /complete/) { print $0 }}' > ${outdir}/accessions.tsv
