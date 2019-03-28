#!/bin/bash
set -e
set -u
set -o pipefail

outdir=${1}

> ${outdir}/rmlstloci.txt
mkdir -p ${outdir}/blastdbs

fastas=($(find ${outdir}/ -maxdepth 1 -mindepth 1 -type f -name '*.tfa' -printf '%f\n'))
for fasta in ${fastas[@]}; do
    locus=$(echo ${fasta} | cut -f1 -d'.')
    echo ${locus} >> rmlstloci.txt
    makeblastdb -dbtype nucl -in ${outdir}/${locus}.tfa -out ${outdir}/blastdbs/${locus}db
done
