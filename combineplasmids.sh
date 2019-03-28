#!/bin/bash
set -e
set -u
set -o pipefail

#use fastas in output_plasmids directories for inhouse data (any contig that is plasmid or putative_plasmid; unlike output_characterisedplasmids whihc is plasmids from characterised samples); contactenate


cat ${1}/accessions_final.fa > ${1}/accessions_final_inhouse.fa
cat '/well/bag/orlek/3rd_task/data/mankpc/fastqs_longread/finalgenomes/output_plasmids/allplasmids.fa' >> ${1}/accessions_final_inhouse.fa
cat '/well/bag/orlek/3rd_task/data/apha/fastqs_longread/finalgenomes/output_plasmids/allplasmids.fa' >> ${1}/accessions_final_inhouse.fa
