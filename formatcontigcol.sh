#!/bin/bash
set -e
set -u
set -o pipefail


#bash formatcontigcol filenamefiltered writefile_contigs
paste <(cat $1 | cut -f1 | cut -f1 -d'|') <(cat $1) > $2
