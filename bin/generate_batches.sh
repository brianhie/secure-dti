#!/bin/bash

zcat data/stitch/9606.protein_chemical.links.v5.0.tsv.gz | \
    tail -n+2 | \
    awk '$3 >= 400' | \
    shuf \
        > "$1"data/stitch/large_interactions.txt

cat data/stitch/large_interactions.txt | \
    sed 's/CID[sm]/CIDs/' | \
    sort -u -k1,1 -k2,2 | \
    shuf \
        > "$1"data/stitch/large_interactions_uniq.txt

python bin/generate_data.py \
       "$1"data/stitch/large_interactions_uniq.txt \
       "$1"data/stitch/filtered_chem.txt \
       "$1"data/stitch/ensp_to_bitvec.txt

mkdir -p "$1"data/batch

split -l 20000 "$1"data/X.txt "$1"data/batch/X
split -l 20000 "$1"data/y.txt "$1"data/batch/y
