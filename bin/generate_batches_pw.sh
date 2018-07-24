#!/bin/bash

zcat data/stitch/9606.protein_chemical.links.v5.0.tsv.gz | \
    tail -n+2 | \
    awk '$3 >= 400' | \
    sed 's/CID[sm]/CIDs/' | \
    sort -u -k1,1 -k2,2 | \
    shuf \
         > data/stitch/large_interactions_uniq.txt

python bin/generate_data_pw.py \
       data/stitch/large_interactions_uniq.txt \
       data/stitch/filtered_chem.txt \
       data/stitch/ensp_to_bitvec.txt

mkdir -p data/batch_pw

split -l 20000 data/X.txt data/batch_pw/X
split -l 20000 data/y.txt data/batch_pw/y
