ls $1 | \
    grep -v masked | \
    grep X | \
    sed 's/X//' | \
    head -n65 \
         > $1/train_suffixes.txt

ls $1 | \
    grep -v masked | \
    grep X | \
    sed 's/X//' | \
    tail -n+66 \
         > $1/test_suffixes.txt
