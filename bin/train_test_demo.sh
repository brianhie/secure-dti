head -n150 data/X.txt > data/batch_pw/Xtrain
head -n150 data/y.txt > data/batch_pw/ytrain
tail -n+151 data/X.txt > data/batch_pw/Xtest
tail -n+151 data/y.txt > data/batch_pw/ytest
echo train > data/batch_pw/train_suffixes.txt
echo test > data/batch_pw/test_suffixes.txt
