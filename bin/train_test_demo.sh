head -n150 data/X.txt > demo_data/batch_pw/Xtrain
head -n150 data/y.txt > demo_data/batch_pw/ytrain
tail -n+151 data/X.txt > demo_data/batch_pw/Xtest
tail -n+151 data/y.txt > demo_data/batch_pw/ytest
echo train > demo_data/batch_pw/train_suffixes.txt
echo test > demo_data/batch_pw/test_suffixes.txt
