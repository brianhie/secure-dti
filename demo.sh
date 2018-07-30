#!/bin/bash

# Example file for general protocol workflow on one machine.
# Tasks should be split up over multiple machines according
# to the description in the README.

# Download and unpack data.
#wget http://secure-dti.csail.mit.edu/demo_data.tar.gz
#tar xvf data.tar.gz

# Generate data.
#sh bin/generate_batches_pw.sh
#
## Split into train and test.
#sh bin/train_test_demo.sh
#
cd mpc/code/
#
## Generate keys for secret sharing.
#./bin/GenerateKey ../key/P0_P1.key
#./bin/GenerateKey ../key/P0_P2.key
#./bin/GenerateKey ../key/P1_P2.key
#./bin/GenerateKey ../key/P1_P3.key
#./bin/GenerateKey ../key/P2_P3.key
#./bin/GenerateKey ../key/global.key
#
## Mask data.
#./bin/ShareData 0 ../par/demo.par.0.txt &
#./bin/ShareData 1 ../par/demo.par.1.txt &
#./bin/ShareData 2 ../par/demo.par.2.txt &
#./bin/ShareData 3 ../par/demo.par.3.txt ../../data/batch_pw/ &
#wait

# Run protocol.
./bin/TrainSecureDTI 0 ../par/demo.par.0.txt &
./bin/TrainSecureDTI 1 ../par/demo.par.1.txt &
./bin/TrainSecureDTI 2 ../par/demo.par.2.txt &
wait

python bin/evaluate_demo.py
