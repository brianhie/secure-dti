#!/bin/bash

# Example file for general protocol workflow on one machine.
# Tasks should be split up over multiple machines according
# to the description in the README.

# Generate data.
sh bin/generate_batches_pw.sh demo_

# Split into train and test.
sh bin/train_test_demo.sh

cd mpc/code/

# Generate keys for secret sharing.
./bin/GenerateKey ../key/P0_P1.key
./bin/GenerateKey ../key/P0_P2.key
./bin/GenerateKey ../key/P1_P2.key
./bin/GenerateKey ../key/P1_P3.key
./bin/GenerateKey ../key/P2_P3.key
./bin/GenerateKey ../key/global.key

# Mask data.
./bin/ShareData 0 ../par/demo.par.0.txt &
./bin/ShareData 1 ../par/demo.par.1.txt &
./bin/ShareData 2 ../par/demo.par.2.txt &
./bin/ShareData 3 ../par/demo.par.3.txt ../../demo_data/batch_pw/ &
wait

# Wait for system ports to settle down.
sleep 100

# Run protocol.
./bin/TrainSecureDTI 0 ../par/demo.par.0.txt &
./bin/TrainSecureDTI 1 ../par/demo.par.1.txt &
./bin/TrainSecureDTI 2 ../par/demo.par.2.txt &
wait

cd ../..

# Evaluate neural network performance in plaintext.
python bin/evaluate_demo.py
