#!/bin/bash

# Example file for general protocol workflow on one machine.
# Tasks should be split up over multiple machines according
# to the description in the README.

# Download and unpack data.
wget http://secure-dti.csail.mit.edu/data.tar.gz
tar xvf data.tar.gz

# Generate data.
sh bin/generate_batches.sh # Randomly sample negative set (Secure DTI-A).
sh bin/generate_batches_pw.sh # Control for degree (Secure DTI-B).

cd mpc/

# Generate keys for secret sharing.
code/bin/GenerateKey key/P0_P1.key
code/bin/GenerateKey key/P0_P2.key
code/bin/GenerateKey key/P1_P2.key
code/bin/GenerateKey key/P1_P3.key
code/bin/GenerateKey key/P2_P3.key
code/bin/GenerateKey key/global.key

cd code/

# Mask data.
./bin/ShareData 0 ../par/test.par.0.txt &
./bin/ShareData 1 ../par/test.par.1.txt &
./bin/ShareData 2 ../par/test.par.2.txt &
./bin/ShareData 3 ../par/test.par.3.txt ../../data/batch/ &
wait

# Run protocol.
./bin/TrainSecureDTI 0 ../par/test.par.0.txt &
./bin/TrainSecureDTI 1 ../par/test.par.1.txt &
./bin/TrainSecureDTI 2 ../par/test.par.2.txt &
wait
