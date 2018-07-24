#!/bin/bash

# Example file for general protocol workflow on one machine.
# Tasks should be split up over multiple machines according
# as described in the README.

# Generate data.
sh bin/generate_batches_pw.sh
sh bin/generate_batches.sh

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
./bin/MaskDnn 0 ../par/test.par.0.txt &
./bin/MaskDnn 1 ../par/test.par.1.txt &
./bin/MaskDnn 2 ../par/test.par.2.txt &
./bin/MaskDnn 3 ../par/test.par.3.txt ../../data/batch/ &
wait

# Run protocol.
./bin/DnnClient 0 ../par/test.par.0.txt &
./bin/DnnClient 1 ../par/test.par.1.txt &
./bin/DnnClient 2 ../par/test.par.2.txt &
wait
