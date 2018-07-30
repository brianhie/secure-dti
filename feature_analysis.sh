#!/bin/bash

# This script runs a protocol that demonstrates how the collaborating entities
# can perform additional computation to better interpret the prediction
# results from our Secure DTI model, while maintaining the same
# privacy guarantees as the main protocol.
#
# Note that this script runs everything on a single machine; in practice
# the computation should be split among multiple computing entities.
#
# The output of this script reproduces Table S3 of our manuscript.
#

cd mpc/code/

# Run protocol.
./bin/AnalyzeFeature 0 ../par/test.par.0.txt &
./bin/AnalyzeFeature 1 ../par/test.par.1.txt &
./bin/AnalyzeFeature 2 ../par/test.par.2.txt &
wait
