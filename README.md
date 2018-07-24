# Secure DTI

## This package contains the client software for an end-to-end
## multiparty computation (MPC) protocol for secure DTI prediction.

**Dependencies:**

clang++ compiler (3.9; https://clang.llvm.org/)
GMP library (6.1.2; https://gmplib.org/)
libssl-dev package (1.0.2g-1ubuntu11.2)
NTL package (10.3.0; http://www.shoup.net/ntl/)

**Notes on NTL:**

We made a few modifications to NTL for our random streams.
Copy the contents of `code/NTL_mod/` into the NTL source
directory as follows before compiling it:

Copy `ZZ.cpp` into `NTL_PACKAGE_DIR/src/`
Copy `ZZ.h` into `NTL_PACKAGE_DIR/include/NTL/`

In addition, we recommend setting `NTL_THREAD_BOOST=on`
during the configuration of NTL to enable thread-boosting.

**Compilation:**

First, update the paths in `code/Makefile`:

`CPP` points to the clang++ compiler executable.
`INCPATHS` points to the header files of installed libraries.
`LDPATH` contains the `.a` files of installed libraries.

To compile, run `make` inside the `code/` directory.

This will create three executables of interest:

`bin/GenerateKey`
`bin/MaskDnn`
`bin/DnnClient`

**How to run:**

Our MPC protocol consists of four entities: SP, CP0, CP1, and CP2.

Note we treat SP as a single entity that holds all input features
and labels, but it is straightforward to generalize this setup
to the collaborative scenario where multiple SPs securely share
their data with the CPs, or function as the CPs themselves.

An instance of the client program is created for each involved
party on different machines, where the ID of the corresponding
party is provided as an input argument: `0=CP0`, `1=CP1`, `2=CP2`,
and `3=SP`.

These multiple instances of the client will interact over the
network to jointly carry out the MPC protocol.

#### Step 0: Generate DTI Features and Labels

As described in the paper, we encode a drug-target interaction as a
binary vector which is fed into a neural network model that predicts
whether the corresponding drug-target pair is interactive (output
value of 1) or not (0).

The `bin/generate_batches.sh` helper script will take drug-target
interactions and map them to the corresponding bitvectors. The data
is split up into file batches of 20,000 features each.

#### Step 1: Setup Shared Random Keys

Secure communication channels needed for the overall protocol are:

CP0 <-> CP1, CP0 <-> CP2, CP1 <-> CP2, CP1 <-> SP, CP2 <-> SP

Use GenerateKey to obtain a secret key for each pair, which should
be named:

P0_P1.key, P0_P2.key, P1_P2.key, P1_P3.key, P2_P3.key

The syntax for running `GenerateKey` is as follows:

./GenerateKey out.key

In addition, generate `global.key` and share it with all parties.

We provide pre-generated keys in the `key/` directory. In practice,
each party should have the keys for only the channels
they are involved in.

#### Step 2: Setup Parameters

We provide example parameter settings in:

par/test.par.PARTY_ID.txt

For more information about each parameter, please consult `code/param.h`
and Supplementary Information of our publication.

For a test run, update the following parameters and leave the rest:

PORT_*
IP_ADDR_*

#### Step 3: Setup Input Data 

On the machine where the SP instance will be running, the data set
should be available in plaintext. The data should have been split up
into multiple file batches as by `bin/generate_batches.sh`. The
directory prefix of these file batches can be specified using the
`FEATURES_FILE` and `LABELS_FILE` parameters.

Each row of the `FEATURES_FILE` should contain the feature vector for
a drug-target pair. The `LABELS_FILE` will label the drug-target pair
in the corresponding row as interactive (1) or non-interactive (0).

#### Step 4: Initial Data Sharing

On the respective machines, `cd` into `mpc/code/` and run `MaskDnn`
for each party in the following order:

CP0: `bin/MaskDnn 0 ../par/test.par.0.txt`
CP1: `bin/MaskDnn 1 ../par/test.par.1.txt`
CP2: `bin/MaskDnn 2 ../par/test.par.2.txt`
SP:  `bin/MaskDnn 3 ../par/test.par.3.txt ../test_data/`

During this step, SP computes the secret shares with CP1 and CP2.
The resulting shares are stored in the same directory as the
features and the labels, and can be transmitted using a file
sharing protocol such as FTP.

#### Step 5: Secure DTI Prediction

On the respective machines, `cd` into `mpc/code/` and run `DnnClient`
for each party (excluding SP) in the following order:

CP0: `bin/DnnClient 0 ../par/test.par.0.txt`
CP1: `bin/DnnClient 1 ../par/test.par.1.txt`
CP2: `bin/DnnClient 2 ../par/test.par.2.txt`

As the model runs, it will output the current model parameters in
plaintext in the cache/ directory. This can be modified to only
output the model parameters at the end of model training.

**Contact for questions:**

Brian Hie, brianhie@mit.edu  
Hoon Cho, hhcho@mit.edu
