# Secure DTI Demo

This tutorial walks through the steps outlined in `demo.sh` and gives an overview of all stages of the Secure DTI pipeline.

## Data set

We created a small drug-target interaction data set based on three different drug targets from the [STITCH data set](http://stitch.embl.de/): ABL1, a tyrosine kinase target for chronic myeloid leukemia; HMGCR, the rate limiting enzyme in cholesterol synthesis; and OPRM1, the mu 1 opioid receptor.

For each of the above targets, we picked 33 chemicals (34 for ABL1) with the most confident interactions to get to 100 total positive interactions. This data is put in `demo_data`:

```
demo_data/
└── stitch
    ├── 9606.protein_chemical.links.v5.0.tsv.gz
    ├── ensp_to_bitvec.txt
    └── filtered_chem.txt
```

## Prepare data for neural network training

We want to balance our positive interactions with an equal number of drug-target pairs that do not interact (our negative set). To do this, we can randomly pair drugs and targets and assume they do not interact (which is true most of the time).

To generate a better negative set, we can also preserve the degree distributions of the drugs/targets, so that a drug or a target is equally represented in the negative set as it is in the positive set. This helps prevent a model from, e.g., simply learning to label as positive any drug-target pair with a drug that is frequently seen in the positive set but rarely seen in the negative set.

We can generate such a negative set with the command:
```
sh bin/generate_batches_pw.sh demo_
```

This will create files `demo_data/X.txt` and `demo_data/y.txt` with the positive interactions randomly shuffled and equally balanced with negative interactions that preserve the representation of drugs and targets. Note that the naive uniform sampling approach is also provided in ``bin/generate_batches_pw.sh``.

Next, we split our data into training and test data using the command:
```
sh bin/train_test_demo.sh
```

This will treat the first 150 drug-target pairs as the training set (75% of all pairs) and the last 50 pairs as the test set (25%). These files can be found in `demo_data/batch_pw/` and are called `Xtrain`, `ytrain`, `Xtest`, and `ytest`.


## Secret share the input data

We are ready to start with the multiparty computation (MPC) portion of the pipeline! Change directories into `mpc/code/` and generate random keys for secret sharing:
```
cd mpc/code/
./bin/GenerateKey ../key/P0_P1.key
./bin/GenerateKey ../key/P0_P2.key
./bin/GenerateKey ../key/P1_P2.key
./bin/GenerateKey ../key/P1_P3.key
./bin/GenerateKey ../key/P2_P3.key
./bin/GenerateKey ../key/global.key
```

These keys are used to secure the communication between pairs of computing parties and also to set up shared pseudorandom number generators that make our protocols more efficient. Technically, each of these keys should be generated and owned by only the relevant computing parties, but for demonstration purposes we are just keeping them on one machine for now.

Now, we share the underlying drug-target pairs from the SP to the CPs:
```
./bin/ShareData 0 ../par/demo.par.0.txt &
./bin/ShareData 1 ../par/demo.par.1.txt &
./bin/ShareData 2 ../par/demo.par.2.txt &
./bin/ShareData 3 ../par/demo.par.3.txt ../../demo_data/batch_pw/ &
wait
```

Where the id 3 refers to SP, 1 and 2 are CP1 and CP2 that handle the main computation, and 0 refers to CP0, an auxiliary party we introduce for precomputation. In practice, we may have multiple parties contributing data (perhaps the CPs themselves), but here we work with a single data contributor (i.e., SP) for simplicity.

The process for the SP generates "masked" files in `../../demo_data/batch_pw/` that contain the blinded drug-target pairs that can be sent to CP1 using some kind of file transfer protocol (e.g., `ftp`). On the other hand, CP2 reconstructs its share of the data using the seed files saved in `mpc/cache/` (in this case, `test_seedtrain.bin`). Again, these should technically be distributed across multiple machines, but in this demo they are all on the same machine, which allows us to skip this file transfer step.

## Train the model

We are now ready to train the neural network! That can be done with the commands:
```
./bin/TrainSecureDTI 0 ../par/demo.par.0.txt &
./bin/TrainSecureDTI 1 ../par/demo.par.1.txt &
./bin/TrainSecureDTI 2 ../par/demo.par.2.txt &
wait
```

This will train a simple feedforward neural network via stochsatic gradient descent to classify drug-target interactions from non-interactive drug-target pairs, while keeping the underlying data private. On our machine, this process completes in around 15 minutes.

## Evaluate the model

When model training is done, it will save the parameters in plaintext to `mpc/cache/`, which can then be used to evaluate the classification performance of the algorithm on both the training set and the test set. In practice, the trained model may be saved as secret shares to allow even the prediction of novel interactions to be performed in a privacy-preserving manner. 

Change back into the main directory and evaluate the model with the commands:
```
cd ../..
python bin/evaluate_demo.py
```

This will output something like:
```
Training accuracy:
2018-07-29 23:40:09.027171 | Batch accuracy: 1.0
2018-07-29 23:40:09.029643 | Accuracy: 1.00
2018-07-29 23:40:09.030162 | F1: 1.00
2018-07-29 23:40:09.032121 | Precision: 1.00
2018-07-29 23:40:09.032924 | Recall: 1.00
2018-07-29 23:40:09.033697 | ROC AUC: 1.00
2018-07-29 23:40:09.035296 | Avg. precision: 1.00
Testing accuracy:
2018-07-29 23:40:09.657741 | Batch accuracy: 0.96
2018-07-29 23:40:09.658088 | Accuracy: 0.96
2018-07-29 23:40:09.658368 | F1: 0.96
2018-07-29 23:40:09.659097 | Precision: 0.93
2018-07-29 23:40:09.659831 | Recall: 1.00
2018-07-29 23:40:09.660555 | ROC AUC: 0.97
2018-07-29 23:40:09.661240 | Avg. precision: 0.96
```

The report gives various metrics for gauging classification accuracy. The values might vary a little bit depending on the (pseudo)random assignment of drug-target pairs to the training or test sets, but in general the training and testing accuracy should be close to 100%.

In the above case, the trained model achieves near-perfect classification on the test data. The model was able to learn from the training data, while keeping everything private throughout the entire training procedure!

This concludes our demonstration of Secure DTI pipeline on a small data set. The main README gives instructions on how to run the pipeline on the full STITCH data set with many more drug-target pairs. Our pipeline can also be easily applied to other data sets.
