#include <algorithm>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <map>
#include <random>
#include <sstream>

#include "connect.h"
#include "mpc.h"
#include "protocol.h"
#include "util.h"
#include "NTL/ZZ_p.h"

using namespace NTL;
using namespace std;

void reveal(ZZ_p X, string fname, MPCEnv& mpc) {
  mpc.RevealSym(X);
  double X_double;
  fstream fs;
  fs.open(fname.c_str(), ios::out);
  X_double = FPToDouble(X, Param::NBIT_K, Param::NBIT_F);
  fs << X_double;
  fs.close();
}

void reveal(Vec<ZZ_p> X, string fname, MPCEnv& mpc) {
  mpc.RevealSym(X);
  Vec<double> X_double;
  fstream fs;
  fs.open(fname.c_str(), ios::out);
  FPToDouble(X_double, X, Param::NBIT_K, Param::NBIT_F);
  for (int i = 0; i < X.length(); i++) {
    fs << X_double[i] << '\t';
  }
  fs.close();
}

void reveal(Mat<ZZ_p> X, string fname, MPCEnv& mpc) {
  mpc.RevealSym(X);
  Mat<double> X_double;
  fstream fs;
  fs.open(fname.c_str(), ios::out);
  FPToDouble(X_double, X, Param::NBIT_K, Param::NBIT_F);
  for (int i = 0; i < X.NumRows(); i++) {
    for (int j = 0; j < X.NumCols(); j++) {
      fs << X_double[i][j] << '\t';
    }
    fs << endl;
  }
  fs.close();
}

bool read_matrix(Mat<ZZ_p>& matrix, ifstream& ifs, string fname,
                 size_t n_rows, size_t n_cols, MPCEnv& mpc) {
  ifs.open(fname.c_str(), ios::in | ios::binary);
  if (!ifs.is_open()) {
    tcout() << "Could not open : " << fname << endl;
    return false;
  }
  mpc.ReadFromFile(matrix, ifs, n_rows, n_cols);
  ifs.close();
  return true;
}  

void initialize_model(vector<Mat<ZZ_p> >& W, vector<Vec<ZZ_p> >& b,
                      vector<Mat<ZZ_p> >& dW, vector<Vec<ZZ_p> >& db,
                      vector<Mat<ZZ_p> >& vW, vector<Vec<ZZ_p> >& vb,
                      int pid, MPCEnv& mpc) {
  /* Random number generator for Gaussian noise
     initialization of weight matrices. */
  std::default_random_engine generator (0);
  std::normal_distribution<double> distribution (0.0, 0.01);

  for (int l = 0; l < Param::N_HIDDEN + 1; l++) {
    Mat<ZZ_p> W_layer, dW_layer, vW_layer;
    Vec<ZZ_p> b_layer, db_layer, vb_layer;

    /* Handle case with 0 hidden layers. */
    if (Param::N_HIDDEN == 0 && l >= 1) {
      break;
    } else if (Param::N_HIDDEN == 0 && l == 0) {
      W_layer.SetDims(Param::FEATURE_RANK, Param::N_CLASSES - 1);
      b_layer.SetLength(Param::N_CLASSES - 1);
    
    /* Set dimensions of the input layer. */
    } else if (l == 0) {
      W_layer.SetDims(Param::FEATURE_RANK, Param::N_NEURONS);
      b_layer.SetLength(Param::N_NEURONS);

    /* Set dimensions of the output layer. */
    } else if (l == Param::N_HIDDEN) {
      W_layer.SetDims(Param::N_NEURONS, Param::N_CLASSES - 1);
      b_layer.SetLength(Param::N_CLASSES - 1);
      
    /* Set dimensions of the hidden layers. */
    } else {
      W_layer.SetDims(Param::N_NEURONS, Param::N_NEURONS);
      b_layer.SetLength(Param::N_NEURONS);
    }
    
    dW_layer.SetDims(W_layer.NumRows(), W_layer.NumCols());
    Init(vW_layer, W_layer.NumRows(), W_layer.NumCols());
    
    db_layer.SetLength(b_layer.length());
    Init(vb_layer, b_layer.length());
     
    Mat<ZZ_p> W_r;
    Vec<ZZ_p> b_r;
    if (pid == 2) {
      /* CP2 will have real data minus random data. */
      /* Initialize weight matrix with Gaussian noise. */
      for (int i = 0; i < W_layer.NumRows(); i++) {
        for (int j = 0; j < W_layer.NumCols(); j++) {
          double noise = distribution(generator);
          DoubleToFP(W_layer[i][j], noise, Param::NBIT_K, Param::NBIT_F);
        }
      }
      
      Init(b_layer, b_layer.length());

      /* Blind the data. */
      mpc.SwitchSeed(1);
      mpc.RandMat(W_r, W_layer.NumRows(), W_layer.NumCols());
      mpc.RandVec(b_r, b_layer.length());
      mpc.RestoreSeed();
      W_layer -= W_r;
      b_layer -= b_r;
      
    } else if (pid == 1) {
      /* CP1 will just have the random data. */
      mpc.SwitchSeed(2);
      mpc.RandMat(W_r, W_layer.NumRows(), W_layer.NumCols());
      mpc.RandVec(b_r, b_layer.length());
      mpc.RestoreSeed();
      W_layer = W_r;
      b_layer = b_r;
    }
    
    W.push_back(W_layer);
    dW.push_back(dW_layer);
    vW.push_back(vW_layer);
    b.push_back(b_layer);
    db.push_back(db_layer);
    vb.push_back(vb_layer);
  }
}

void gradient_descent(Mat<ZZ_p>& X, Mat<ZZ_p>& y,
                      vector<Mat<ZZ_p> >& W, vector<Vec<ZZ_p> >& b,
                      vector<Mat<ZZ_p> >& dW, vector<Vec<ZZ_p> >& db,
                      vector<Mat<ZZ_p> >& vW, vector<Vec<ZZ_p> >& vb,
                      vector<Mat<ZZ_p> >& act, vector<Mat<ZZ_p> >& relus,
                      int epoch, int pid, MPCEnv& mpc) {
  if (pid == 2)
    tcout() << "Epoch: " << epoch << endl;
      
  /************************
   * Forward propagation. *
   ************************/
  for (int l = 0; l < Param::N_HIDDEN; l++) {

    if (pid == 2)
      tcout() << "Forward prop, multiplication." << endl;
    /* Multiply weight matrix. */
    Mat<ZZ_p> activation;
    if (l == 0) {
      mpc.MultMat(activation, X, W[l]);
    } else {
      mpc.MultMat(activation, act[l-1], W[l]);
    }
    mpc.Trunc(activation);

    /* Add bias term; */
    for (int i = 0; i < activation.NumRows(); i++) {
      activation[i] += b[l];
    }

    if (pid == 2)
      tcout() << "Forward prop, ReLU." << endl;
    /* Apply ReLU non-linearity. */
    Mat<ZZ_p> relu;
    mpc.IsPositive(relu, activation);
    Mat<ZZ_p> after_relu;
    assert(activation.NumRows() == relu.NumRows());
    assert(activation.NumCols() == relu.NumCols());
    mpc.MultElem(after_relu, activation, relu);
    /* Note: Do not call Trunc() here because IsPositive()
       returns a secret shared integer, not a fixed point.*/
    
    // TODO: Implement dropout.

    /* Save activation for backpropagation. */
    act.push_back(after_relu);
    relus.push_back(relu);
  }

  
  /**************************
   * Evaluate class scores. *
   **************************/
  if (pid == 2)
    tcout() << "Score computation." << endl;
  Mat<ZZ_p> scores;
  if (Param::N_HIDDEN == 0) {
    mpc.MultMat(scores, X, W.back());
  } else {
    mpc.MultMat(scores, act.back(), W.back());
  }
  mpc.Trunc(scores);

  /* Add bias term; */
  for (int i = 0; i < scores.NumRows(); i++) {
    scores[i] += b.back();
  }

  Mat<ZZ_p> dscores;
  if (Param::LOSS == "hinge") {
    /* Scale y to be -1 or 1. */
    y *= 2;
    if (pid == 2) {
      for (int i = 0; i < y.NumRows(); i++) {
        for (int j = 0; j < y.NumCols(); j++) {
          y[i][j] -= DoubleToFP(1, Param::NBIT_K, Param::NBIT_F);
        }
      }
    }

    /* Compute 1 - y * scores. */
    y *= -1;
    Mat<ZZ_p> mod_scores;
    mpc.MultElem(mod_scores, y, scores);
    mpc.Trunc(mod_scores);
    if (pid == 2) {
      for (int i = 0; i < mod_scores.NumRows(); i++) {
        for (int j = 0; j < mod_scores.NumCols(); j++) {
          mod_scores[i][j] += DoubleToFP(1, Param::NBIT_K, Param::NBIT_F);
        }
      }
    }

    /* Compute hinge loss and derivative. */
    Mat<ZZ_p> hinge;
    mpc.IsPositive(hinge, mod_scores);
    mpc.MultElem(dscores, y, hinge);
    /* Note: No need to not call Trunc(). */
    
  } else {
    /* Compute derivative of the scores using MSE loss. */
    dscores = scores - y;
  }
  
  ZZ_p norm_examples;
  DoubleToFP(norm_examples, 1. / ((double) X.NumRows()),
             Param::NBIT_K, Param::NBIT_F);
  dscores *= norm_examples;
  mpc.Trunc(dscores);

  if (pid > 0 && epoch % 500 == 0) {
    reveal(scores, cache(pid, "scores_epoch" + to_string(epoch)), mpc);
  }

  /*********************
   * Back propagation. *
   *********************/
  Mat<ZZ_p> dhidden = dscores;
  for (int l = Param::N_HIDDEN; l >= 0; l--) {
    
    if (pid == 2)
      tcout() << "Back prop, multiplication." << endl;
    /* Compute derivative of weights. */
    Init(dW[l], W[l].NumRows(), W[l].NumCols());
    Mat<ZZ_p> X_T;
    if (l == 0) {
      X_T = transpose(X);
    } else {
      X_T = transpose(act.back());
      act.pop_back();
    }
    mpc.MultMat(dW[l], X_T, dhidden);
    mpc.Trunc(dW[l]);
  
    /* Add regularization term to weights. */
    ZZ_p REG;
    DoubleToFP(REG, Param::REG, Param::NBIT_K, Param::NBIT_F);
    Mat<ZZ_p> reg = W[l] * REG;
    mpc.Trunc(reg);
    dW[l] += reg;

    /* Compute derivative of biases. */
    Init(db[l], b[l].length());
    for (int i = 0; i < dhidden.NumRows(); i++) {
      db[l] += dhidden[i];
    }

    if (l > 0) {
      /* Compute backpropagated activations. */
      Mat<ZZ_p> dhidden_new, W_T;
      W_T = transpose(W[l]);
      mpc.MultMat(dhidden_new, dhidden, W_T);
      mpc.Trunc(dhidden_new);

      if (pid == 2)
        tcout() << "Back prop, ReLU." << endl;
      /* Apply derivative of ReLU. */
      Init(dhidden, dhidden_new.NumRows(), dhidden_new.NumCols());
      mpc.MultElem(dhidden, dhidden_new, relus.back());
      /* Note: No need to not call Trunc().*/
      relus.pop_back();
    }
  }

  assert(act.size() == 0);
  assert(relus.size() == 0);

  if (pid == 2)
    tcout() << "Momentum update." << endl;
  /* Update the model using Nesterov momentum. */
  /* Compute constants that update various parameters. */
  ZZ_p MOMENTUM = DoubleToFP(Param::MOMENTUM,
                             Param::NBIT_K, Param::NBIT_F);
  ZZ_p MOMENTUM_PLUS1 = DoubleToFP(Param::MOMENTUM + 1,
                                   Param::NBIT_K, Param::NBIT_F);
  ZZ_p LEARN_RATE = DoubleToFP(Param::LEARN_RATE,
                               Param::NBIT_K, Param::NBIT_F);

  for (int l = 0; l < Param::N_HIDDEN + 1; l++) {
    /* Update the weights. */
    Mat<ZZ_p> vW_prev = vW[l];
    vW[l] = (MOMENTUM * vW[l]) - (LEARN_RATE * dW[l]);
    mpc.Trunc(vW[l]);
    Mat<ZZ_p> W_update = (-MOMENTUM * vW_prev) + (MOMENTUM_PLUS1 * vW[l]);
    mpc.Trunc(W_update);
    W[l] += W_update;

    /* Update the biases. */
    Vec<ZZ_p> vb_prev = vb[l];
    vb[l] = (MOMENTUM * vb[l]) - (LEARN_RATE * db[l]);
    mpc.Trunc(vb[l]);
    Vec<ZZ_p> b_update = (-MOMENTUM * vb_prev) + (MOMENTUM_PLUS1 * vb[l]);
    mpc.Trunc(b_update);
    b[l] += b_update;
  }
}

void load_X_y(string suffix, Mat<ZZ_p>& X, Mat<ZZ_p>& y,
              int pid, MPCEnv& mpc) {
  if (pid == 0)
    /* Matrices must also be initialized even in CP0,
       but they do not need to be filled. */
    return;
  ifstream ifs;
  
  /* Load seed for CP1. */
  if (pid == 1) {
    string fname = "../cache/test_seed" + suffix + ".bin";
    ifs.open(fname.c_str(), ios::binary);
    if (!ifs.is_open()) {
      tcout() << "Error: could not open " << fname << endl;
      return;
    }
    mpc.ImportSeed(20, ifs);
    ifs.close();
  }

  if (pid == 2) {
    /* In CP2, read in blinded matrix. */
    tcout() << "reading in " << Param::FEATURES_FILE << suffix << endl;
    if (!read_matrix(X, ifs, Param::FEATURES_FILE + suffix + "_masked.bin",
                     X.NumRows(), X.NumCols(), mpc))
      return;
    
    tcout() << "reading in " << Param::LABELS_FILE << suffix << endl;
    if (!read_matrix(y, ifs, Param::LABELS_FILE + suffix + "_masked.bin",
                     y.NumRows(), y.NumCols(), mpc))
      return;

  } else if (pid == 1) {
    /* In CP1, use seed to regenerate blinding factors.
       These need to be generated in the same order as the
       original blinding factors! */
    mpc.SwitchSeed(20);
    mpc.RandMat(X, X.NumRows(), X.NumCols());
    mpc.RandMat(y, y.NumRows(), y.NumCols());
    mpc.RestoreSeed();
  }
}

void model_update(Mat<ZZ_p>& X, Mat<ZZ_p>& y,
                  vector<Mat<ZZ_p> >& W, vector<Vec<ZZ_p> >& b,
                  vector<Mat<ZZ_p> >& dW, vector<Vec<ZZ_p> >& db,
                  vector<Mat<ZZ_p> >& vW, vector<Vec<ZZ_p> >& vb,
                  vector<Mat<ZZ_p> >& act, vector<Mat<ZZ_p> >& relus,
                  int& epoch, int pid, MPCEnv& mpc) {

  /* Round down number of batches in file. */
  int batches_in_file = X.NumRows() / Param::BATCH_SIZE;
  Mat<ZZ_p> X_batch;
  Mat<ZZ_p> y_batch;
  X_batch.SetDims(Param::BATCH_SIZE, X.NumCols());
  y_batch.SetDims(Param::BATCH_SIZE, y.NumCols());

  vector<int> random_idx(X.NumRows());
  iota(random_idx.begin(), random_idx.end(), 0);
  random_shuffle(random_idx.begin(), random_idx.end());

  for (int i = 0; i < batches_in_file; i++) {

    /* Scan matrix (pre-shuffled) to get batch. */
    int base_j = i * Param::BATCH_SIZE;
    for (int j = base_j;
         j < base_j + Param::BATCH_SIZE && j < X.NumRows();
         j++) {
      X_batch[j - base_j] = X[random_idx[j]];
      y_batch[j - base_j] = y[random_idx[j]];
    }
    
    /* Do one round of mini-batch gradient descent. */
    gradient_descent(X_batch, y_batch,
                     W, b, dW, db, vW, vb, act, relus,
                     epoch, pid, mpc);

    /* Save state every 500 epochs. */
    if (epoch % 500 == 0) {
      for (int l = 0; l < Param::N_HIDDEN + 1; l++) {
        Mat<ZZ_p> W_out;
        Init(W_out, W[l].NumRows(), W[l].NumCols());
        W_out += W[l];
        reveal(W_out, cache(pid, "W" + to_string(l) + "_" +
                            to_string(epoch)), mpc);
        
        Vec<ZZ_p> b_out;
        Init(b_out, b[l].length());
        b_out += b[l];
        reveal(b_out, cache(pid, "b" + to_string(l) + "_" +
                            to_string(epoch)), mpc);
      }
    }

    /* Update reference to training epoch. */
    epoch++;

    if (epoch >= Param::MAX_EPOCHS) {
      break;
    }
  }
}

bool dti_protocol(MPCEnv& mpc, int pid) {
  /* Initialize threads. */
  SetNumThreads(Param::NUM_THREADS);
  tcout() << AvailableThreads() << " threads created" << endl;

  /* Initialize model and data structures. */
  tcout() << "Initializing model." << endl;
  vector<Mat<ZZ_p> > W, dW, vW, act, relus;
  vector<Vec<ZZ_p> > b, db, vb;
  initialize_model(W, b, dW, db, vW, vb, pid, mpc);

  srand(0);  /* Seed 0 to have deterministic testing. */

  /* Create list of training file suffixes. */
  vector<string> suffixes;
  suffixes = load_suffixes(Param::TRAIN_SUFFIXES);

  /* Initialize data matries. */
  Mat<ZZ_p> X, y;
  X.SetDims(Param::N_FILE_BATCH, Param::FEATURE_RANK);
  y.SetDims(Param::N_FILE_BATCH, Param::N_CLASSES - 1);
  string suffix = suffixes[rand() % suffixes.size()];  
  load_X_y(suffix, X, y, pid, mpc);
  
  /* Do gradient descent over multiple training epochs. */
  for (int epoch = 0; epoch < Param::MAX_EPOCHS;
       /* model_update() updates epoch. */) {
    /* Do model updates and file reads in parallel. */
    model_update(X, y, W, b, dW, db, vW, vb, act, relus,
                 epoch, pid, mpc);

    suffix = suffixes[rand() % suffixes.size()];

    load_X_y(suffix, X, y, pid, mpc);
  }
  
  if (pid > 0) {
    for (int l = 0; l < Param::N_HIDDEN + 1; l++) {
      Mat<ZZ_p> W_out;
      Init(W_out, W[l].NumRows(), W[l].NumCols());
      W_out += W[l];
      reveal(W_out, cache(pid, "W" + to_string(l) + "_final"), mpc);
        
      Vec<ZZ_p> b_out;
      Init(b_out, b[l].length());
      b_out += b[l];
      reveal(b_out, cache(pid, "b" + to_string(l) + "_final"), mpc);
    }
  }
  
  return true;
}

int main(int argc, char** argv) {
  if (argc < 3) {
    tcout() << "Usage: TrainSecureDTI party_id param_file" << endl;
    return 1;
  }

  string pid_str(argv[1]);
  int pid;
  if (!Param::Convert(pid_str, pid, "party_id") || pid < 0 || pid > 2) {
    tcout() << "Error: party_id should be 0, 1, or 2" << endl;
    return 1;
  }

  if (!Param::ParseFile(argv[2])) {
    tcout() << "Could not finish parsing parameter file" << endl;
    return 1;
  }

  vector< pair<int, int> > pairs;
  pairs.push_back(make_pair(0, 1));
  pairs.push_back(make_pair(0, 2));
  pairs.push_back(make_pair(1, 2));

  /* Initialize MPC environment */
  MPCEnv mpc;
  if (!mpc.Initialize(pid, pairs)) {
    tcout() << "MPC environment initialization failed" << endl;
    return 1;
  }

  bool success = dti_protocol(mpc, pid);

  // This is here just to keep P0 online until the end for data transfer
  // In practice, P0 would send data in advance before each phase and go offline
  if (pid == 0) {
    mpc.ReceiveBool(2);
  } else if (pid == 2) {
    mpc.SendBool(true, 0);
  }

  mpc.CleanUp();

  if (success) {
    tcout() << "Protocol successfully completed" << endl;
    return 0;
  } else {
    tcout() << "Protocol abnormally terminated" << endl;
    return 1;
  }
}
