#include "connect.h"
#include "mpc.h"
#include "protocol.h"
#include "util.h"
#include "NTL/ZZ_p.h"

#include <cstdlib>
#include <fstream>
#include <map>
#include <iostream>
#include <sstream>

using namespace NTL;
using namespace std;

/* This is a protocol for demonstrating how the collaborating entities
 * can perform additional computation to better interpret the prediction
 * results from the Secure DTI model while maintaining the same
 * privacy guarantees as the main protocol.
 *
 * After the initial training and testing phases of the Secure DTI
 * pipeline, the entities select top few predictions they are
 * interested in and securely evaluate the prediction score for 
 * modified feature vectors where individual features are masked
 * one-at-a-time to observe their impact on the prediction score.
 *
 * We demonstrate this using a pre-trained model and our top predictions
 * from the STITCH database as reported in our manuscript (Table 1a).
 * We mask each of the annotated domains of the target protein and
 * reveal the perturbed prediction score.
 * 
 */

bool feature_analysis_protocol(MPCEnv& mpc, int pid) {

  SetNumThreads(Param::NUM_THREADS);
  cout << AvailableThreads() << " threads created" << endl;

  bool debug = false;

  int cdim = 1024; // chemical feature dimension
  int pdim = 5879; // protein feature dimension
  int nfeat = cdim + pdim; // total # of features

  // Network size
  int nunit = 250; 
  int nlayer = 2;

  // Model parameter files
  string model_prefix = "../../models/stitch_a/";

  // Drug-target pairs of interest
  string predfile = "../../models/stitch_a/top12_pred_features.txt";

  ofstream ofs;
  ifstream ifs;
  Vec<ZZ_p> tmp_vec;

  // Model parameter variables
  Vec< Mat<ZZ_p> > W;
  W.SetLength(nlayer + 1);
  Vec< Vec<ZZ_p> > b;
  b.SetLength(nlayer + 1);

  // Read in the saved parameters by layer
  for (int i = 0; i <= nlayer; i++) {
    int indim = (i == 0) ? nfeat : nunit;
    int outdim = (i == nlayer) ? 1 : nunit;

    W[i].SetDims(outdim, indim);
    b[i].SetLength(outdim);

    // Read from saved file. We provided the parameters in plaintext,
    // but here we revert them back to secret shares
    // to simulate the scenario where the model is kept private
    // throughout this analysis
    if (pid == 1) {
      ostringstream oss;
      oss << model_prefix << "W" << i << ".txt";
      string Wfile = oss.str();

      oss.str("");
      oss << model_prefix << "b" << i << ".txt";
      string bfile = oss.str();

      ifs.open(Wfile.c_str());
      double val;
      for (int rowi = 0; rowi < indim; rowi++) {
        for (int coli = 0; coli < outdim; coli++) {
          ifs >> val;
          W[i][coli][rowi] = DoubleToFP(val, Param::NBIT_K, Param::NBIT_F);
        }
      }
      ifs.close();

      ifs.open(bfile.c_str());
      for (int j = 0; j < outdim; j++) {
        ifs >> val;
        b[i][j] = DoubleToFP(val, Param::NBIT_K, Param::NBIT_F);
      }
      ifs.close();

      Mat<ZZ_p> Wr;
      Vec<ZZ_p> br;

      mpc.SwitchSeed(2);
      mpc.RandMat(Wr, outdim, indim);
      mpc.RandVec(br, outdim);
      mpc.RestoreSeed();

      W[i] -= Wr;
      b[i] -= br;

    } else if (pid == 2) {

      mpc.SwitchSeed(1);
      mpc.RandMat(W[i], outdim, indim);
      mpc.RandVec(b[i], outdim);
      mpc.RestoreSeed();

    }

    if (debug) {
      cout << "W" << i << ":" << endl;
      mpc.PrintFP(W[i], 3, 3);
      cout << "b" << i << ":" << endl;
      mpc.PrintFP(b[i], 3);
    }
    
  }

  // Beaver partition the weight matrices for fast reuse
  Vec< Mat<ZZ_p> > Wr, Wm;
  Wr.SetLength(nlayer + 1);
  Wm.SetLength(nlayer + 1);
  for (int i = 0; i <= nlayer; i++) {
    mpc.BeaverPartition(Wr[i], Wm[i], W[i]);
  }

  // Go through the data instances of interest and analyze them.
  // For simplicity, we use the input features provided in a plaintext file
  // to simulate the secret shares for CP1/2. In practice, each data 
  // owner will choose which drug-target pairs to investigate
  // further and secret share the necessary data for evaluation,
  // where individual features are set to zero. Which feature is set to
  // zero is not revealed to CP1/2, neither how many perturbations are
  // performed for a particular data instance. These aspects are reported in
  // the protocol below for demonstration purposes only.
  //  
  ifs.open(predfile.c_str());
  int i = 0;
  string line;
  while (ifs >> line) { // read until end of file
    i++;

    if (pid == 2) {
      cout << "[" << i << "]" << endl;
    }

    Vec<bool> xraw;
    xraw.SetLength(nfeat);
    for (int j = 0; j < nfeat; j++) {
      xraw[j] = line[j] == '1';
    }

    Vec<ZZ_p> x;
    x.SetLength(nfeat);

    // Go through each domain annotation and evaluate
    // the impact of excluding it
    // p = -1 means use all domains
    double base_score = 0;
    for (int p = -1; p < pdim; p++) {
      int o = cdim + p;
      if (!xraw[o] && p >= 0) continue;

      // Simulate secret shares of input
      if (pid == 1) {
        for (int j = 0; j < nfeat; j++) {
          x[j] = ZZ_p(xraw[j] ? 1 : 0);
          if (j == o && p >= 0) x[j] = 0; // Exclude domain
        }

        Vec<ZZ_p> xr;
        mpc.SwitchSeed(2);
        mpc.RandVec(xr, nfeat);
        mpc.RestoreSeed();

        x -= xr;
      } else if (pid == 2) {
        mpc.SwitchSeed(1);
        mpc.RandVec(x, nfeat);
        mpc.RestoreSeed();
      }

      // Evaluate model
      Vec<ZZ_p> z = x; 
      Vec<ZZ_p> zr, zm;
      for (int l = 0; l <= nlayer; l++) {
        int outdim = (l == nlayer) ? 1 : nunit;

        mpc.BeaverPartition(zr, zm, z);
        Init(z, outdim);
        mpc.BeaverMult(z, Wr[l], Wm[l], zr, zm);
        mpc.BeaverReconstruct(z);

        if (l > 0) {
          mpc.Trunc(z);
        }

        z += b[l];

        if (l < nlayer) { // ReLU activation
          Vec<ZZ_p> c;
          mpc.IsPositive(c, z);
          mpc.MultElem(z, c, z);
        }
      }

      // Output prediction score
      mpc.RevealSym(z);
      double score = FPToDouble(z[0], Param::NBIT_K, Param::NBIT_F);
      if (pid == 2) {
        if (p < 0) {
          base_score = score;
          cout << "Full score: " << score << endl;
        } else {
          cout << "- Domain " << p+1 << ": " << score << " (" << 100*(score-base_score)/base_score << "%)" << endl;
        }
      }
    }
  }
    
  ifs.close();

  return true;
}

int main(int argc, char** argv) {
  if (argc < 3) {
    cout << "Usage: AnalyzeFeature party_id param_file" << endl;
    return 1;
  }

  string pid_str(argv[1]);
  int pid;
  if (!Param::Convert(pid_str, pid, "party_id") || pid < 0 || pid > 2) {
    cout << "Error: party_id should be 0, 1, or 2" << endl;
    return 1;
  }

  if (!Param::ParseFile(argv[2])) {
    cout << "Could not finish parsing parameter file" << endl;
    return 1;
  }

  vector< pair<int, int> > pairs;
  pairs.push_back(make_pair(0, 1));
  pairs.push_back(make_pair(0, 2));
  pairs.push_back(make_pair(1, 2));

  /* Initialize MPC environment */
  MPCEnv mpc;
  if (!mpc.Initialize(pid, pairs)) {
    cout << "MPC environment initialization failed" << endl;
    return 1;
  }

  bool success = feature_analysis_protocol(mpc, pid);

  // This is here just to keep P0 online until the end for data transfer
  // In practice, P0 would send data in advance before each phase and go offline
  if (pid == 0) {
    mpc.ReceiveBool(2);
  } else if (pid == 2) {
    mpc.SendBool(true, 0);
  }

  mpc.CleanUp();

  if (success) {
    cout << "Protocol successfully completed" << endl;
    return 0;
  } else {
    cout << "Protocol abnormally terminated" << endl;
    return 1;
  }
}
