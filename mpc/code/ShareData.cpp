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

bool mask_matrix(string data_dir, MPCEnv& mpc, string name,
                 size_t n_rows, size_t n_cols) {
  /* Open file. */
  string fname = data_dir + name;
  ifstream fin(fname.c_str());
  if (!fin.is_open()) {
    tcout() << "Error: could not open " << fname << endl;
    return false;
  }

  /* Read in matrix. */
  Mat<ZZ_p> matrix;
  Init(matrix, n_rows, n_cols);
  
  string line;
  int i = 0;
  while(getline(fin, line)) {
    if (i % 1000 == 0) {
      tcout() << "Reading line " << i << endl;
    }
    for (int j = 0; j < n_cols; j++) {
      double val = ((double) line[j]) - 48;
      ZZ_p val_fp;
      DoubleToFP(val_fp, val, Param::NBIT_K, Param::NBIT_F);
      matrix[i][j] = val_fp;
    }
    i++;
  }
  
  if (i != n_rows) {
    tcout() << "Error: Invalid number of rows: " << i << endl;
    return false;
  }
  fin.close();

  /* Mask matrix. */
  Mat<ZZ_p> mask;
  mpc.RandMat(mask, n_rows, n_cols);
  matrix -= mask; /* Masked `matrix' should be sent to CP2. */

  /* Save data to file. */
  fstream fs;
  fs.open((data_dir + name + "_masked.bin").c_str(),
          ios::out | ios::binary);
  mpc.WriteToFile(matrix, fs);
  fs.close();
  tcout() << "Finished writing to file." << endl;

  return true;
}

bool mask_data(string data_dir, MPCEnv& mpc) {
  vector<string> suffixes;
  suffixes = load_suffixes(Param::TRAIN_SUFFIXES);
  
  mpc.SwitchSeed(1); /* Use CP1's seed. */

  fstream fs;
  string fname;
  for (int i = 0; i < suffixes.size(); i++) {
    /* Save seed state to file for each batch. */
    fname = cache(1, "seed" + suffixes[i]);
    fs.open(fname.c_str(), ios::out | ios::binary);
    if (!fs.is_open()) {
      tcout() << "Error: could not open " << fname << endl;
      return false;
    }
    mpc.ExportSeed(fs);
    fs.close();

    /* Write batch to file. */
    if (!mask_matrix(data_dir, mpc, "X" + suffixes[i],
                     Param::N_FILE_BATCH, Param::FEATURE_RANK))
      return false;
  
    if (!mask_matrix(data_dir, mpc, "y" + suffixes[i],
                     Param::N_FILE_BATCH, Param::N_CLASSES - 1))
      return false;
  }
  
  mpc.RestoreSeed();
  
  return true;
}

int main(int argc, char* argv[]) {
  if (argc < 3) {
    tcout() << "Usage: ShareData party_id param_file [data_dir (for P3/SP)]" << endl;
    return 1;
  }

  /* Load party id. */
  string pid_str(argv[1]);
  int pid;
  if (!Param::Convert(pid_str, pid, "party_id") || pid < 0 || pid > 3) {
    tcout() << "Error: party_id should be 0, 1, 2, or 3" << endl;
    return 1;
  }

  /* Load values in parameter file. */
  if (!Param::ParseFile(argv[2])) {
    tcout() << "Could not finish parsing parameter file" << endl;
    return 1;
  }

  /* Load data directory name. */
  string data_dir;
  if (pid == 3) {
    if (argc < 4) {
      tcout() << "Error: for P3/SP, data directory should be provided as the last argument" << endl;
      return 1;
    }
    data_dir = argv[3];
    if (data_dir[data_dir.size() - 1] != '/') {
      data_dir += "/";
    }
    tcout() << "Data directory: " << data_dir << endl;
  }

  /* Initalize MPC environment. */
  vector< pair<int, int> > pairs;
  pairs.push_back(make_pair(0, 1));
  pairs.push_back(make_pair(0, 2));
  pairs.push_back(make_pair(1, 2));
  pairs.push_back(make_pair(1, 3));
  pairs.push_back(make_pair(2, 3));
  MPCEnv mpc;
  if (!mpc.Initialize(pid, pairs)) {
    tcout() << "MPC environment initialization failed" << endl;
    return 1;
  }

  /* Mask the data and save to file. */
  bool success = true;  
  if (pid == 3) {
    success = mask_data(data_dir, mpc);
    if (!success) {
      tcout() << "Data masking failed." << endl;
    } else {
      tcout() << "Party 3 done streaming data." << endl;
    }
    mpc.SendBool(true, 2);

  } else if (pid == 2) {
    /* Keep CP2 alive until SP has shared data with it. */
    mpc.ReceiveBool(3);
    success = true;

  } else if (pid == 1) {
    /* CP1 needs to save its seed so it can reconstruct
       the masked data.*/
    string fname = cache(pid, "initial_seed");
    fstream fs;
    fs.open(fname.c_str(), ios::out | ios::binary);
    if (!fs.is_open()) {
      tcout() << "Error: could not open " << fname << endl;
      return false;
    }
    mpc.SwitchSeed(3);
    mpc.ExportSeed(fs);
    mpc.RestoreSeed();    
    fs.close();
    success = true;
  }

  /* Keep party 0 online until end of data masking. */
  if (pid == 0) {
    mpc.ReceiveBool(2);
  } else if (pid == 2) {
    mpc.SendBool(true, 0);
  }

  mpc.CleanUp();

  if (success) {
    tcout() << "Protocol successfully completed." << endl;
    return 0;
  } else {
    tcout() << "Protocol abnormally terminated." << endl;
    return 1;
  }
}
