#include "param.h"
#include "util.h"
#include <fstream>
#include <sstream>

template<class T>
bool Param::Convert(string s, T &var, string name) {
  istringstream iss(s);
  iss >> std::skipws >> var >> std::skipws;
  if (iss.tellg() != -1) {
    tcout() << "Parameter parse error: " << name << endl;
    return false;
  }
  return true;
}

bool Param::ParseFile(const char *param_file) {
  ifstream pfile(param_file);
  if (!pfile.is_open()) {
    tcout() << "Failed to open parameter file: " << param_file << endl;
    return false;
  } else {
    tcout() << "Using parameters in " << param_file << endl;
  }

  string k, v;
  while (pfile >> std::skipws >> k) {
    if (k[0] == '#') { // Comments
      getline(pfile, v);
      continue;
    } 

    getline(pfile, v);

    bool ret;
    if (k == "PORT_P0_P1") {
      ret = Convert(v, Param::PORT_P0_P1, k);
    } else if (k == "PORT_P0_P2") {
      ret = Convert(v, Param::PORT_P0_P2, k);
    } else if (k == "PORT_P1_P2") {
      ret = Convert(v, Param::PORT_P1_P2, k);
    } else if (k == "PORT_P1_P3") {
      ret = Convert(v, Param::PORT_P1_P3, k);
    } else if (k == "PORT_P2_P3") {
      ret = Convert(v, Param::PORT_P2_P3, k);
    } else if (k == "IP_ADDR_P0") {
      ret = Convert(v, Param::IP_ADDR_P0, k);
    } else if (k == "IP_ADDR_P1") {
      ret = Convert(v, Param::IP_ADDR_P1, k);
    } else if (k == "IP_ADDR_P2") {
      ret = Convert(v, Param::IP_ADDR_P2, k);
    } else if (k == "KEY_PATH") {
      ret = Convert(v, Param::KEY_PATH, k);
    } else if (k == "NBIT_K") {
      ret = Convert(v, Param::NBIT_K, k);
    } else if (k == "NBIT_F") {
      ret = Convert(v, Param::NBIT_F, k);
    } else if (k == "NBIT_V") {
      ret = Convert(v, Param::NBIT_V, k);
    } else if (k == "BASE_P") {
      ret = Convert(v, Param::BASE_P, k);
    } else if (k == "DIV_MAX_N") {
      ret = Convert(v, Param::DIV_MAX_N, k);
    } else if (k == "PAR_THRES") {
      ret = Convert(v, Param::PAR_THRES, k);
    } else if (k == "NUM_THREADS") {
      ret = Convert(v, Param::NUM_THREADS, k);
    } else if (k == "MPC_BUF_SIZE") {
      ret = Convert(v, Param::MPC_BUF_SIZE, k);
    } else if (k == "ITER_PER_EVAL") {
      ret = Convert(v, Param::ITER_PER_EVAL, k);
    } else if (k == "OUTPUT_FILE_PREFIX") {
      ret = Convert(v, Param::OUTPUT_FILE_PREFIX, k);
    } else if (k == "LOG_FILE") {
      ret = Convert(v, Param::LOG_FILE, k);
    } else if (k == "CACHE_FILE_PREFIX") {
      ret = Convert(v, Param::CACHE_FILE_PREFIX, k);
      if (v[v.size() - 1] != '/') v += "/";
    } else if (k == "PROFILER") {
      ret = Convert(v, Param::PROFILER, k);
    } else if (k == "FEATURE_RANK") {
      ret = Convert(v, Param::FEATURE_RANK, k);
    } else if (k == "MAX_EPOCHS") {
      ret = Convert(v, Param::MAX_EPOCHS, k);
    } else if (k == "MOMENTUM") {
      ret = Convert(v, Param::MOMENTUM, k);
    } else if (k == "LEARN_RATE") {
      ret = Convert(v, Param::LEARN_RATE, k);
    } else if (k == "N_FILE_BATCH") {
      ret = Convert(v, Param::N_FILE_BATCH, k);
    } else if (k == "N_CLASSES") {
      ret = Convert(v, Param::N_CLASSES, k);
    } else if (k == "FEATURES_FILE") {
      ret = Convert(v, Param::FEATURES_FILE, k);
    } else if (k == "LABELS_FILE") {
      ret = Convert(v, Param::LABELS_FILE, k);
    } else if (k == "TRAIN_SUFFIXES") {
      ret = Convert(v, Param::TRAIN_SUFFIXES, k);
    } else if (k == "TEST_SUFFIXES") {
      ret = Convert(v, Param::TEST_SUFFIXES, k);
    } else if (k == "N_CLASSES") {
      ret = Convert(v, Param::N_CLASSES, k);
    } else if (k == "N_HIDDEN") {
      ret = Convert(v, Param::N_HIDDEN, k);
    } else if (k == "N_NEURONS") {
      ret = Convert(v, Param::N_NEURONS, k);
    } else if (k == "LEARN_RATE") {
      ret = Convert(v, Param::LEARN_RATE, k);
    } else if (k == "ANNEAL") {
      ret = Convert(v, Param::ANNEAL, k);
    } else if (k == "ANNEAL_FREQ") {
      ret = Convert(v, Param::ANNEAL_FREQ, k);
    } else if (k == "REG") {
      ret = Convert(v, Param::REG, k);
    } else if (k == "DROPOUT") {
      ret = Convert(v, Param::DROPOUT, k);
    } else if (k == "LOSS") {
      ret = Convert(v, Param::LOSS, k);
    } else if (k == "BATCH_SIZE") {
      ret = Convert(v, Param::BATCH_SIZE, k);
    } else if (k == "DEBUG") {
      ret = Convert(v, Param::DEBUG, k);
    } else {
      tcout() << "Unknown parameter: " << k << endl;
      ret = false;
    }

    if (!ret) {
      return false;
    }
  }

  return true;
}

int Param::PORT_P0_P1 = 8000;
int Param::PORT_P0_P2 = 8001;
int Param::PORT_P1_P2 = 8000;
int Param::PORT_P1_P3 = 8001;
int Param::PORT_P2_P3 = 8000;

string Param::IP_ADDR_P0 = "128.00.00.101";
string Param::IP_ADDR_P1 = "128.00.00.102";
string Param::IP_ADDR_P2 = "128.00.00.103";

string Param::KEY_PATH = "../key/";

int Param::NBIT_K = 60;
int Param::NBIT_F = 45;
int Param::NBIT_V = 64;

string Param::BASE_P = "1461501637330902918203684832716283019655932542929";

uint64_t Param::MPC_BUF_SIZE = 1000000;

int Param::ITER_PER_EVAL= 5;
string Param::OUTPUT_FILE_PREFIX = "../out/test";
string Param::CACHE_FILE_PREFIX = "../cache/test";
string Param::LOG_FILE = "../log/log.txt";

long Param::FEATURE_RANK = 100;
double Param::MOMENTUM = 0.9;
double Param::LEARN_RATE = 0.01;
int Param::N_CLASSES = 2;
int Param::N_HIDDEN = 2;
bool Param::ANNEAL = false;
int Param::ANNEAL_FREQ = 400;
double Param::REG = 0.001;
double Param::DROPOUT = 0;
string Param::LOSS = "hinge";
int Param::N_FILE_BATCH = 20000;
int Param::N_NEURONS = 250;
int Param::BATCH_SIZE = 50;
int Param::MAX_EPOCHS = 20000;

string Param::FEATURES_FILE = "../test_data/features_masked.bin";
string Param::LABELS_FILE = "../test_data/labels_masked.bin";
string Param::TRAIN_SUFFIXES = "../test_data/train_suffixes.txt";
string Param::TEST_SUFFIXES = "../test_data/test_suffixes.txt";

long Param::DIV_MAX_N = 100000;
long Param::PITER_BATCH_SIZE = 100;
long Param::PAR_THRES = 50;
long Param::NUM_THREADS = 20;

bool Param::PROFILER = true;
bool Param::DEBUG = false;
