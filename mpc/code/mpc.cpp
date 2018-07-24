#include "mpc.h"
#include "param.h"
#include "assert.h"
#include "primes.h"
#include <NTL/mat_ZZ_p.h>
#include <NTL/mat_ZZ.h>
#include <NTL/ZZ.h>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>

using namespace NTL;
using namespace std;

bool MPCEnv::Initialize(int pid, vector< pair<int, int> > &pairs) {
  tcout() << "Initializing MPC environment" << endl;

  /* Set base prime for the finite field */
  ZZ base_p = conv<ZZ>(Param::BASE_P.c_str());
  ZZ_p::init(base_p);

  this->pid = pid;
  this->clock_start = chrono::steady_clock::now();
  debug = false;

  if (!SetupChannels(pairs)) {
    tcout() << "MPCEnv::Initialize: failed to initialize communication channels" << endl;
    return false;
  }

  if (!SetupPRGs(pairs)) {
    tcout() << "MPCEnv::Initialize: failed to initialize PRGs" << endl;
    return false;
  }
  
  SetSeed(prg.find(pid)->second);
  cur_prg_pid = pid;

  primes.SetLength(3);
  primes[0] = ZZ_p::modulus();

  bool found1 = false;
  bool found2 = false;
  long thres1 = Param::NBIT_K / 2;
  long thres2 = ((long) ceil(sqrt((double) NumBits(ZZ_p::modulus())))) + 1;

  long ind = -1;
  long maxind = sizeof(PRIME_LIST) / sizeof(PRIME_LIST[0]);
  while ((!found1 || !found2) && ++ind < maxind) {
    long p = PRIME_LIST[ind];
    if (!found1 && p > thres1) {
      found1 = true;
      primes[1] = ZZ(p);
    }
    if (!found2 && p > thres2) {
      found2 = true;
      primes[2] = ZZ(p);
    }
  }

  if (!found1 || !found2) {
    tcout() << "Failed to find suitable small primes" << endl;
    return false;
  }

  tcout() << "Small base primes selected: " << primes[1] << "(>" << thres1 << ") and "
       << primes[2] << "(>" << thres2 << ")" << endl;

  ZZ_bytes.SetLength(primes.length());
  ZZ_bits.SetLength(primes.length());
  ZZ_per_buf.SetLength(primes.length());
  for (int i = 0; i < primes.length(); i++) {
    ZZ_bytes[i] = NumBytes(primes[i]);
    ZZ_bits[i] = NumBits(primes[i]);
    ZZ_per_buf[i] = (uint64_t) (Param::MPC_BUF_SIZE / ZZ_bytes[i]);
  }

  assert(ZZ_bytes[0] <= Param::MPC_BUF_SIZE); // buffer should contain at least one ZZ_p

  pstate.push("");

  buf = (unsigned char *) malloc(Param::MPC_BUF_SIZE + GCM_AUTH_TAG_LEN);
  if (buf == NULL) {
    tcout() << "Fail to allocate MPC buffer" << endl;
    exit(1);
  } else {
    tcout() << "Allocated MPC buffer of size " << Param::MPC_BUF_SIZE << endl;
  }

  tcout() << "Number of bytes per ZZ_p: " << ZZ_bytes[0] << endl;

  tcout() << "Setting up lookup tables" << endl;

  table_cache.SetLength(2);
  table_type_ZZ.SetLength(2);
  table_field_index.SetLength(2);
  lagrange_cache.SetLength(2);

  Mat<ZZ_p> table;

  // Table 0
  table.SetDims(1, 2);
  if (pid > 0) {
    table[0][0] = 1;
    table[0][1] = 0;
  }
  table_type_ZZ[0] = true;
  table_cache[0] = table;
  table_field_index[0] = 2;

  // Table 1
  int half_len = Param::NBIT_K / 2;
  table.SetDims(2, half_len + 1);
  if (pid > 0) {
    for (int i = 0; i < half_len + 1; i++) {
      if (i == 0) {
        table[0][i] = 1;
        table[1][i] = 1;
      } else {
        table[0][i] = table[0][i - 1] * 2;
        table[1][i] = table[1][i - 1] * 4;
      }
    }

  }
  table_type_ZZ[1] = true;
  table_cache[1] = table;
  table_field_index[1] = 1;

  /*
  // Table 3: parameters (intercept, slope) for piecewise-linear approximation of
  //          negative log-sigmoid function
  table.SetDims(2, 64);
  if (pid > 0) {
    ifstream ifs;
    ifs.open("sigmoid_approx.txt");
    if (!ifs.is_open()) {
      tcout() << "Error opening sigmoid_approx.txt" << endl;
      clear(table);
    }
    for (int i = 0; i < table.NumCols(); i++) {
      double intercept, slope;
      ifs >> intercept >> slope;

      ZZ_p fp_intercept, fp_slope;
      DoubleToFP(fp_intercept, intercept, Param::NBIT_K, Param::NBIT_F);
      DoubleToFP(fp_slope, slope, Param::NBIT_K, Param::NBIT_F);

      table[0][i] = fp_intercept;
      table[1][i] = fp_slope;
    }
    ifs.close();
  }
  table_type_ZZ[3] = false;
  table_cache[3] = table;
  */

  for (int cid = 0; cid < table_cache.length(); cid++) {
    long nrow = table_cache[cid].NumRows();
    long ncol = table_cache[cid].NumCols();
    bool index_by_ZZ = table_type_ZZ[cid];
    if (index_by_ZZ) {
      lagrange_cache[cid].SetDims(nrow, 2 * ncol);
    } else {
      lagrange_cache[cid].SetDims(nrow, ncol);
    }

    if (pid > 0) {
      tcout() << "Lagrange interpolation for Table " << cid << " ... ";
      for (int i = 0; i < nrow; i++) {
        Vec<long> x;
        Vec<ZZ_p> y;
        if (index_by_ZZ) {
          x.SetLength(2 * ncol);
          y.SetLength(2 * ncol);
        } else {
          x.SetLength(ncol);
          y.SetLength(ncol);
        }
        for (int j = 0; j < ncol; j++) {
          x[j] = j + 1;
          y[j] = table_cache[cid][i][j];
          if (index_by_ZZ) {
            x[j + ncol] = x[j] + conv<long>(primes[table_field_index[cid]]);
            y[j + ncol] = table_cache[cid][i][j];
          }
        }

        lagrange_interp(lagrange_cache[cid][i], x, y);
      }
      cout << "done" << endl;
    }
  }

  return true;
}

bool MPCEnv::SetupChannels(vector< pair<int, int> > &pairs) {
  for (int i = 0; i < pairs.size(); i++) {
    int p1 = pairs[i].first;
    int p2 = pairs[i].second;

    if (p1 != pid && p2 != pid) {
      continue;
    }

    int port = 8000;
    if (p1 == 0 && p2 == 1) {
      port = Param::PORT_P0_P1;
    } else if (p1 == 0 && p2 == 2) {
      port = Param::PORT_P0_P2;
    } else if (p1 == 1 && p2 == 2) {
      port = Param::PORT_P1_P2;
    } else if (p1 == 1 && p2 == 3) {
      port = Param::PORT_P1_P3;
    } else if (p1 == 2 && p2 == 3) {
      port = Param::PORT_P2_P3;
    }

    ostringstream oss;
    oss << Param::KEY_PATH << "P" << p1 << "_P" << p2 << ".key";
    string key_file = oss.str();

    int pother = p1 + p2 - pid;
    sockets.insert(map<int, CSocket>::value_type(pother, CSocket()));

    if (p1 == pid) {
      if (!OpenChannel(sockets[pother], port)) {
        tcout() << "Failed to connect with P" << pother << endl;
        return false;
      }
    } else {
      string ip_addr;
      if (pother == 0) {
        ip_addr = Param::IP_ADDR_P0;
      } else if (pother == 1) {
        ip_addr = Param::IP_ADDR_P1;
      } else if (pother == 2) {
        ip_addr = Param::IP_ADDR_P2;
      }

      if (!Connect(sockets[pother], ip_addr.c_str(), port)) {
        tcout() << "Failed to connect with P" << pother << endl;
        return false;
      }
    }

    if (!sockets[pother].SetKey(key_file)) {
      tcout() << "Failed to establish a secure channel with P" << pother << endl;
      return false;
    }

    tcout() << "Established a secure channel with P" << pother << endl;
  }

  tcout() << "Network setup complete" << endl;
  return true;
}

bool MPCEnv::SetupPRGs(vector< pair<int, int> > &pairs) {
  int key_len = NTL_PRG_KEYLEN; // from NTL
  unsigned char key[NTL_PRG_KEYLEN + GCM_AUTH_TAG_LEN];

  /* Internal PRG */
  int bytes = randread(key, key_len);
  if (bytes != key_len) {
    tcout() << "Failed to generate an internal PRG key" << endl;
    return false;
  }

  prg.insert(map<int, RandomStream>::value_type(pid, NewRandomStream(key)));
  
  /* Global PRG */
  ifstream ifs;
  string key_file = Param::KEY_PATH + "global.key";
  ifs.open(key_file.c_str(), ios::binary);
  if (!ifs.is_open()) {
    tcout() << "Failed to open global PRG key file: " << key_file << endl;
    return false;
  }

  ifs.read((char *)key, PRF_KEY_BYTES);
  if (ifs.gcount() != PRF_KEY_BYTES) {
    tcout() << "Failed to read " << PRF_KEY_BYTES << " bytes from global key file: " << key_file << endl;
    return false;
  }
  ifs.close();

  AESStream aes(key);
  aes.get(key, key_len);

  prg.insert(map<int, RandomStream>::value_type(-1, NewRandomStream(key)));

  /* Shared PRG (pairwise) */
  for (int i = 0; i < pairs.size(); i++) {
    int p1 = pairs[i].first;
    int p2 = pairs[i].second;

    if (p1 != pid && p2 != pid) {
      continue;
    }

    int pother = p1 + p2 - pid;

    if (p1 == pid) {
      bytes = randread(key, key_len);
      if (bytes != key_len) {
        tcout() << "Failed to generate a shared PRG key" << endl;
        return false;
      }

      prg.insert(map<int, RandomStream>::value_type(pother, NewRandomStream(key)));
      sockets[pother].SendSecure(key, key_len);
    } else {
      sockets[pother].ReceiveSecure(key, key_len);
      prg.insert(map<int, RandomStream>::value_type(pother, NewRandomStream(key)));
    }

    tcout() << "Shared PRG with P" << pother << " initialized" << endl;
  }

  tcout() << "PRG setup complete" << endl;
  return true;
}


void MPCEnv::CleanUp() {
  tcout() << "Closing sockets ... ";
  for (map<int, CSocket>::iterator it = sockets.begin(); it != sockets.end(); ++it) {
    CloseChannel(it->second);
  }
  cout << "done." << endl;
}

void MPCEnv::ProfilerResetTimer() {
  if (!Param::PROFILER) return;

  vector<uint64_t> stat(5, 0);

  chrono::time_point<chrono::steady_clock> clock_end = chrono::steady_clock::now();

  stat[0] = chrono::duration_cast<chrono::milliseconds>(clock_end - clock_start).count();

  int ind = 1;
  for (int p = 0; p < 3; p++) {
    if (p == pid) continue;
    map<int, CSocket>::iterator it = sockets.find(p);
    stat[ind] = it->second.GetBytesSent();
    stat[ind+1] = it->second.GetBytesReceived();
    it->second.ResetStats();
    ind += 2;
  }

  string state = pstate.top();
  if (state != "") {
    map<string, int>::iterator it = ptable_index.find(state);
    if (it != ptable_index.end()) {
      for (int i = 0; i < stat.size(); i++) {
        ptable[it->second].second[i] += stat[i];
      }
    } else {
      ptable.push_back(make_pair(state, stat));
      ptable_index[state] = ptable.size() - 1;
    }
  }

  clock_start = clock_end;
}

void MPCEnv::ProfilerPushState(string desc) {
  if (!Param::PROFILER) return;

  assert(desc != "");

  ProfilerResetTimer();

  string full_desc;
  if (pstate.top() == "") {
    full_desc = desc;
  } else {
    full_desc = pstate.top() + ">" + desc;
  }

  pstate.push(full_desc);
}

void MPCEnv::ProfilerPopState(bool write) {
  if (!Param::PROFILER) return;

  assert(pstate.top() != "");

  ProfilerResetTimer();

  pstate.pop();

  if (write) {
    ProfilerWriteToFile();
  }
}

void MPCEnv::ProfilerWriteToFile() {
  if (!Param::PROFILER) return;

  ProfilerResetTimer();
  
  // Open log file
  logfs.open(Param::LOG_FILE.c_str());
  if (!logfs.is_open()) {
    tcout() << "Fail to open the log file: " << Param::LOG_FILE << endl;
    exit(1);
  }

  ostringstream oss;

  // log file header
  oss << "Desc\tTime(ms)";
  for (int p = 0; p < 3; p++) {
    if (p == pid) continue;
    oss << "\tTo_" << p;
    oss << "\tFrom_" << p;
  }
  oss << endl;

  for (int i = 0; i < ptable.size(); i++) {
    oss << ptable[i].first;
    for (int j = 0; j < ptable[i].second.size(); j++) {
      oss << "\t" << ptable[i].second[j];
    }
    oss << endl;
  }

  logfs << oss.str();
  logfs.flush();
  logfs.close();
}

/*
void MPCEnv::LogisticRegression(ZZ_p& nll, ZZ_p& b0, Vec<ZZ_p>& b, Mat<ZZ_p>& x, Vec<ZZ_p>& y) {
  size_t n = x.NumCols();
  size_t p = x.NumRows();
  assert(y.length() == n);

  // Initialize
  b0 = 0;
  b.SetLength(p);
  clear(b);

  Mat<ZZ_p> xr, xm;
  BeaverPartition(xr, xm, x);

  Vec<ZZ_p> yneg = -y;
  if (pid == 1) {
    for (int i = 0; i < n; i++) {
      yneg[i] += 1;
    }
  }

  Vec<ZZ_p> yneg_r, yneg_m;
  BeaverPartition(yneg_r, yneg_m, yneg);

  ZZ_p fp_memory; 
  DoubleToFP(fp_memory, 0.9, Param::NBIT_K, Param::NBIT_F);

  ZZ_p fp_n_inv; 
  DoubleToFP(fp_n_inv, 1 / ((double) n), Param::NBIT_K, Param::NBIT_F);

  ZZ_p fp_one; 
  DoubleToFP(fp_one, 1, Param::NBIT_K, Param::NBIT_F);

  ZZ_p step0(0);
  Vec<ZZ_p> step;
  step.SetLength(p);
  clear(step);

  // Gradient descent (with momentum)
  for (int it = 0; it < 100; it++) {
    tcout() << "Iter " << it << ": ";

    Mat<ZZ_p> br, bm;
    br.SetDims(1, p);
    bm.SetDims(1, p);
    BeaverPartition(br[0], bm[0], b);

    Mat<ZZ_p> h;
    h.SetDims(1, n);
    clear(h);

    BeaverMultMat(h, br, bm, xr, xm);
    BeaverReconstruct(h);
    Trunc(h);

    for (int i = 0; i < n; i++) {
      h[0][i] += b0;
    }

    Vec<ZZ_p> s, s_grad;
    NegLogSigmoid(s, s_grad, h[0]);

    // Compute gradient
    ZZ_p d0(0);

    for (int i = 0; i < n; i++) {
      s_grad[i] += yneg[i] * fp_one;
      d0 += s_grad[i];
    }

    Vec<ZZ_p> s_grad_r, s_grad_m;
    BeaverPartition(s_grad_r, s_grad_m, s_grad);

    Vec<ZZ_p> d;
    d.SetLength(p);
    clear(d);

    BeaverMult(d, xr, xm, s_grad_r, s_grad_m);
    BeaverReconstruct(d);
    Trunc(d);

    step0 = step0 * fp_memory - d0 * fp_n_inv;
    for (int i = 0; i < p; i++) {
      step[i] = step[i] * fp_memory - d[i] * fp_n_inv;
    }

    Trunc(step0);
    Trunc(step);

    b0 = b0 + step0;
    for (int i = 0; i < p; i++) {
      b[i] += step[i];
    }

    Vec<ZZ_p> bvec;
    bvec.SetLength(p + 1);
    bvec[0] = b0;
    for (int i = 0; i < p; i++) {
      bvec[i + 1] = b[i];
    }
    PrintFP(bvec);
  }

  // Compute final NLL
  Mat<ZZ_p> br, bm;
  br.SetDims(1, p);
  bm.SetDims(1, p);
  BeaverPartition(br[0], bm[0], b);

  Mat<ZZ_p> h;
  h.SetDims(1, n);
  clear(h);

  BeaverMultMat(h, br, bm, xr, xm);
  BeaverReconstruct(h);
  Trunc(h);

  for (int i = 0; i < n; i++) {
    h[0][i] += b0;
  }

  Mat<ZZ_p> hr, hm;
  BeaverPartition(hr, hm, h);

  Vec<ZZ_p> yh;
  yh.SetLength(1);
  clear(yh);

  BeaverMult(yh, hr, hm, yneg_r, yneg_m);
  BeaverReconstruct(yh);

  nll = yh[0];

  Vec<ZZ_p> s, s_grad;
  NegLogSigmoid(s, s_grad, h[0]);

  for (int i = 0; i < n; i++) {
    nll += s[i];
  }

  tcout() << "NLL: ";
  PrintFP(nll);
}

void MPCEnv::NegLogSigmoid(Vec<ZZ_p>& b, Vec<ZZ_p>& b_grad, Vec<ZZ_p>& a) {
  size_t n = a.length();

  int depth = 6;

  Vec<ZZ_p> cur = a; // copy

  Vec<ZZ_p> a_ind;
  a_ind.SetLength(a.length());
  clear(a_ind);

  double step = 4;

  for (int i = 0; i < depth; i++) {
    Vec<ZZ_p> cur_sign;
    IsPositive(cur_sign, cur);

    ZZ_p index_step(1 << (depth - 1 - i));

    for (int j = 0; j < n; j++) {
      a_ind[j] += cur_sign[j] * index_step;
    }

    cur_sign *= 2;
    if (pid == 1) {
      for (int j = 0; j < n; j++) {
        cur_sign[j] -= 1;
      }
    }

    ZZ_p step_fp;
    DoubleToFP(step_fp, step, Param::NBIT_K, Param::NBIT_F);

    for (int j = 0; j < n; j++) {
      cur[j] -= step_fp * cur_sign[j];
    }

    step /= 2;
  }

  // Make indices 1-based
  if (pid == 1) {
    for (int j = 0; j < n; j++) {
      a_ind[j]++;
    }
  }

  // Fetch piecewise linear approx parameters
  Mat<ZZ_p> param;
  TableLookup(param, a_ind, 3);

  MultElem(b, param[1], a);
  Trunc(b);

  if (pid > 0) {
    for (int j = 0; j < n; j++) {
      b[j] += param[0][j];
    }
  }

  b_grad = param[1];
}
*/

void MPCEnv::InnerProd(Vec<ZZ_p>& c, Mat<ZZ_p>& a) {
  if (debug) tcout() << "InnerProd: " << a.NumRows() << ", " << a.NumCols() << endl;

  Mat<ZZ_p> ar, am;
  BeaverPartition(ar, am, a);

  Init(c, a.NumRows());
  for (int i = 0; i < a.NumRows(); i++) {
    BeaverInnerProd(c[i], ar[i], am[i]);
  }

  BeaverReconstruct(c);
}
void MPCEnv::InnerProd(ZZ_p& c, Vec<ZZ_p>& a) {
  if (debug) tcout() << "InnerProd: " << a.length() << endl;

  Vec<ZZ_p> ar, am;
  BeaverPartition(ar, am, a);
  BeaverInnerProd(c, ar, am);
  BeaverReconstruct(c);
}

void MPCEnv::Householder(Vec<ZZ_p>& v, Vec<ZZ_p>& x) {
  if (debug) tcout() << "Householder: " << x.length() << endl;

  int n = x.length();

  Vec<ZZ_p> xr, xm;
  BeaverPartition(xr, xm, x);

  Vec<ZZ_p> xdot;
  Init(xdot, 1);
  BeaverInnerProd(xdot[0], xr, xm);
  BeaverReconstruct(xdot);
  Trunc(xdot);

  Vec<ZZ_p> xnorm, dummy;
  FPSqrt(xnorm, dummy, xdot);

  Vec<ZZ_p> x1;
  x1.SetLength(1);
  x1[0] = x[0];

  Vec<ZZ_p> x1sign;
  IsPositive(x1sign, x1);

  x1sign *= 2;
  if (pid == 1) {
    x1sign[0] -= 1;
  }

  Vec<ZZ_p> shift;
  MultElem(shift, xnorm, x1sign);

  ZZ_p sr, sm;
  BeaverPartition(sr, sm, shift[0]);

  ZZ_p dot_shift(0);
  BeaverMult(dot_shift, xr[0], xm[0], sr, sm);
  BeaverReconstruct(dot_shift);
  Trunc(dot_shift);

  Vec<ZZ_p> vdot;
  vdot.SetLength(1);
  if (pid > 0) {
    vdot[0] = 2 * (xdot[0] + dot_shift);
  }

  Vec<ZZ_p> vnorm_inv;
  FPSqrt(dummy, vnorm_inv, vdot);

  ZZ_p invr, invm;
  BeaverPartition(invr, invm, vnorm_inv[0]);
 
  Vec<ZZ_p> vr, vm;
  if (pid > 0) {
    vr = xr;
    vr[0] += sr;
  } else {
    vr.SetLength(n);
  }
  vm = xm;
  vm[0] += sm;

  Init(v, n);
  BeaverMult(v, vr, vm, invr, invm);
  BeaverReconstruct(v);
  Trunc(v);
}

void MPCEnv::QRFactSquare(Mat<ZZ_p>& Q, Mat<ZZ_p>& R, Mat<ZZ_p>& A) {
  if (debug) tcout() << "QRFactSquare: " << A.NumRows() << ", " << A.NumCols() << endl;

  assert(A.NumRows() == A.NumCols());

  int n = A.NumRows();
  R.SetDims(n, n);
  if (pid > 0) {
    clear(R);
  }

  Mat<ZZ_p> Ap;
  if (pid == 0) {
    Ap.SetDims(n, n);
  } else {
    Ap = A;
  }

  ZZ_p one;
  DoubleToFP(one, 1, Param::NBIT_K, Param::NBIT_F);

  for (int i = 0; i < n - 1; i++) {
    Mat<ZZ_p> v;
    v.SetDims(1, Ap.NumCols());
    Householder(v[0], Ap[0]);

    Mat<ZZ_p> vt;
    if (pid == 0) {
      vt.SetDims(Ap.NumCols(), 1);
    } else {
      transpose(vt, v);
    }

    Mat<ZZ_p> P;
    MultMat(P, vt, v);
    Trunc(P);
    if (pid > 0) {
      P *= -2;
      if (pid == 1) {
        for (int j = 0; j < P.NumCols(); j++) {
          P[j][j] += one;
        }
      }
    }

    Mat<ZZ_p> B;
    if (i == 0) {
      Q = P;
      MultMat(B, Ap, P);
      Trunc(B);
    } else {
      Mat<ZZ_p> Qsub;
      Qsub.SetDims(n - i, n);
      if (pid > 0) {
        for (int j = 0; j < n - i; j++) {
          Qsub[j] = Q[j+i];
        }
      }

      Vec< Mat<ZZ_p> > left;
      Vec< Mat<ZZ_p> > right;
      left.SetLength(2);
      right.SetLength(2);
      left[0] = P;
      right[0] = Qsub;
      left[1] = Ap;
      right[1] = P;

      Vec< Mat<ZZ_p> > prod;
      MultMatParallel(prod, left, right);
      // TODO: parallelize Trunc
      Trunc(prod[0]);
      Trunc(prod[1]);

      if (pid > 0) {
        for (int j = 0; j < n - i; j++) {
          Q[j+i] = prod[0][j];
        }
        B = prod[1];
      } else {
        B.SetDims(n - i, n - i);
      }
    }

    if (pid > 0) {
      for (int j = 0; j < n - i; j++) {
        R[i+j][i] = B[j][0];
      }
      if (i == n - 2) {
        R[n-1][n-1] = B[1][1];
      }

      Ap.SetDims(n - i - 1, n - i - 1);
      for (int j = 0; j < n - i - 1; j++) {
        for (int k = 0; k < n - i - 1; k++) {
          Ap[j][k] = B[j+1][k+1];
        }
      }
    } else {
      Ap.SetDims(n - i - 1, n - i - 1);
    }
  }
}

void MPCEnv::OrthonormalBasis(Mat<ZZ_p>& Q, Mat<ZZ_p>& A) {
  if (debug) tcout() << "OrthonormalBasis: " << A.NumRows() << ", " << A.NumCols() << endl;

  assert(A.NumCols() >= A.NumRows());

  int c = A.NumRows();
  int n = A.NumCols();

  Vec< Vec<ZZ_p> > v_list;
  v_list.SetLength(c);

  Mat<ZZ_p> Ap;
  if (pid == 0) {
    Ap.SetDims(c, n);
  } else {
    Ap = A;
  }

  ZZ_p one;
  DoubleToFP(one, 1, Param::NBIT_K, Param::NBIT_F);

  for (int i = 0; i < c; i++) {
    Mat<ZZ_p> v;
    v.SetDims(1, Ap.NumCols());
    Householder(v[0], Ap[0]);

    if (pid == 0) {
      v_list[i].SetLength(Ap.NumCols());
    } else {
      v_list[i] = v[0];
    }

    Mat<ZZ_p> vt;
    if (pid == 0) {
      vt.SetDims(Ap.NumCols(), 1);
    } else {
      transpose(vt, v);
    }

    Mat<ZZ_p> Apv;
    MultMat(Apv, Ap, vt);
    Trunc(Apv);

    Mat<ZZ_p> B;
    MultMat(B, Apv, v);
    Trunc(B);
    if (pid > 0) {
      B *= -2;
      B += Ap;
    }

    Ap.SetDims(B.NumRows() - 1, B.NumCols() - 1);
    if (pid > 0) {
      for (int j = 0; j < B.NumRows() - 1; j++) {
        for (int k = 0; k < B.NumCols() - 1; k++) {
          Ap[j][k] = B[j+1][k+1];
        }
      }
    }
  }

  Q.SetDims(c, n);
  if (pid > 0) {
    clear(Q);
    if (pid == 1) {
      for (int i = 0; i < c; i++) {
        Q[i][i] = one;
      }
    }
  }

  for (int i = c - 1; i >= 0; i--) {
    Mat<ZZ_p> v;
    v.SetDims(1, v_list[i].length());
    if (pid > 0) {
      v[0] = v_list[i];
    }

    Mat<ZZ_p> vt;
    if (pid == 0) {
      vt.SetDims(v.NumCols(), 1);
    } else {
      transpose(vt, v);
    }

    Mat<ZZ_p> Qsub;
    Qsub.SetDims(c, n - i);
    if (pid > 0) {
      for (int j = 0; j < c; j++) {
        for (int k = 0; k < n - i; k++) {
          Qsub[j][k] = Q[j][k+i];
        }
      }
    }

    Mat<ZZ_p> Qv;
    MultMat(Qv, Qsub, vt);
    Trunc(Qv);

    Mat<ZZ_p> Qvv;
    MultMat(Qvv, Qv, v);
    Trunc(Qvv);
    if (pid > 0) {
      Qvv *= -2;
    }

    if (pid > 0) {
      for (int j = 0; j < c; j++) {
        for (int k = 0; k < n - i; k++) {
          Q[j][k+i] += Qvv[j][k];
        }
      }
    }
  }
}

void MPCEnv::Tridiag(Mat<ZZ_p>& T, Mat<ZZ_p>& Q, Mat<ZZ_p>& A) {
  if (debug) tcout() << "Tridiag: " << A.NumRows() << ", " << A.NumCols() << endl;

  assert(A.NumRows() == A.NumCols());
  assert(A.NumRows() > 2);

  int n = A.NumRows();

  ZZ_p one;
  DoubleToFP(one, 1, Param::NBIT_K, Param::NBIT_F);

  Q.SetDims(n, n);
  T.SetDims(n, n);
  if (pid > 0) {
    clear(Q);
    clear(T);
    if (pid == 1) {
      for (int i = 0; i < n; i++) {
        Q[i][i] = one;
      }
    }
  }

  Mat<ZZ_p> Ap;
  if (pid == 0) {
    Ap.SetDims(n, n);
  } else {
    Ap = A;
  }

  for (int i = 0; i < n - 2; i++) {
    Vec<ZZ_p> x;
    x.SetLength(Ap.NumCols() - 1);
    if (pid > 0) {
      for (int j = 0; j < Ap.NumCols() - 1; j++) {
        x[j] = Ap[0][j+1];
      }
    }

    Mat<ZZ_p> v;
    v.SetDims(1, x.length());
    Householder(v[0], x);

    Mat<ZZ_p> vt;
    if (pid == 0) {
      vt.SetDims(x.length(), 1);
    } else {
      transpose(vt, v);
    }

    Mat<ZZ_p> vv;
    MultMat(vv, vt, v);
    Trunc(vv);

    Mat<ZZ_p> P;
    P.SetDims(Ap.NumCols(), Ap.NumCols());
    if (pid > 0) {
      P[0][0] = (pid == 1) ? one : ZZ_p(0);
      for (int j = 1; j < Ap.NumCols(); j++) {
        for (int k = 1; k < Ap.NumCols(); k++) {
          P[j][k] = -2 * vv[j-1][k-1];
          if (pid == 1 && j == k) {
            P[j][k] += one;
          }
        }
      }
    }

    // TODO: parallelize? (minor improvement)
    Mat<ZZ_p> PAp;
    MultMat(PAp, P, Ap);
    Trunc(PAp);

    Mat<ZZ_p> B;
    MultMat(B, PAp, P);
    Trunc(B);

    Mat<ZZ_p> Qsub;
    Qsub.SetDims(n, n - i);
    if (pid > 0) {
      for (int j = 0; j < n; j++) {
        for (int k = 0; k < n - i; k++) {
          Qsub[j][k] = Q[j][k+i];
        }
      }
    }

    MultMat(Qsub, Qsub, P);
    Trunc(Qsub);
    if (pid > 0) {
      for (int j = 0; j < n; j++) {
        for (int k = 0; k < n - i; k++) {
          Q[j][k+i] = Qsub[j][k];
        }
      }
    }

    if (pid > 0) {
      T[i][i] = B[0][0];
      T[i+1][i] = B[1][0];
      T[i][i+1] = B[0][1];
      if (i == n - 3) {
        T[i+1][i+1] = B[1][1];
        T[i+1][i+2] = B[1][2];
        T[i+2][i+1] = B[2][1];
        T[i+2][i+2] = B[2][2];
      }
    }

    Ap.SetDims(B.NumRows() - 1, B.NumCols() - 1);
    if (pid > 0) {
      for (int j = 0; j < B.NumRows() - 1; j++) {
        for (int k = 0; k < B.NumCols() - 1; k++) {
          Ap[j][k] = B[j+1][k+1];
        }
      }
    }
  }
}

void MPCEnv::EigenDecomp(Mat<ZZ_p>& V, Vec<ZZ_p>& L, Mat<ZZ_p>& A) {
  if (debug) tcout() << "EigenDecomp: " << A.NumRows() << ", " << A.NumCols() << endl;

  assert(A.NumRows() == A.NumCols());
  int n = A.NumRows();

  L.SetLength(n);
  clear(L);

  Mat<ZZ_p> Ap, Q;
  Tridiag(Ap, Q, A);

  if (pid == 0) {
    V.SetDims(n, n);
  } else {
    transpose(V, Q);
  }

  for (int i = n - 1; i >= 1; i--) {
    tcout() << "EigenDecomp: " << i << "-th eigenvalue" << endl;
    for (int it = 0; it < Param::ITER_PER_EVAL; it++) {
      ZZ_p shift = Ap[i][i];
      if (pid > 0) {
        for (int j = 0; j < Ap.NumCols(); j++) {
          Ap[j][j] -= shift;
        }
      }

      Mat<ZZ_p> R;
      QRFactSquare(Q, R, Ap);

      MultMat(Ap, Q, R);
      Trunc(Ap);

      if (pid > 0) {
        for (int j = 0; j < Ap.NumCols(); j++) {
          Ap[j][j] += shift;
        }
      }

      Mat<ZZ_p> Vsub;
      Vsub.SetDims(i + 1, n);
      if (pid > 0) {
        for (int j = 0; j < i + 1; j++) {
          Vsub[j] = V[j];
        }
      }

      MultMat(Vsub, Q, Vsub);
      Trunc(Vsub);

      if (pid > 0) {
        for (int j = 0; j < i + 1; j++) {
          V[j] = Vsub[j];
        }
      }
    }

    L[i] = Ap[i][i];
    if (i == 1) {
      L[0] = Ap[0][0];
    }

    Mat<ZZ_p> Ap_copy = Ap;
    Ap.SetDims(i, i);
    if (pid > 0) {
      for (int j = 0; j < i; j++) {
        for (int k = 0; k < i; k++) {
          Ap[j][k] = Ap_copy[j][k];
        }
      }
    }
  }
  tcout() << "EigenDecomp: complete" << endl;
}

void MPCEnv::LessThanBitsAux(Vec<ZZ>& c, Mat<ZZ>& a, Mat<ZZ>& b, int public_flag, int fid) {
  if (debug) tcout() << "LessThanBitsAux: " << a.NumRows() << ", " << a.NumCols() << endl;
 
  assert(a.NumRows() == b.NumRows());
  assert(a.NumCols() == b.NumCols());
 
  int n = a.NumRows();
  int L = a.NumCols();
 
  /* Calculate XOR */
  Mat<ZZ> x;
  x.SetDims(n, L);
 
  if (public_flag == 0) {
    MultElem(x, a, b, fid);
    if (pid > 0) {
      x = a + b - 2 * x;
      Mod(x, fid);
    }
  } else if (pid > 0) {
    mul_elem(x, a, b);
    x = a + b - 2 * x;
    if (pid == 2) {
      x -= (public_flag == 1) ? a : b;
    }
    Mod(x, fid);
  }
 
  Mat<ZZ> f;
  PrefixOr(f, x, fid);
  x.kill();
 
  if (pid > 0) {
    for (int i = 0; i < n; i++) {
      for (int j = L - 1; j >= 1; j--) {
        f[i][j] -= f[i][j - 1];
      }
    }
    Mod(f, fid);
  }
 
  if (public_flag == 2) {
    c.SetLength(n);
 
    if (pid > 0) {
      for (int i = 0; i < n; i++) {
        c[i] = 0;
        for (int j = 0; j < L; j++) {
          c[i] += f[i][j] * b[i][j];
        }
      }
      Mod(c, fid);
    }
  } else {
    //TODO: optimize
    Vec< Mat<ZZ> > f_arr, b_arr;
    f_arr.SetLength(n);
    b_arr.SetLength(n);
    for (int i = 0; i < n; i++) {
      f_arr[i].SetDims(1, L);
      b_arr[i].SetDims(L, 1);
    }
 
    if (pid > 0) {
      for (int i = 0; i < n; i++) {
        f_arr[i][0] = f[i];
        for (int j = 0; j < L; j++) {
          b_arr[i][j][0] = b[i][j];
        }
      }
    }
 
    Vec< Mat<ZZ> > c_arr;
    MultMatParallel(c_arr, f_arr, b_arr, fid);
 
    c.SetLength(n);
    if (pid > 0) {
      for (int i = 0; i < n; i++) {
        c[i] = c_arr[i][0][0];
      }
    }
  }
}

void MPCEnv::LessThan(Vec<ZZ_p>& c, Vec<ZZ_p>& a, Vec<ZZ_p>& b) {
  Vec<ZZ_p> a_cpy;
  a_cpy.SetLength(a.length());
  if (pid > 0) {
    for (int i = 0; i < a.length(); i++) {
      a_cpy[i] = a[i] - b[i];
    }
  }

  // a - b >= 0?
  IsPositive(c, a_cpy);
  FlipBit(c);
}

void MPCEnv::LessThanPublic(Vec<ZZ_p>& c, Vec<ZZ_p>& a, ZZ_p bpub) {
  Vec<ZZ_p> a_cpy;
  a_cpy.SetLength(a.length());
  if (pid > 0) {
    for (int i = 0; i < a.length(); i++) {
      a_cpy[i] = a[i];
      if (pid == 1) {
        a_cpy[i] -= bpub;
      }
    }
  }

  // a - b >= 0?
  IsPositive(c, a_cpy);
  FlipBit(c);
}

void MPCEnv::IsPositive(Mat<ZZ_p>& b, Mat<ZZ_p>& a) {
  b.SetDims(a.NumRows(), a.NumCols());
  for (int i = 0; i < a.NumRows(); i++) {
    if ((pid) == 69)
      tcout() << "Row " << i << endl;
    IsPositive(b[i], a[i]);
  }
}

// Failure probability of 1 / BASE_P
// Base field index 2
void MPCEnv::IsPositive(Vec<ZZ_p>& b, Vec<ZZ_p>& a) {
  if (debug) tcout() << "IsPositive: " << a.length() << endl; 

  int n = a.length();
  int nbits = ZZ_bits[0];
  int fid = 2;

  if ((pid) == 69)
    tcout() << "Transfering data" << endl;
  Vec<ZZ_p> r;
  Mat<ZZ> r_bits;
  if (pid == 0) {
    RandVec(r, n);
    NumToBits(r_bits, r, nbits);

    SwitchSeed(1);
    Vec<ZZ_p> r_mask;
    Mat<ZZ> r_bits_mask;
    RandVec(r_mask, n);
    RandMat(r_bits_mask, n, nbits, fid);
    RestoreSeed();

    r -= r_mask;
    r_bits -= r_bits_mask;
    Mod(r_bits, fid);

    SendVec(r, 2);
    SendMat(r_bits, 2, fid);
  } else if (pid == 2) {
    ReceiveVec(r, 0, n);
    ReceiveMat(r_bits, 0, n, nbits, fid);
  } else {
    SwitchSeed(0);
    RandVec(r, n);
    RandMat(r_bits, n, nbits, fid);
    RestoreSeed();
  }

  if ((pid) == 69)
    tcout() << "Compute c" << endl;
  Vec<ZZ_p> c;
  if (pid == 0) {
    c.SetLength(n);
  } else {
    c = 2 * a + r;
  }

  if ((pid) == 69)
    tcout() << "Reveal c" << endl;
  RevealSym(c);

  if ((pid) == 69)
    tcout() << "c to c_bits" << endl;
  Mat<ZZ> c_bits;
  if (pid == 0) {
    c_bits.SetDims(n, nbits);
  } else {
    NumToBits(c_bits, c, nbits);
  }

  if ((pid) == 69)
    tcout() << "Less than bits public" << endl;
  // Incorrect result if r = 0, which happens with probaility 1 / BASE_P
  Vec<ZZ> no_overflow;
  LessThanBitsPublic(no_overflow, r_bits, c_bits, fid);

  if ((pid) == 69)
    tcout() << "Compute c_xor_r" << endl;
  Vec<ZZ> c_xor_r;
  c_xor_r.SetLength(n);
  if (pid > 0) {
    for (int i = 0; i < n; i++) {
      c_xor_r[i] = r_bits[i][nbits-1] - 2 * c_bits[i][nbits-1] * r_bits[i][nbits-1];
      if (pid == 1) {
        c_xor_r[i] += c_bits[i][nbits-1];
      }
    }
    Mod(c_xor_r, fid);
  }

  if ((pid) == 69)
    tcout() << "Compute lsb" << endl;
  Vec<ZZ> lsb;
  MultElem(lsb, c_xor_r, no_overflow, fid);
  if (pid > 0) {
    lsb *= 2;
    lsb -= no_overflow + c_xor_r;
    if (pid == 1) {
      for (int i = 0; i < n; i++) {
        lsb[i] += 1;
      }
    }
    Mod(lsb, fid);
  }

  if ((pid) == 69)
    tcout() << "Mod lsb" << endl;
  // 0, 1 -> 1, 2
  if (pid == 1) {
    for (int i = 0; i < n; i++) {
      lsb[i] += 1;
    }
    Mod(lsb, fid);
  }

  if ((pid) == 69)
    tcout() << "Table lookup" << endl;
  Mat<ZZ_p> b_mat;
  TableLookup(b_mat, lsb, 0, fid);

  b = b_mat[0];
  if ((pid) == 69)
    tcout() << "Done" << endl;
}

void MPCEnv::FlipBit(Vec<ZZ_p>& b, Vec<ZZ_p>& a) {
  if (debug) tcout() << "FlipBit: " << a.length() << endl;
  if (pid == 0) {
    b.SetLength(a.length());
  } else {
    b = -a;
  }

  if (pid == 1) {
    for (int i = 0; i < b.length(); i++) {
      b[i] += 1;
    }
  }
}

// Assumes Param::NBIT_K - NBIT_F is even
void MPCEnv::FPSqrt(Vec<ZZ_p>& b, Vec<ZZ_p>& b_inv, Vec<ZZ_p>& a) {
  if (debug) tcout() << "FPSqrt: " << a.length() << endl; 

  int n = a.length();

  if (n > Param::DIV_MAX_N) {
    int nbatch = ceil(n / ((double) Param::DIV_MAX_N));
    b.SetLength(n);
    b_inv.SetLength(n);
    for (int i = 0; i < nbatch; i++) {
      tcout() << "FPSqrt on large vector: " << i + 1 << "/" << nbatch << endl;
      int start = Param::DIV_MAX_N * i;
      int end = start + Param::DIV_MAX_N;
      if (end > n) {
        end = n;
      }
      int batch_size = end - start;
      Vec<ZZ_p> a_copy;
      a_copy.SetLength(batch_size);
      for (int j = 0; j < batch_size; j++) {
        a_copy[j] = a[start + j];
      }
      Vec<ZZ_p> b_copy, b_inv_copy;
      FPSqrt(b_copy, b_inv_copy, a_copy);
      for (int j = 0; j < batch_size; j++) {
        b[start + j] = b_copy[j];
        b_inv[start + j] = b_inv_copy[j];
      }
    }
    return;
  }

  // TODO: Currently using the same # iter as division -- possibly need to update
  int niter = 2 * ceil(log2(((double) Param::NBIT_K) / 3.5));

  /* Initial approximation: 1 / sqrt(a_scaled) ~= 2.9581 - 4 * a_scaled + 2 * a_scaled^2 */
  Vec<ZZ_p> s, s_sqrt;
  NormalizerEvenExp(s, s_sqrt, a);

  Vec<ZZ_p> a_scaled;
  MultElem(a_scaled, a, s);
  Trunc(a_scaled, Param::NBIT_K, Param::NBIT_K - Param::NBIT_F);

  Vec<ZZ_p> a_scaled_sq;
  MultElem(a_scaled_sq, a_scaled, a_scaled);
  Trunc(a_scaled_sq);

  Vec<ZZ_p> scaled_est;
  if (pid == 0) {
    scaled_est.SetLength(n);
  } else {
    scaled_est = - 4 * a_scaled + 2 * a_scaled_sq;
    if (pid == 1) {
      ZZ_p coeff;
      DoubleToFP(coeff, 2.9581, Param::NBIT_K, Param::NBIT_F);
      for (int i = 0; i < n; i++) {
        scaled_est[i] += coeff;
      }
    }
  }

  Vec< Mat<ZZ_p> > h_and_g;
  h_and_g.SetLength(2);
  h_and_g[0].SetDims(1, n);
  h_and_g[1].SetDims(1, n);

  MultElem(h_and_g[0][0], scaled_est, s_sqrt);
  // Our scaled initial approximation (scaled_est) has bit length <= NBIT_F + 2
  // and s_sqrt is at most NBIT_K/2 bits, so their product is at most NBIT_K/2 + NBIT_F + 2
  Trunc(h_and_g[0], Param::NBIT_K/2 + Param::NBIT_F + 2, ((Param::NBIT_K - Param::NBIT_F) / 2) + 1); 

  h_and_g[1][0] = h_and_g[0][0] * 2;
  MultElem(h_and_g[1][0], h_and_g[1][0], a);
  Trunc(h_and_g[1]);

  ZZ_p onepointfive;
  DoubleToFP(onepointfive, 1.5, Param::NBIT_K, Param::NBIT_F);

  for (int it = 0; it < niter; it++) {
    Mat<ZZ_p> r;
    MultElem(r, h_and_g[0], h_and_g[1]);
    Trunc(r);
    r = -r;
    if (pid == 1) {
      for (int i = 0; i < n; i++) {
        r[0][i] += onepointfive;
      }
    }

    Vec< Mat<ZZ_p> > r_dup;
    r_dup.SetLength(2);
    r_dup[0] = r;
    r_dup[1] = r;

    MultElemParallel(h_and_g, h_and_g, r_dup);
    // TODO: write a version of Trunc with parallel processing
    Trunc(h_and_g[0]);
    Trunc(h_and_g[1]);
  }

  b_inv = 2 * h_and_g[0][0];
  b = h_and_g[1][0];
}

void MPCEnv::FPDiv(Vec<ZZ_p>& c, Vec<ZZ_p>& a, Vec<ZZ_p>& b) {
  if (debug) tcout() << "FPDiv: " << a.length() << endl; 

  assert(a.length() == b.length());

  int n = a.length();
  if (n > Param::DIV_MAX_N) {
    int nbatch = ceil(n / ((double) Param::DIV_MAX_N));
    c.SetLength(n);
    for (int i = 0; i < nbatch; i++) {
      int start = Param::DIV_MAX_N * i;
      int end = start + Param::DIV_MAX_N;
      if (end > n) {
        end = n;
      }
      int batch_size = end - start;

      tcout() << "FPDiv on large vector: " << i + 1 << "/" << nbatch << ", n = " << batch_size << endl;

      Vec<ZZ_p> a_copy, b_copy;
      a_copy.SetLength(batch_size);
      b_copy.SetLength(batch_size);
      for (int j = 0; j < batch_size; j++) {
        a_copy[j] = a[start + j];
        b_copy[j] = b[start + j];
      }
      Vec<ZZ_p> c_copy;
      FPDiv(c_copy, a_copy, b_copy);
      for (int j = 0; j < batch_size; j++) {
        c[start + j] = c_copy[j];
      }
    }
    return;
  }

  int niter = 2 * ceil(log2(((double) Param::NBIT_K) / 3.5)) + 1;

  /* Initial approximation: 1 / x_scaled ~= 5.9430 - 10 * x_scaled + 5 * x_scaled^2 */
  Vec<ZZ_p> s, s_sqrt;
  NormalizerEvenExp(s, s_sqrt, b);

  Vec<ZZ_p> b_scaled;
  MultElem(b_scaled, b, s);
  Trunc(b_scaled, Param::NBIT_K, Param::NBIT_K - Param::NBIT_F);

  Vec<ZZ_p> b_scaled_sq;
  MultElem(b_scaled_sq, b_scaled, b_scaled);
  Trunc(b_scaled_sq);

  Vec<ZZ_p> scaled_est;
  if (pid == 0) {
    scaled_est.SetLength(n);
  } else {
    scaled_est = - 10 * b_scaled + 5 * b_scaled_sq;
    if (pid == 1) {
      ZZ_p coeff;
      DoubleToFP(coeff, 5.9430, Param::NBIT_K, Param::NBIT_F);
      AddScalar(scaled_est, coeff);
    }
  }

  Vec<ZZ_p> w;
  MultElem(w, scaled_est, s);
  // scaled_est has bit length <= NBIT_F + 2, and s has bit length <= NBIT_K
  // so the bit length of w is at most NBIT_K + NBIT_F + 2
  Trunc(w, Param::NBIT_K + Param::NBIT_F + 2, Param::NBIT_K - Param::NBIT_F);

  Vec<ZZ_p> x;
  MultElem(x, w, b);
  Trunc(x);

  ZZ_p one;
  IntToFP(one, 1, Param::NBIT_K, Param::NBIT_F);

  x *= -1;
  if (pid == 1) {
    for (int i = 0; i < x.length(); i++) {
      x[i] += one;
    }
  }

  Vec<ZZ_p> y;
  MultElem(y, a, w);
  Trunc(y);

  for (int i = 0; i < niter; i++) {
    Vec<ZZ_p> xr, xm, yr, ym;
    BeaverPartition(xr, xm, x);
    BeaverPartition(yr, ym, y);

    Vec<ZZ_p> xpr = xr;
    if (pid > 0) {
      AddScalar(xpr, one);
    }

    Init(x, n);
    Init(y, n);

    BeaverMultElem(y, yr, ym, xpr, xm);
    BeaverMultElem(x, xr, xm, xr, xm);
    BeaverReconstruct(x);
    BeaverReconstruct(y);

    Trunc(x);
    Trunc(y);
  }

  if (pid == 1) {
    for (int i = 0; i < x.length(); i++) {
      x[i] += one;
    }
  }

  MultElem(c, y, x);
  Trunc(c);
}

void MPCEnv::Trunc(Mat<ZZ_p>& a, int k, int m) {
  if (debug) tcout() << "Trunc: " << a.NumRows() << ", " << a.NumCols() << endl; 

  Mat<ZZ_p> r;
  Mat<ZZ_p> r_low;
  if (pid == 0) {
    RandMatBits(r, a.NumRows(), a.NumCols(), k + Param::NBIT_V);

    r_low.SetDims(a.NumRows(), a.NumCols());
    for (int i = 0; i < a.NumRows(); i++) {
      for (int j = 0; j < a.NumCols(); j++) {
        r_low[i][j] = conv<ZZ_p>(trunc_ZZ(rep(r[i][j]), m));
      }
    }

    Mat<ZZ_p> r_mask;
    Mat<ZZ_p> r_low_mask;
    SwitchSeed(1);
    RandMat(r_mask, a.NumRows(), a.NumCols());
    RandMat(r_low_mask, a.NumRows(), a.NumCols());
    RestoreSeed();

    r -= r_mask;
    r_low -= r_low_mask;

    SendMat(r, 2);
    SendMat(r_low, 2);
  } else if (pid == 2) {
    ReceiveMat(r, 0, a.NumRows(), a.NumCols());
    ReceiveMat(r_low, 0, a.NumRows(), a.NumCols());
  } else {
    SwitchSeed(0);
    RandMat(r, a.NumRows(), a.NumCols());
    RandMat(r_low, a.NumRows(), a.NumCols());
    RestoreSeed();
  }

  Mat<ZZ_p> c;
  if (pid > 0) {
    c = a + r;
  } else {
    c.SetDims(a.NumRows(), a.NumCols());
  }

  RevealSym(c);
  
  Mat<ZZ_p> c_low;
  c_low.SetDims(a.NumRows(), a.NumCols());
  if (pid > 0) {
    for (int i = 0; i < a.NumRows(); i++) {
      for (int j = 0; j < a.NumCols(); j++) {
        c_low[i][j] = conv<ZZ_p>(trunc_ZZ(rep(c[i][j]), m));
      }
    }
  }

  if (pid > 0) {
    a += r_low;
    if (pid == 1) {
      a -= c_low;
    }

    ZZ_p twoinvm;
    map<int, ZZ_p>::iterator it = invpow_cache.find(m);
    if (it == invpow_cache.end()) {
      ZZ_p two(2);
      ZZ_p twoinv;
      inv(twoinv, two);
      power(twoinvm, twoinv, m);
      invpow_cache[m] = twoinvm;
    } else {
      twoinvm = it->second;
    }

    a *= twoinvm;
  }
}

void MPCEnv::PrefixOr(Mat<ZZ>& b, Mat<ZZ>& a, int fid) {
  if (debug) tcout() << "PrefixOr: " << a.NumRows() << ", " << a.NumCols() << endl;

  int n = a.NumRows();

  /* Find next largest squared integer */
  int L = (int) ceil(sqrt((double) a.NumCols()));
  int L2 = L * L;

  assert(primes[fid] > L + 1);

  /* Zero-pad to L2 bits */
  Mat<ZZ> a_padded;
  a_padded.SetDims(n, L2);

  if (pid > 0) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < L2; j++) {
        if (j < L2 - a.NumCols())
          a_padded[i][j] = 0;
        else
          a_padded[i][j] = a[i][j - L2 + a.NumCols()];
      }
    }
  }

  Reshape(a_padded, n * L, L);
  
  Vec<ZZ> x;
  FanInOr(x, a_padded, fid);

  Mat<ZZ> xpre;
  xpre.SetDims(n * L, L);
  
  if (pid > 0) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < L; j++) {
        int xpi = L * i + j;
        for (int k = 0; k < L; k++) {
          xpre[xpi][k] = (k <= j) ? x[L * i + k] : ZZ(0);
        }
      }
    }
  }

  Vec<ZZ> y;
  FanInOr(y, xpre, fid);
  xpre.kill();

  Vec< Mat<ZZ> > f; // f is a concatenation of n 1-by-L matrices
  f.SetLength(n);
  for (int i = 0; i < n; i++) {
    f[i].SetDims(1, L);
  }

  if (pid > 0) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < L; j++) {
        if (j == 0) {
          f[i][0][j] = x[L * i];
        } else {
          f[i][0][j] = y[L * i + j] - y[L * i + j - 1];
        }
      }
      Mod(f[i], fid);
    }
  }
  x.kill();

  Vec< Mat<ZZ> > tmp;
  tmp.SetLength(n);
  for (int i = 0; i < n; i++) {
    tmp[i].SetDims(L, L);
  }

  if (pid > 0) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < L; j++) {
        tmp[i][j] = a_padded[L * i + j];
      }
    }
  }
  a_padded.kill();

  Vec< Mat<ZZ> > c;
  MultMatParallel(c, f, tmp, fid); // c is a concatenation of n 1-by-L matrices
  tmp.kill();

  Mat<ZZ> cpre;
  cpre.SetDims(n * L, L);
  if (pid > 0) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < L; j++) {
        int cpi = L * i + j;
        for (int k = 0; k < L; k++) {
          cpre[cpi][k] = (k <= j) ? c[i][0][k] : ZZ(0);
        }
      }
    }
  }
  c.kill();

  Vec<ZZ> bdot_vec;
  FanInOr(bdot_vec, cpre, fid);
  cpre.kill();

  Vec< Mat<ZZ> > bdot;
  bdot.SetLength(n);
  for (int i = 0; i < n; i++) {
    bdot[i].SetDims(1, L);
  }

  if (pid > 0) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < L; j++) {
        bdot[i][0][j] = bdot_vec[L * i + j];
      }
    }
  }
  bdot_vec.kill();

  for (int i = 0; i < n; i++) {
    Reshape(f[i], L, 1);
  }

  Vec< Mat<ZZ> > s;
  MultMatParallel(s, f, bdot, fid);
  bdot.kill();

  b.SetDims(n, a.NumCols());
  if (pid > 0) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < a.NumCols(); j++) {
        int j_pad = L2 - a.NumCols() + j; 

        int il = (int) (j_pad / L);
        int jl = j_pad - il * L;

        b[i][j] = s[i][il][jl] + y[L * i + il] - f[i][il][0];
      }
    }
  }
  Mod(b, fid);
  s.kill();
  y.kill();
  f.kill();
}

void MPCEnv::FanInOr(Vec<ZZ>& b, Mat<ZZ>& a, int fid) {
  if (debug) tcout() << "FanInOr: " << a.NumRows() << ", " << a.NumCols() << endl;

  int n = a.NumRows();
  int d = a.NumCols();

  Vec<ZZ> a_sum;
  a_sum.SetLength(n);

  if (pid > 0) {
    for (int i = 0; i < n; i++) {
      a_sum[i] = (pid == 1) ? 1 : 0;
      for (int j = 0; j < d; j++) {
        a_sum[i] += a[i][j];
      }
    }
    Mod(a_sum, fid);
  }

  Mat<ZZ> coeff;
  coeff.SetDims(1, d + 1);
  pair<int, int> key = make_pair(d + 1, fid);
  if (or_lagrange_cache.find(key) != or_lagrange_cache.end()) {
    coeff[0] = or_lagrange_cache[key];
  } else {
    Vec<ZZ> y;
    y.SetLength(d + 1);
    for (int i = 0; i < d + 1; i++) {
      y[i] = (i == 0) ? 0 : 1;
    }
    lagrange_interp_simple(coeff[0], y, fid); // OR function
    or_lagrange_cache[key] = coeff[0];
  }

  Mat<ZZ> bmat;
  EvaluatePoly(bmat, a_sum, coeff, fid);
  b = bmat[0];
}

void MPCEnv::ShareRandomBits(Vec<ZZ_p>& r, Mat<ZZ>& rbits, int k, int n, int fid) {
  if (debug) tcout() << "ShareRandomBits: " << n << endl;

  if (pid == 0) {
    RandVecBits(r, n, k + Param::NBIT_V);
    NumToBits(rbits, r, k);

    Vec<ZZ_p> r_mask;
    Mat<ZZ> rbits_mask;

    SwitchSeed(1);
    RandVec(r_mask, n);
    RandMat(rbits_mask, n, k, fid);
    RestoreSeed();

    r -= r_mask;
    rbits -= rbits_mask;
    Mod(rbits, fid);

    SendVec(r, 2);
    SendMat(rbits, 2, fid);
  } else if (pid == 2) {
    ReceiveVec(r, 0, n);
    ReceiveMat(rbits, 0, n, k, fid);
  } else {
    SwitchSeed(0);
    RandVec(r, n);
    RandMat(rbits, n, k, fid);
    RestoreSeed();
  }
}

void MPCEnv::TableLookup(Mat<ZZ_p>& b, Vec<ZZ_p>& a, int table_id) {
  if (debug) tcout() << "TableLookup: " << a.length() << endl; 

  assert(!table_type_ZZ[table_id]);

  EvaluatePoly(b, a, lagrange_cache[table_id]);
}

void MPCEnv::TableLookup(Mat<ZZ_p>& b, Vec<ZZ>& a, int table_id, int fid) {
  if (debug) tcout() << "TableLookup: " << a.length() << endl; 

  assert(table_type_ZZ[table_id]);
  assert(table_field_index[table_id] == fid);

  int s = table_cache[table_id].NumCols();
  int n = a.length();

  Vec<ZZ_p> a_exp;
  a_exp.SetLength(n);
  if (pid > 0) {
    for (int i = 0; i < n; i++) {
      a_exp[i] = conv<ZZ_p>(a[i]);
    }
  }

  if (debug) tcout() << "Evaluating polynomial" << endl;
  if (debug) tcout() << s << ", " << lagrange_cache[table_id].NumCols() << endl;

  EvaluatePoly(b, a_exp, lagrange_cache[table_id]);
}

// Base field index 1
void MPCEnv::NormalizerEvenExp(Vec<ZZ_p>& b, Vec<ZZ_p>& b_sqrt, Vec<ZZ_p>& a) {
  if (debug) tcout() << "NormalizerEvenExp: " << a.length() << endl;

  int n = a.length();
  int fid = 1;

  Vec<ZZ_p> r; 
  Mat<ZZ> rbits;
  ShareRandomBits(r, rbits, Param::NBIT_K, n, fid);

  Vec<ZZ_p> e;
  if (pid == 0) {
    e.SetLength(n);
  } else {
    e = a + r;
  }
  r.kill();

  RevealSym(e);

  Mat<ZZ> ebits;
  if (pid == 0) {
    ebits.SetDims(n, Param::NBIT_K);
  } else {
    NumToBits(ebits, e, Param::NBIT_K);
  }
  e.kill();

  Vec<ZZ> c;
  LessThanBitsPublic(c, rbits, ebits, fid);
  if (pid > 0) {
    c = -c;
    if (pid == 1) {
      for (int i = 0; i < n; i++) {
        c[i] += 1;
      }
    }
    Mod(c, fid);
  }

  Mat<ZZ> ep;
  ep.SetDims(n, Param::NBIT_K + 1);
  if (pid > 0) {
    for (int i = 0; i < n; i++) {
      ep[i][0] = c[i];
      for (int j = 1; j < Param::NBIT_K + 1; j++) {
        ep[i][j] = (1 - 2 * ebits[i][j-1]) * rbits[i][j-1];
        if (pid == 1) {
          ep[i][j] += ebits[i][j-1];
        }
      }
    }
    Mod(ep, fid);
  }
  c.kill();

  Mat<ZZ> E;
  PrefixOr(E, ep, fid);
  ep.kill();

  Mat<ZZ> tpneg;
  tpneg.SetDims(n, Param::NBIT_K);
  if (pid > 0) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < Param::NBIT_K; j++) {
        tpneg[i][j] = E[i][j] - rbits[i][j] * (1 - ebits[i][j]);
      }
    }
    Mod(tpneg, fid);
  }
  E.kill();

  Mat<ZZ> Tneg;
  PrefixOr(Tneg, tpneg, fid);
  tpneg.kill();

  int half_len = Param::NBIT_K / 2;

  Mat<ZZ> efir, rfir;
  efir.SetDims(n, Param::NBIT_K);
  rfir.SetDims(n, Param::NBIT_K);
  if (pid > 0) {
    mul_elem(efir, ebits, Tneg);
    Mod(efir, fid);
  }
  MultElem(rfir, rbits, Tneg, fid);
  ebits.kill();
  rbits.kill();

  Vec<ZZ> double_flag;
  LessThanBits(double_flag, efir, rfir, fid);
  efir.kill();
  rfir.kill();
  
  Mat<ZZ> odd_bits, even_bits;
  odd_bits.SetDims(n, half_len);
  even_bits.SetDims(n, half_len);
  if (pid > 0) {
    for (int i = 0; i < n; i++) {
      for (int j = 0; j < half_len; j++) {
        odd_bits[i][j] = (pid == 1) ? (1 - Tneg[i][2*j+1]) : -Tneg[i][2*j+1];
        if ((2 * j + 2) < Param::NBIT_K) {
          even_bits[i][j] = (pid == 1) ? (1 - Tneg[i][2*j+2]) : -Tneg[i][2*j+2];
        } else {
          even_bits[i][j] = 0;
        }
      }
    }
    Mod(odd_bits, fid);
    Mod(even_bits, fid);
  }
  Tneg.kill();

  Vec<ZZ> odd_bit_sum, even_bit_sum;
  Init(odd_bit_sum, n);
  Init(even_bit_sum, n);
  for (int i = 0; i < n; i++) {
    for (int j = 0; j < half_len; j++) {
      odd_bit_sum[i] += odd_bits[i][j];
      even_bit_sum[i] += even_bits[i][j];
    }
    if (pid == 1) {
      odd_bit_sum[i] += 1;
      even_bit_sum[i] += 1;
    }
  }
  Mod(odd_bit_sum, fid);
  Mod(even_bit_sum, fid);
  odd_bits.kill();
  even_bits.kill();

  // If double_flag = true, then use odd_bits, otherwise use even_bits

  Vec<ZZ> diff;
  if (pid == 0) {
    diff.SetLength(n);
  } else {
    diff = odd_bit_sum - even_bit_sum;
    Mod(diff, fid);
  }
  MultElem(diff, double_flag, diff, fid);
  double_flag.kill();

  Vec<ZZ> chosen_bit_sum;
  if (pid == 0) {
    chosen_bit_sum.SetLength(n);
  } else {
    chosen_bit_sum = even_bit_sum + diff;
    Mod(chosen_bit_sum, fid);
  }
  odd_bit_sum.kill();
  even_bit_sum.kill();
  diff.kill();

  Mat<ZZ_p> b_mat;
  TableLookup(b_mat, chosen_bit_sum, 1, fid);

  if (pid > 0) {
    b_sqrt = b_mat[0];
    b = b_mat[1];
  } else {
    b_sqrt.SetLength(n);
    b.SetLength(n);
  }
}

void MPCEnv::ReadFromFile(ZZ_p& a, ifstream& ifs) {
  Vec<ZZ_p> avec;
  if (pid > 0) {
    Read(avec, ifs, 1);
    a = avec[0];
  }
}
void MPCEnv::ReadFromFile(Vec<ZZ_p>& a, ifstream& ifs, int n) {
  if (pid > 0) {
    Read(a, ifs, n);
  } else {
    a.SetLength(n);
  }
}
void MPCEnv::ReadFromFile(Mat<ZZ_p>& a, ifstream& ifs, int nrow, int ncol) {
  if (pid > 0) {
    Read(a, ifs, nrow, ncol);
  } else {
    a.SetDims(nrow, ncol);
  }
}
void MPCEnv::WriteToFile(ZZ_p& a, fstream& ofs) {
  if (pid > 0) {
    Vec<ZZ_p> avec;
    avec.SetLength(1);
    avec[0] = a;
    Write(avec, ofs);
  }
}
void MPCEnv::WriteToFile(Vec<ZZ_p>& a, fstream& ofs) {
  if (pid > 0) {
    Write(a, ofs);
  }
}
void MPCEnv::WriteToFile(Mat<ZZ_p>& a, fstream& ofs) {
  if (pid > 0) {
    Write(a, ofs);
  }
}

void MPCEnv::Write(Vec<ZZ_p>& a, fstream& ofs) {
  Mat<ZZ_p> a_copy;
  a_copy.SetDims(1, a.length());
  a_copy[0] = a;
  Write(a_copy, ofs);
}

void MPCEnv::Read(Vec<ZZ_p>& a, ifstream& ifs, int n) {
  Mat<ZZ_p> tmp;
  Read(tmp, ifs, 1, n);
  a = tmp[0];
}

void MPCEnv::ReadWithFilter(Vec<ZZ_p>& a, ifstream& ifs, Vec<ZZ_p>& filt) {
  assert(ifs.is_open());
  a.SetLength(filt.length());

  unsigned char *buf_ptr = buf;
  uint64_t stored_in_buf = 0;

  for (int i = 0; i < filt.length(); i++) {
    if (filt[i] != 1) {
      uint64_t count = 0;
      int k = i;
      while (k < filt.length() && filt[k] != 1) {
        k++;
        count++;
      }
      ifs.ignore(count * ZZ_bytes[0]);
      i += count - 1;
    } else {
      if (stored_in_buf == 0) {
        uint64_t count = 0;
        int k = i;
        while (k < filt.length() && filt[k] == 1) {
          k++;
          count++;
        }

        if (count > ZZ_per_buf[0]) {
          count = ZZ_per_buf[0];
        }

        ifs.read((char *)buf, count * ZZ_bytes[0]);
        stored_in_buf += count;
        buf_ptr = buf;
      }

      a[i] = conv<ZZ_p>(ZZFromBytes(buf_ptr, ZZ_bytes[0]));
      buf_ptr += ZZ_bytes[0];
      stored_in_buf--;
    }
  }
}

void MPCEnv::Write(Mat<ZZ_p>& a, fstream& ofs) {
  assert(ofs.is_open());

  unsigned char *buf_ptr = buf;
  uint64_t stored_in_buf = 0;
  for (int i = 0; i < a.NumRows(); i++) {
    for (int j = 0; j < a.NumCols(); j++) {
      if (stored_in_buf == ZZ_per_buf[0]) {
        ofs.write((const char *)buf, ZZ_bytes[0] * stored_in_buf);
        stored_in_buf = 0;
        buf_ptr = buf;
      }

      BytesFromZZ(buf_ptr, rep(a[i][j]), ZZ_bytes[0]);
      stored_in_buf++;
      buf_ptr += ZZ_bytes[0];
    }
  }

  if (stored_in_buf > 0) {
    ofs.write((const char *)buf, ZZ_bytes[0] * stored_in_buf);
  }
}

void MPCEnv::SkipData(ifstream& ifs, int n) {
  if (pid > 0) {
    assert(ifs.is_open());
    ifs.ignore(n * ZZ_bytes[0]);
  }
}

void MPCEnv::SkipData(ifstream& ifs, int nrows, int ncols) {
  if (pid > 0) {
    assert(ifs.is_open());
    for (int i = 0; i < nrows; i++) {
      ifs.ignore(ncols * ZZ_bytes[0]);
    }
  }
}

void MPCEnv::Read(Mat<ZZ_p>& a, ifstream& ifs, int nrows, int ncols) {
  assert(ifs.is_open());

  a.SetDims(nrows, ncols);
  unsigned char *buf_ptr = buf;
  uint64_t stored_in_buf = 0;
  uint64_t remaining = nrows * ncols;
  for (int i = 0; i < a.NumRows(); i++) {
    for (int j = 0; j < a.NumCols(); j++) {
      if (stored_in_buf == 0) {
        uint64_t count;
        if (remaining < ZZ_per_buf[0]) {
          count = remaining;
        } else {
          count = ZZ_per_buf[0];
        }
        ifs.read((char *)buf, count * ZZ_bytes[0]);
        stored_in_buf += count;
        remaining -= count;
        buf_ptr = buf;
      }

      a[i][j] = conv<ZZ_p>(ZZFromBytes(buf_ptr, ZZ_bytes[0]));
      buf_ptr += ZZ_bytes[0];
      stored_in_buf--;
    }
  }
}

void MPCEnv::SendInt(int num, int to_pid) {
  tcout() << "SendInt called: num(" << num << "), to_pid(" << to_pid << ")" << endl;
  *((int *)buf) = num;
  sockets.find(to_pid)->second.Send(buf, sizeof(int));
}

int MPCEnv::ReceiveInt(int from_pid) {
  tcout() << "ReceiveInt called: from_pid(" << from_pid << ")" << endl;
  sockets.find(from_pid)->second.Receive(buf, sizeof(int));
  return *((int *)buf);
}

void MPCEnv::SendBool(bool flag, int to_pid) {
  tcout() << "SendBool called: flag(" << flag << "), to_pid(" << to_pid << ")" << endl;
  *((bool *)buf) = flag;
  sockets.find(to_pid)->second.Send(buf, sizeof(bool));
}

bool MPCEnv::ReceiveBool(int from_pid) {
  tcout() << "ReceiveBool called: from_pid(" << from_pid << ")" << endl;
  sockets.find(from_pid)->second.Receive(buf, sizeof(bool));
  return *((bool *)buf);
}

void MPCEnv::SwitchSeed(int pid) {
  prg.find(cur_prg_pid)->second = GetCurrentRandomStream();
  SetSeed(prg.find(pid)->second);
  cur_prg_pid = pid;
}

void MPCEnv::ExportSeed(fstream& ofs, int pid) {
  assert(ofs.is_open());

  RandomStream rs = prg.find(pid)->second;
  rs.serialize(buf);

  ofs.write((const char *)buf, RandomStream::numBytes());
}

void MPCEnv::ExportSeed(fstream& ofs) {
  assert(ofs.is_open());

  RandomStream rs = GetCurrentRandomStream();
  rs.serialize(buf);

  ofs.write((const char *)buf, RandomStream::numBytes());
}

void MPCEnv::ImportSeed(int newid, ifstream& ifs) {
  assert(ifs.is_open());

  ifs.read((char *)buf, RandomStream::numBytes());

  RandomStream rs((const unsigned char *)buf, true);

  pair<map<int,RandomStream>::iterator,bool> ret;
  ret = prg.insert(pair<int, RandomStream>(newid, rs));
  if (!ret.second) { // ID exists already
    ret.first->second = rs;
  }
}

void MPCEnv::BeaverReadFromFile(Mat<ZZ_p>& ar, Mat<ZZ_p>& am, ifstream& ifs, int nrow, int ncol) {
  if (pid > 0) {
    Read(ar, ifs, nrow, ncol);
  } else {
    ar.SetDims(nrow, ncol);
  }
  Read(am, ifs, nrow, ncol);
}

void MPCEnv::BeaverReadFromFile(Vec<ZZ_p>& ar, Vec<ZZ_p>& am, ifstream& ifs, int n) {
  if (pid > 0) {
    Read(ar, ifs, n);
  } else {
    ar.SetLength(n);
  }
  Read(am, ifs, n);
}

void MPCEnv::BeaverReadFromFileWithFilter(Vec<ZZ_p>& ar, Vec<ZZ_p>& am, ifstream& ifs, Vec<ZZ_p>& filt) {
  if (pid > 0) {
    ReadWithFilter(ar, ifs, filt);
  } else {
    ar.SetLength(filt.length());
  }
  ReadWithFilter(am, ifs, filt);
}

void MPCEnv::BeaverWriteToFile(Vec<ZZ_p>& ar, Vec<ZZ_p>& am, fstream& ofs) {
  if (pid > 0) {
    Write(ar, ofs);
  }
  Write(am, ofs);
}

void MPCEnv::BeaverWriteToFile(Mat<ZZ_p>& ar, Mat<ZZ_p>& am, fstream& ofs) {
  if (pid > 0) {
    Write(ar, ofs);
  }
  Write(am, ofs);
}

void MPCEnv::BeaverMultElem(Vec<ZZ_p>& ab, Vec<ZZ_p>& ar, Vec<ZZ_p>& am, Vec<ZZ_p>& br, Vec<ZZ_p>& bm, int fid) {
  if (pid == 0) {
    Vec<ZZ_p> ambm;
    mul_elem(ambm, am, bm);
    ab += ambm;
  } else {

    ZZ_pContext context;
    context.save();

    NTL_GEXEC_RANGE(ab.length() > Param::PAR_THRES, ab.length(), first, last)

    context.restore();

    for (int i = first; i < last; i++) {
      ab[i] += ar[i] * bm[i];
      ab[i] += am[i] * br[i];
      if (pid == 1) {
        ab[i] += ar[i] * br[i];
      }
    }

    NTL_GEXEC_RANGE_END
  }
}

void MPCEnv::BeaverMult(Mat<ZZ_p>& ab, Mat<ZZ_p>& ar, Mat<ZZ_p>& am, Mat<ZZ_p>& br, Mat<ZZ_p>& bm, bool elem_wise, int fid) {
  if (pid == 0) {
    Mat<ZZ_p> ambm;
    if (elem_wise) {
      mul_elem(ambm, am, bm);
    } else {
      mul(ambm, am, bm);
    }
    ab += ambm;
  } else {
    if (elem_wise) {

      ZZ_pContext context;
      context.save();

      NTL_GEXEC_RANGE(ab.NumRows() > Param::PAR_THRES, ab.NumRows(), first, last)

      context.restore();

      for (int i = first; i < last; i++) {
        for (int j = 0; j < ab.NumCols(); j++) {
          ab[i][j] += ar[i][j] * bm[i][j];
          ab[i][j] += am[i][j] * br[i][j];
          if (pid == 1) {
            ab[i][j] += ar[i][j] * br[i][j];
          }
        }
      }

      NTL_GEXEC_RANGE_END

    } else {
      ab += ar * bm;
      ab += am * br;
      if (pid == 1) {
        ab += ar * br;
      }
    }
  }
}

void MPCEnv::BeaverMultElem(Vec<ZZ>& ab, Vec<ZZ>& ar, Vec<ZZ>& am, Vec<ZZ>& br, Vec<ZZ>& bm, int fid) {
  if (pid == 0) {
    Vec<ZZ> ambm;
    mul_elem(ambm, am, bm);
    ab += ambm;
  } else {
    NTL_GEXEC_RANGE(ab.length() > Param::PAR_THRES, ab.length(), first, last)

    for (int i = first; i < last; i++) {
      ab[i] += ar[i] * bm[i];
      ab[i] += am[i] * br[i];
      if (pid == 1) {
        ab[i] += ar[i] * br[i];
      }
    }

    NTL_GEXEC_RANGE_END
  }

  Mod(ab, fid);
}

void MPCEnv::BeaverMult(Mat<ZZ>& ab, Mat<ZZ>& ar, Mat<ZZ>& am, Mat<ZZ>& br, Mat<ZZ>& bm, bool elem_wise, int fid) {
  if (pid == 0) {
    Mat<ZZ> ambm;
    if (elem_wise) {
      mul_elem(ambm, am, bm);
    } else {
      mul(ambm, am, bm);
    }
    ab += ambm;
  } else {
    if (elem_wise) {
      NTL_GEXEC_RANGE(ab.NumRows() > Param::PAR_THRES, ab.NumRows(), first, last)

      for (int i = first; i < last; i++) {
        for (int j = 0; j < ab.NumCols(); j++) {
          ab[i][j] += ar[i][j] * bm[i][j];
          ab[i][j] += am[i][j] * br[i][j];
          if (pid == 1) {
            ab[i][j] += ar[i][j] * br[i][j];
          }
        }
      }

      NTL_GEXEC_RANGE_END

    } else {
      ab += ar * bm;
      ab += am * br;
      if (pid == 1) {
        ab += ar * br;
      }
    }
  }

  Mod(ab, fid);
}

void MPCEnv::Quicksort(Mat<ZZ_p>& a, int column) {
  Quicksort(a, column, 0, -1, -1);
}

void MPCEnv::Quicksort(Mat<ZZ_p>& a, int column,
                       int depth, int start, int end) {
  if (start == -1) {
    start = 0;
  }
  if (end == -1) {
    end = a.NumRows() - 1;
  }
  for (int i = 0; i < depth; i++) {
    cout << "\t";
  }
  cout << start << " " << end << endl;

  if (end - start + 1 > 1) {
    int p;
    QuicksortPartition(p, a, column, start, end);
    Quicksort(a, column, depth + 1, start, p - 1);
    Quicksort(a, column, depth + 1, p + 1, end);
  }
}

void MPCEnv::QuicksortPartition(int& p, Mat<ZZ_p>& a, int column,
                                int start, int end) {
  int i = start;
  for (int j = start; j < end - 1; j++) {
    Vec<ZZ_p> c, a_j, a_end;
    c.SetLength(1);
    a_j.SetLength(1);
    a_j[0] = a[j][column];
    a_end.SetLength(1);
    a_end[0] = a[end][column];
    
    LessThan(c, a_j, a_end);
    RevealSym(c);
    double lt = FPToDouble(c[0], Param::NBIT_K, Param::NBIT_F);
    if (lt) {
      swap(a[i], a[j]);
      i++;
    }
  }
  p = i;
  swap(a[p], a[end]);
}
