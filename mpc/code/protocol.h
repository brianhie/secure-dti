#ifndef __PROTOCOL_H_
#define __PROTOCOL_H_

#include "mpc.h"
#include "util.h"
#include <vector>
#include <NTL/mat_ZZ_p.h>
#include <NTL/mat_ZZ.h>
#include <NTL/ZZ.h>
#include <NTL/BasicThreadPool.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include <chrono>

using namespace NTL;
using namespace std;

using msec = chrono::milliseconds;
using get_time = chrono::steady_clock;

#define ABS(a) (((a)<0)?-(a):(a))

auto clock_start = get_time::now();

void tic() {
  clock_start = get_time::now();
}

int toc() {
  auto clock_end = get_time::now();
  int duration = chrono::duration_cast<msec>(clock_end - clock_start).count();
  tcout() << "Elapsed time is " << duration / 1000.0 << " secs" << endl;
  return duration;
}

string cache(int pid, string desc) {
  ostringstream oss;
  oss << Param::CACHE_FILE_PREFIX << "_" << desc << ".bin";
  return oss.str();
}

string cache(int pid, int index) {
  ostringstream oss;
  oss << Param::CACHE_FILE_PREFIX << "_" << index << ".bin";
  return oss.str();
}

string outname(string desc) {
  ostringstream oss;
  oss << Param::OUTPUT_FILE_PREFIX << "_" << desc << ".txt";
  return oss.str();
}

bool measure_network_speed_protocol(MPCEnv& mpc, int pid) {
  Vec<ZZ_p> a;
  Init(a, 10000000);
  uint64_t tot = 0;
  int niter = 10;
  fstream fs(cache(2, "test"), ios::out | ios::binary);
  for (int i = 0; i < niter; i++) {
    tic();
    /*
    if (pid == 1) {
      mpc.SendVec(a, 2);
    } else if (pid == 2) {
      mpc.ReceiveVec(a, 1, a.length());
    }
    */
    mpc.WriteToFile(a, fs);
    int d = toc();
    tot += d;
  }

  double nbytes = a.length() * NumBytes(ZZ_p::modulus());
  double secs = tot / ((double)niter) / 1000;
  tcout() << nbytes << endl;
  tcout() << secs << endl;
  tcout() << nbytes / secs << endl;

  return true;
}

bool div_test(MPCEnv& mpc, int pid) {
  Vec<ZZ_p> a, b, c;
  Init(a, 543653);
  ZZ_p fp_one = DoubleToFP(1, Param::NBIT_K, Param::NBIT_F);
  AddScalar(a, fp_one);
  tic();
  b = 2 * a;
  mpc.FPDiv(c, a, b);
  toc();
  return true;
}

bool unit_test(MPCEnv& mpc, int pid) {
  ZZ_p x, y, z;
  Vec<ZZ_p> xv, yv, zv, wv;
  Vec<double> xdv, ydv, zdv, wdv;
  Mat<ZZ_p> xm, ym, zm;
  ZZ a, b, c;
  Vec<ZZ> av, bv, cv, dv;
  Mat<ZZ> am, bm, cm, dm;
  double d;
  double eps = 1e-6;

  tcout() << "[Fixed-point ZZ_p <-> Double conversion] ";
  x = DoubleToFP(3.141592653589793238462643383279, Param::NBIT_K, Param::NBIT_F);
  d = ABS(FPToDouble(x, Param::NBIT_K, Param::NBIT_F) - 3.141592653589793238462643383279);
  if (pid > 0) {
    assert(d < eps);
    tcout() << "Success";
  }
  tcout() << endl;

  tcout() << "[FP multiplcation] ";
  Init(xv, 3); Init(yv, 3);
  if (pid == 2) {
    xv[0] = DoubleToFP(1.34, Param::NBIT_K, Param::NBIT_F);
    xv[1] = DoubleToFP(100.3, Param::NBIT_K, Param::NBIT_F);
    xv[2] = DoubleToFP(-0.304, Param::NBIT_K, Param::NBIT_F);
    yv[0] = DoubleToFP(-0.001, Param::NBIT_K, Param::NBIT_F);
    yv[1] = DoubleToFP(303, Param::NBIT_K, Param::NBIT_F);
    yv[2] = DoubleToFP(-539, Param::NBIT_K, Param::NBIT_F);
  }
  mpc.MultElem(zv, xv, yv);
  mpc.Trunc(zv);
  mpc.RevealSym(zv);

  FPToDouble(zdv, zv, Param::NBIT_K, Param::NBIT_F);
  if (pid > 0) {
    tcout() << zdv << endl;
    assert(ABS(zdv[0] - (-0.00134)) < eps);
    assert(ABS(zdv[1] - (30390.9)) < eps);
    assert(ABS(zdv[2] - (163.856)) < eps);
    tcout() << "Success";
  }
  tcout() << endl;

  tcout() << "[Powers]" << endl;;
  Init(xv, 5);
  if (pid == 1) {
    xv[0] = 0;
    xv[1] = 1;
    xv[2] = 2;
    xv[3] = 3;
    xv[4] = 4;
  }

  mpc.Powers(ym, xv, 3);
  mpc.Print(ym, cout);

  tcout() << "[FanInOr]" << endl;
  Init(am, 5, 3);
  if (pid == 1) {
    am[0][0] = 0; am[0][1] = 1; am[0][2] = 1;
    am[1][0] = 0; am[1][1] = 1; am[1][2] = 0;
    am[2][0] = 0; am[2][1] = 0; am[2][2] = 0;
    am[3][0] = 0; am[3][1] = 0; am[3][2] = 1;
    am[4][0] = 1; am[4][1] = 1; am[4][2] = 1;
  }
  mpc.FanInOr(bv, am, 2);
  mpc.Print(bv, 100, 2);

  tcout() << "[LessThanBitsPublic]" << endl;
  Init(bm, 5, 3);
    bm[0][0] = 0; bm[0][1] = 0; bm[0][2] = 1;
    bm[1][0] = 0; bm[1][1] = 1; bm[1][2] = 1;
    bm[2][0] = 0; bm[2][1] = 0; bm[2][2] = 1;
    bm[3][0] = 1; bm[3][1] = 0; bm[3][2] = 1;
    bm[4][0] = 0; bm[4][1] = 1; bm[4][2] = 1;
  mpc.LessThanBitsPublic(cv, am, bm, 2);
  mpc.Print(cv, 100, 2);

  tcout() << "[TableLookup]" << endl;
  Init(av, 5);
  if (pid == 1) {
    for (int i = 0; i < 5; i++) {
      av[i] = i+1;
    }
  }
  mpc.TableLookup(ym, av, 1, 1);
  mpc.Print(ym);

  tcout() << "[FP normalization] ";
  yv[0] *= -1;
  yv[2] *= -1;
  mpc.NormalizerEvenExp(zv, wv, yv);

  mpc.MultElem(zv, yv, zv);
  mpc.Trunc(zv, Param::NBIT_K, Param::NBIT_K - Param::NBIT_F);
  mpc.RevealSym(zv);
  FPToDouble(zdv, zv, Param::NBIT_K, Param::NBIT_F);
  tcout() << zdv;
  tcout() << endl;

  tcout() << "[FP sqrt] ";
  Init(xv, 3);
  if (pid == 2) {
    xv[0] = DoubleToFP(0.001, Param::NBIT_K, Param::NBIT_F);
    xv[1] = DoubleToFP(303, Param::NBIT_K, Param::NBIT_F);
    xv[2] = DoubleToFP(539, Param::NBIT_K, Param::NBIT_F);
  }
  mpc.FPSqrt(yv, zv, xv);
  mpc.PrintFP(yv);
  mpc.PrintFP(zv);
  tcout() << endl;

  tcout() << "[FP division] ";
  Init(xv, 3); Init(yv, 3);
  if (pid == 2) {
    xv[0] = DoubleToFP(1.34, Param::NBIT_K, Param::NBIT_F);
    xv[1] = DoubleToFP(100.3, Param::NBIT_K, Param::NBIT_F);
    xv[2] = DoubleToFP(-0.304, Param::NBIT_K, Param::NBIT_F);
    yv[0] = DoubleToFP(0.001, Param::NBIT_K, Param::NBIT_F);
    yv[1] = DoubleToFP(303, Param::NBIT_K, Param::NBIT_F);
    yv[2] = DoubleToFP(539, Param::NBIT_K, Param::NBIT_F);
  }
  mpc.FPDiv(zv, xv, yv);
  mpc.RevealSym(zv);

  FPToDouble(zdv, zv, Param::NBIT_K, Param::NBIT_F);
  if (pid > 0) {
    tcout() << zdv << endl;
    assert(ABS(zdv[0] - (1340.000000000000000)) < eps);
    assert(ABS(zdv[1] - (0.331023102310231)) < eps);
    assert(ABS(zdv[2] - (-0.000564007421150)) < eps);
    tcout() << "Success";
  }
  tcout() << endl;

  tcout() << "[Householder] ";
  mpc.Householder(yv, xv);
  mpc.PrintFP(xv);
  mpc.PrintFP(yv);

  tcout() << "[Eigendecomp]";
  Init(xm, 5, 5);
  if (pid == 2) {
    xm[0][0] = DoubleToFP(1.34, Param::NBIT_K, Param::NBIT_F);
    xm[0][1] = DoubleToFP(0, Param::NBIT_K, Param::NBIT_F);
    xm[0][2] = DoubleToFP(-3, Param::NBIT_K, Param::NBIT_F);
    xm[0][3] = DoubleToFP(5, Param::NBIT_K, Param::NBIT_F);
    xm[0][4] = DoubleToFP(0.003, Param::NBIT_K, Param::NBIT_F);
    xm[1][0] = DoubleToFP(10, Param::NBIT_K, Param::NBIT_F);
    xm[1][1] = DoubleToFP(1, Param::NBIT_K, Param::NBIT_F);
    xm[1][2] = DoubleToFP(2.2, Param::NBIT_K, Param::NBIT_F);
    xm[1][3] = DoubleToFP(3.33, Param::NBIT_K, Param::NBIT_F);
    xm[1][4] = DoubleToFP(4.444, Param::NBIT_K, Param::NBIT_F);
    xm[2][0] = DoubleToFP(5, Param::NBIT_K, Param::NBIT_F);
    xm[2][1] = DoubleToFP(4, Param::NBIT_K, Param::NBIT_F);
    xm[2][2] = DoubleToFP(3, Param::NBIT_K, Param::NBIT_F);
    xm[2][3] = DoubleToFP(2, Param::NBIT_K, Param::NBIT_F);
    xm[2][4] = DoubleToFP(1, Param::NBIT_K, Param::NBIT_F);
    xm[3][0] = DoubleToFP(0, Param::NBIT_K, Param::NBIT_F);
    xm[3][1] = DoubleToFP(2, Param::NBIT_K, Param::NBIT_F);
    xm[3][2] = DoubleToFP(5, Param::NBIT_K, Param::NBIT_F);
    xm[3][3] = DoubleToFP(1, Param::NBIT_K, Param::NBIT_F);
    xm[3][4] = DoubleToFP(0, Param::NBIT_K, Param::NBIT_F);
    xm[4][0] = DoubleToFP(1, Param::NBIT_K, Param::NBIT_F);
    xm[4][1] = DoubleToFP(2, Param::NBIT_K, Param::NBIT_F);
    xm[4][2] = DoubleToFP(1, Param::NBIT_K, Param::NBIT_F);
    xm[4][3] = DoubleToFP(2, Param::NBIT_K, Param::NBIT_F);
    xm[4][4] = DoubleToFP(1, Param::NBIT_K, Param::NBIT_F);
  }
  mpc.EigenDecomp(ym, yv, xm);
  mpc.PrintFP(yv);
  mpc.PrintFP(ym);
  tcout() << endl;
  //
  // This is here just to keep P0 online until the end for data transfer
  // In practice, P0 would send data in advance before each phase and go offline
  if (pid == 0) {
    mpc.ReceiveBool(2);
  } else if (pid == 2) {
    mpc.SendBool(true, 0);
  }

  mpc.CleanUp();

  return true;
}

#endif
