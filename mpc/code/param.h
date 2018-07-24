#ifndef __PARAM_H__
#define __PARAM_H__

#include <iostream>
using namespace std;

class Param {
  public:
    /* IP addresses */
    static string IP_ADDR_P0;
    static string IP_ADDR_P1;
    static string IP_ADDR_P2;

    /* Ports */
    // The party with smaller ID listens on the port 
    // and the other connects to it. Make sure the firewall 
    // setting of the listener allows these ports.
    static int PORT_P0_P1;
    static int PORT_P0_P2;
    static int PORT_P1_P2;
    static int PORT_P1_P3;
    static int PORT_P2_P3;
    
    /* Directories and files */
    static string KEY_PATH; // encryption keys for secure channels
    static string LOG_FILE; // runtime/communication profiler output file
    static string OUTPUT_FILE_PREFIX; // prefix for output
    static string CACHE_FILE_PREFIX; // prefix for cache files; includes
                                     // input shares, which can be large

    /* Secret sharing parameters (see Supplementary Information) */
    // Require:
    //   NBIT_F < NBIT_K
    //   Both NBIT_K and NBIT_K - NBIT_F are even
    //   log2(BASE_P) > NBIT_K + NBIT_F + NBIT_V + 2
    // Useful resource:
    //   https://primes.utm.edu/lists/2small/200bit.html
    static int NBIT_K; // total bit length of a data value
    static int NBIT_F; // bit length of fractional range
    static int NBIT_V; // security parameter
    static string BASE_P; // base prime for finite field
    
    /* Algorithm parameters */
    static long N_INTERACTIONS;
    static long N_DRUGS;
    static long N_TARGETS;
    static long FEATURE_RANK;
    static int MAX_EPOCHS;
    static double LAMBDA_D;
    static double LAMBDA_T;
    static double MOMENTUM;
    static double LEARN_RATE;
    static int N_CLASSES;
    static int N_HIDDEN;
    static int N_NEURONS;
    static bool ANNEAL;
    static int ANNEAL_FREQ;
    static double REG;
    static double DROPOUT;
    static string LOSS;
    static int BATCH_SIZE;
    static int N_FILE_BATCH;
    static string FEATURES_FILE;
    static string LABELS_FILE;
    static int ITER_PER_EVAL;
    
    /* Software parameters */
    static uint64_t MPC_BUF_SIZE; // size of data buffer in bytes
    static long DIV_MAX_N; // maximum number of parallel divisions/sqrts
    static long PITER_BATCH_SIZE; // size of batch processing of genomes
    static long PAR_THRES; // minimum threshold for thread-boosting
    static long NUM_THREADS; // number of threads for thread-boosting

    /* Global control */
    static bool PROFILER; // turn on/off profiler for runtime/communication
    static bool DEBUG; // debug flag
    
    /* Helper functions */
    template<class T>
    static bool Convert(string s, T &var, string name);
    static bool ParseFile(const char *param_file);
};

#endif 
