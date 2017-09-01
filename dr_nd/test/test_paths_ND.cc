// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
//
// runs test_paths_ND for some random curves of degree 3,4 and 5 and K3 surfaces

   
#include "dr_nd.h"
#include "timing.h"

#include <iostream>

#define NTL_RANGE_CHECK
int main() {
    timestamp_pair pair;
    timestamp_mark(pair);
    SetSeed(to_ZZ(time(NULL)));
    cout<<"test_paths_ND...";
    fflush(stdout);
    for (int64_t i = 0; i < 10; ++i) {
        int64_t p = 0;
        while(p < 43)
            p = GenPrime_long(13);
        int64_t precision = 4;
        for (int64_t n = 2; n < 4; ++n) {
            int64_t degree_bound = 5;
            int64_t N = 5;
            if (n == 3) {
                degree_bound = 4;
                N = 3;
            }
            for (int64_t d = n + 1; d <= degree_bound; d++) {
                zz_p::init(p);
                ZZ_p::init( power_ZZ(p, precision) );
                de_Rham_non_degenerate_local D;
                //cout << p <<" "<<n<<" "<<d<<endl; 
                D = de_Rham_non_degenerate_local(p, precision, n, d, false);
                bool test_success = D.test_paths_ND(N, N);;
                if (not test_success ) {
                    cout << "\n\tFAIL test_paths_ND" <<endl;
                    cout << "\tp = " << p << " n = " << n << " d = " << d << " N = " << N << endl;
                    abort();
                }
            }
        }
    }
    cout << "PASS ";
    timestamp_report(pair);
    cout << endl;
    return 0;
}
