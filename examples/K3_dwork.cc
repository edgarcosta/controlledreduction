// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface_nd.h"
#include "timing.h"
#include "tools.h"
#include <time.h>
#include <cstdlib>

#include <sstream>

using namespace std;
using namespace NTL;

int usage(char name[]){
    cout << "Computes the zeta function of: x^4 + y^4 + z^4 + w^4 + lambda * x*y*z*w == 0" << endl;
    cout << "Usage: " << name << " <p> <lambda>" << endl;
    return 1;
}



int main(int argc, char* argv[]) {
    if(argc != 3) 
        return usage(argv[0]);
    int64_t p = atoi(argv[1]);
    if( not ProbPrime(ZZ(p)) ) {
        cout << p << " is not a prime!" << endl << endl; 
        usage(argv[0]);
    }
    int64_t lambda = atoi(argv[2]);
    int64_t n = 3;
    int64_t d = 4;


    Vec<int64_t> N;
    Vec<int64_t> charpoly_prec;
    int64_t precision;
    default_args(precision, N, charpoly_prec, p, n, d);
    {
        zz_pPush push(p);
        ZZ_pPush push2(power_ZZ(p,precision));

        map< Vec<int64_t>, zz_p , vi64less> f;

        // a very lazy way to create f
        for(int64_t i = 0; i < 4; ++i) {
            Vec<int64_t> v;
            v.SetLength(4, 0);
            v[i] = d;
            f[v] = 1;
        }
        Vec<int64_t> v;
        v.SetLength(4, 1);
        f[v] = lambda;

        cout <<= f;
        cout << endl;

        cout << p << endl;

        if( not isSmooth(f) )
        {
            cout << "The surface x^4 + y^4 + z^4 + w^4 + " << lambda <<" * x*y*z*w is not smooth!" << endl << endl; 
            usage(argv[0]);
        }



        hypersurface_non_degenerate hs(p, precision, f, false);
        Mat<ZZ_p> Frob = hs.frob_matrix_ND(N);
        Mat<ZZ>  Frob_ZZ = conv<Mat<ZZ> >(Frob);
        Vec<ZZ> cp = charpoly_frob(Frob_ZZ, charpoly_prec, p, n - 1);
        cout << cp << endl;
    }

    return 0;
}

