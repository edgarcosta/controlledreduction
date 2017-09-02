// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
//
// the zeta function is the same under a change of variables

#include "hypersurface.h"
#include "timing.h"
#include "tools.h"

#define NTL_RANGE_CHECK
using namespace NTL;

void test_invariance_J(int64_t p, int64_t n, int64_t d)
{
    Vec<int64_t> N;
    Vec<int64_t> charpoly_prec;
    int64_t precision;
    default_args(precision, N, charpoly_prec, p, n, d);
    {
        zz_pPush push(p);
        ZZ_pPush push2(power_ZZ(p,precision));

        hypersurface hs;
        hs = hypersurface(p, precision, n, d, false);
        Mat<ZZ_p> Fp_J = hs.frob_matrix_J(N);
        Mat<ZZ>  Frob_ZZ = conv<Mat<ZZ> >(Fp_J);
        Vec<ZZ> cp = charpoly_frob(Frob_ZZ, charpoly_prec, p, n - 1);

        map< Vec<int64_t>, zz_p , vi64less> new_fbar;
        Mat<zz_p> slmatrix = random_SL_matrix(n+1);
        new_fbar =  change_of_variables((hs.dR)->fbar, slmatrix);

        hypersurface hs2;
        hs2 = hypersurface(p, precision, new_fbar, false);
        Mat<ZZ_p> Fp_J2 = hs.frob_matrix_J(N);
        Mat<ZZ>  Frob_ZZ2 = conv<Mat<ZZ> >(Fp_J2);
        Vec<ZZ> cp2 = charpoly_frob(Frob_ZZ, charpoly_prec, p, n - 1);
        if(cp2 != cp)
        {
            cout << "\n\tFAIL test_invariance_J" <<endl;
            cout << "\tp = " << p << " n = " << n << " d = " << d << " N = " << N << endl;
            cout << "\tslmatrix = " << slmatrix << endl;
            cout << "\tfbar = ";
            cout <<= (hs.dR)->fbar;
            cout << endl;
            cout << "\tnew_fbar = ";
            cout <<= new_fbar;
            cout << endl;
            abort();
        }
    }
}



int main()
{
    timestamp_pair pair;
    timestamp_mark(pair);
    SetSeed(to_ZZ(time(NULL)));
    cout<<"test_invariance_J...";
    fflush(stdout);
    int64_t n = 2;
    for (int64_t d = n + 1; d <= 5; ++d) {
        int64_t plen = 5;
        int64_t runs = 20;
        if (d == 5 ) {
            plen = 4;
            runs = 1;
        }
        
        for (int64_t i = 0; i < runs; ++i) {
            int64_t p = GenPrime_long(plen);
            while ( p <= d )
                p = GenPrime_long(plen);
            test_invariance_J(p, n, d);
        }
    }
/*
    it takes too long to be a regular test
    n = 3;
    int64_t d = 4;
    int64_t p = GenPrime_long(5);
    while ( p <= 23 )
        p = GenPrime_long(5);
    test_invariance_J(p, n, d);
*/
    cout << "PASS ";
    timestamp_report(pair);
    cout << endl;
    return 0;
}





