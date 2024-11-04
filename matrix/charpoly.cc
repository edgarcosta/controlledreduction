// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "matrix.h"

#include <assert.h>
#include <stdint.h>

#include <NTL/matrix.h>
#include <NTL/ZZ.h>

using namespace NTL;

Vec<ZZ> charpoly(const Mat<ZZ> M) {
    int64_t n = M.NumRows();
    assert( n == (int64_t) M.NumCols() );

    int64_t rem;
    Vec<ZZ> result, e, p;
    Mat<ZZ> powerM;
    ZZ sump(0);
    ZZ sumn(0);

    result.SetLength(n+1);
    e.SetLength(n+1);
    p.SetLength(n+1);
    NTL::set(e[0]);
    NTL::set(result[n]);
    p[1] = trace(M);

    powerM = M;

    for( int64_t i = 2; i <=n; ++i) { 
        powerM *= M; //powerM = M^i;
        p[i] = trace(powerM);
    }

    for( int64_t k = 1; k <= n; k++) {
        clear(sumn);
        clear(sump);

        for( int64_t i = 1; i <= k ; i+=2)
            sump += e[k-i]*p[i];

        for( int64_t i = 2; i <= k; i+=2)
            sumn += e[k-i]*p[i];

        rem = DivRem(e[k], sump - sumn, k);
        assert(rem == 0);
        result[n-k] = (k%2 == 0)? e[k]: -e[k];
    }
    return result;
}

