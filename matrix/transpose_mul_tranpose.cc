// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "matrix.h"


#include <NTL/matrix.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>

#include <cstdint>

using namespace std;


template <>
void transpose_mul_transpose<ZZ_p>(Mat<ZZ_p> &result,const Mat<ZZ_p> &A, const Mat<ZZ_p> &B)
{
    int64_t rowsA = A.NumRows();
    int64_t colsA = A.NumCols();
    int64_t colsB = B.NumRows();
    int64_t rowsB = B.NumCols();
    assert(colsA == rowsB);
    result.SetDims(colsB,rowsA);
    int64_t i, j, k;
    ZZ acc, tmp;
    for(j = 1; j <= colsB; j++)
    {
        for(i = 1; i <= rowsA; ++i )
        {
            clear(acc);
            for(k = 1; k <= rowsB; k++)
            {
                mul(tmp,rep(A(i,k)),rep(B(j,k)));
                add(acc,acc,tmp);
            }
            result(j,i) = conv<ZZ_p>(acc);
        }
    }
}

