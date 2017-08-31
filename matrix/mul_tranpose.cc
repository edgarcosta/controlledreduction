// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "matrix.h"

#include <assert.h>
#include <gmp.h>
#include <NTL/matrix.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/ZZ.h>

#include <cstdint>

using namespace NTL;

// specialization of mul_transpose and transpose_mul_transpose
template <>
void mul_transpose<ZZ_p>(Mat<ZZ_p> &result,const Mat<ZZ_p> &A, const Mat<ZZ_p> &B)
{
    int64_t rowsA = A.NumRows();
    int64_t colsA = A.NumCols();
    int64_t colsB = B.NumRows();
    int64_t rowsB = B.NumCols();
    assert(colsA == rowsB);
    result.SetDims(rowsA,colsB);
    int64_t i, j, k;
    ZZ acc, tmp;
    for(i = 1; i <= rowsA; ++i)
    {
        for(j = 1; j <= colsB; j++)
        {
            clear(acc);
            for(k = 1; k <= rowsB; k++)
            {
                mul(tmp,rep(A(i,k)),rep(B(j,k)));
                add(acc,acc,tmp);
            }
            result(i,j) = conv<ZZ_p>(acc);
        }
    }
}


// result = A * B where B is given in column major, ie, as the transpose
void mul_transpose(mpz_t* result, mpz_t* A, int64_t rowsA, int64_t colsA, mpz_t* B, int64_t rowsB, int64_t colsB)
{
    assert(colsA == rowsB);
    int64_t i,j,k;
    mpz_t tmp, acc;
    mpz_inits(tmp, acc,NULL);
    for(i = 0; i < rowsA; ++i)
    {
        for(j = 0; j < colsB; j++)
        {
            mpz_set_ui(acc,0);
            for(k = 0; k < rowsB; k++)
            {
                mpz_mul(tmp, A[i*rowsA + k], B[j*colsB + k]);
                mpz_add(acc,acc,tmp);
            }
            mpz_set(result[i*rowsA + j],acc);
        }
    }
    mpz_clears(tmp, acc,NULL);
}

