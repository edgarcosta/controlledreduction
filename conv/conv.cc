// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include <gmp.h>
#include <flint/fmpz.h>
#include <NTL/matrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/ZZ.h>

#include <cstdint>

using namespace NTL;

//conversion between ZZ and mpz

void conv(mpz_t a,const ZZ &b)
{
    long n = NumBytes(b);
    unsigned char* bytes = new unsigned char[n];
    BytesFromZZ(bytes,b,n);
    mpz_import(a, n, -1, 1, 1, 0, bytes);
    if( sign(b) == -1)
        mpz_neg(a,a);
    delete[] bytes;
}

void conv(fmpz_t a, const ZZ &b)
{
    //me being lazy
    mpz_t tmp;
    mpz_init(tmp);
    conv(tmp, b);
    fmpz_set_mpz(a,tmp);
    mpz_clear(tmp);
}


void conv(ZZ &a,const mpz_t b)
{
    clear(a);
    size_t n = (mpz_sizeinbase( b , 2) + 7) / 8;
    unsigned char* bytes = new unsigned char[n];
    mpz_export(bytes, &n, -1, 1, 1, 0, b);
    ZZFromBytes(a, bytes,n);
    if(mpz_sgn(b) == -1)
        a = -a;
    delete[] bytes;
}

void conv(ZZ &a, const fmpz_t b)
{
    //me being lazy
    mpz_t tmp;
    mpz_init(tmp);
    fmpz_get_mpz(tmp,b);
    conv(a,tmp);
    mpz_clear(tmp);
}



void conv(mpz_t* A, const Mat<ZZ> &B)
{
    int64_t rowsB = B.NumRows();
    int64_t colsB = B.NumCols();
    int64_t i, j;

    for(i = 0; i < rowsB; ++i)
    {
        for(j = 0; j < colsB; j++)
        {
            conv(A[i*rowsB+j], B[i][j]);
        }
    }
}


void conv(Mat<ZZ> &A, mpz_t* B, int64_t rowsB, int64_t colsB)
{
    int64_t i, j;
    A.SetDims(rowsB,colsB);
    for(i = 0; i < rowsB; ++i)
    {
        for(j = 0; j < colsB; j++)
        {
            conv(A[i][j], B[i*rowsB + j]);
        }
    }

}

