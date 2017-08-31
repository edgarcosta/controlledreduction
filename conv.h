// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
// tools.h: header file for miscellaneous routines

#include <assert.h>
#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/nmod_mat.h>
#include <NTL/matrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vector.h>
#include <NTL/ZZ.h>


//conversion between NTL::ZZ, mpz and fmpz
//b -> a

void conv(mpz_t a, const NTL::ZZ &b);
void conv(fmpz_t a, const NTL::ZZ &b);
void conv(NTL::ZZ &a, const mpz_t b);
void conv(NTL::ZZ &a, const fmpz_t b);

//conversion between NTL::Vec/NTL::Mat<T> and it's equivalents over fmpz_t, mpz_t and "nmod"
void conv(mpz_t* A, const NTL::Mat<NTL::ZZ> &B);

template<class T>
void conv(fmpz_mat_t A, const NTL::Mat<T> &B)
{
    assert( B.NumRows() == fmpz_mat_nrows(A));
    assert( B.NumCols() == fmpz_mat_ncols(A));
    for (int64_t i = 0; i < B.NumRows(); ++i)
        for (int64_t j = 0; j < B.NumCols(); j++)
            conv(fmpz_mat_entry(A, i, j), NTL::conv<NTL::ZZ>(B[i][j]));
}

template<class T>
void conv(fmpz * a, const NTL::Vec<T> &b) {
    for(int64_t i = 0; i < b.length(); ++i)
        conv(a + i, NTL::conv<NTL::ZZ>(b[i]));
}

template<class T>
void conv(nmod_mat_t A, const NTL::Mat<T> &B)
{
    assert( B.NumRows() == nmod_mat_nrows(A));
    assert( B.NumCols() == nmod_mat_ncols(A));
    for (int64_t i = 0; i < B.NumRows(); ++i)
        for (int64_t j = 0; j < B.NumCols(); j++)
            nmod_mat_entry(A, i, j) = NTL::conv<ulong>(rep(B[i][j]) % A->mod.n );
}

void conv(NTL::Mat<NTL::ZZ> &A, mpz_t* B, int64_t rowsB, int64_t colsB);

template<class T>
void conv(NTL::Mat<T> &A, const fmpz_mat_t B)
{
    int64_t rowsA = fmpz_mat_nrows(B);
    int64_t colsA = fmpz_mat_ncols(B);
    A.SetDims(rowsA, colsA);
    NTL::ZZ tmp;
    for (int64_t i = 0; i < rowsA; ++i) {
        for (int64_t j = 0; j < colsA; j++) {
            conv(tmp, fmpz_mat_entry(B, i, j));
            conv(A[i][j], tmp);
        }
    }
}

template<class T>
void conv(NTL::Vec<T> &a, const fmpz * b, int64_t len) {
    a.SetLength(len);
    NTL::ZZ tmp;
    for(int64_t i = 0; i < len; ++i) {
        conv(tmp, b + i);
        conv(a[i], tmp);
    }
}


template<class T>
void conv(NTL::Mat<T> &A, const nmod_mat_t B) {
    int64_t rowsA = nmod_mat_nrows(B);
    int64_t colsA = nmod_mat_ncols(B);
    A.SetDims(rowsA, colsA);
    for(int64_t i = 0; i < rowsA; ++i)
        for(int64_t j = 0; j < colsA; j++)
            conv(A[i][j], nmod_mat_entry(B, i, j));
}
