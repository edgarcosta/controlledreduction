// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
// 
// matrix.h: matrix miscellaneous routines

#ifndef CONTROLLED_REDUCTION_MATRIX_H
#define CONTROLLED_REDUCTION_MATRIX_H

#include <assert.h>
#include <stdint.h>

#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <flint/nmod_mat.h>
#include <flint/nmod_vec.h>
#include <gmp.h>
#include <NTL/matrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/mat_ZZ_p.h>
#include <NTL/vector.h>
#include <NTL/ZZ_p.h>
#include <NTL/ZZ.h>

using namespace NTL;



//computes the charpoly using Newton identities of M
Vec<ZZ> charpoly(const Mat<ZZ> M);
//computes the charpoly using Newton identities of M, know bounds of the powersums and the precision of charpoly if we just used the newton identites
Vec<ZZ> charpoly_frob(const Mat<ZZ> M, Vec<int64_t> prec, const int64_t prime, const int64_t dimension);


//Some routines for "Mat"<mpz_t> and "Vec"<mpz_t>
// result = A * B
void mul(mpz_t* result, mpz_t* A, int64_t rowsA, int64_t colsA, mpz_t* B, int64_t rowsB, int64_t colsB);

// result = A * B where B is given in column major, ie, as the transpose
void mul_transpose(mpz_t* result, mpz_t* A, int64_t rowsA, int64_t colsA, mpz_t* B, int64_t rowsB, int64_t colsB);

//result = a * A;
void mul(mpz_t* result, mpz_t a, mpz_t* A, int64_t rowsA, int64_t colsA);

// result = A * v
void mul(mpz_t* result, mpz_t* A, int64_t rowsA, int64_t colsA, mpz_t* v, int64_t v_size);

// result = a * v
void mul(mpz_t* result, mpz_t a, mpz_t* v, int64_t v_length);

// result = a + b;
// supports a = a + b
void mul(mpz_t* result, mpz_t* a, int64_t a_length, mpz_t * b, int64_t b_length);

//result = <u,v>
void mul(mpz_t &result, mpz_t* u, int64_t u_length, mpz_t* v, int64_t v_length);


//Some routines for fmpz_{mat,vec} and nmod_{mat,vec}
//result = A * v
void mul(fmpz * result,const fmpz_mat_t A,const fmpz * v);
//result = <v, w>;
void mul(fmpz_t result,const fmpz * v,const fmpz * w, const int64_t len);

//result = A * v;
void inline mul(mp_ptr result, const nmod_mat_t A, mp_srcptr v, const int nlimbs)
{
    int64_t ar, ac, i ;
    ar = nmod_mat_nrows(A);
    ac = nmod_mat_ncols(A);
    for( i = 0; i < ar; ++i)
        result[i] = _nmod_vec_dot((A->rows)[i], v, ac, A->mod, nlimbs);       
}


// computes result = A * transpose(B)
template<class T>
void mul_transpose(Mat<T> &result,const Mat<T> &A, const Mat<T> &B)
{
    int64_t rowsA = A.NumRows();
    int64_t colsA = A.NumCols();
    int64_t colsB = B.NumRows();
    int64_t rowsB = B.NumCols();
    assert(colsA == rowsB);
    result.SetDims(rowsA,colsB);
    int64_t i, j, k;
    T acc, tmp;
    for(i = 1; i <= rowsA; ++i)
    {
        for(j = 1; j <= colsB; j++)
        {
            clear(acc);
            for(k = 1; k <= rowsB; k++)
            {
                mul(tmp,A(i,k),B(j,k));
                add(acc,acc,tmp);
            }
            result(i,j) = acc;
        }
    }
}


// computes result = A * transpose(B)
template <> void mul_transpose<ZZ_p>(Mat<ZZ_p> &result,const Mat<ZZ_p> &A, const Mat<ZZ_p> &B);
extern template void mul_transpose<ZZ_p>(Mat<ZZ_p> &result,const Mat<ZZ_p> &A, const Mat<ZZ_p> &B);



void inline nmod_mat_sub(nmod_mat_t C, const nmod_mat_t A, const nmod_mat_t B)
{
    slong i;

    if (C->c == 0)
        return;

    for (i = 0; i < C->r; ++i)
    {
        _nmod_vec_sub(C->rows[i], A->rows[i], B->rows[i], C->c, C->mod);
    }
}

// returns a random SL_n matrix
Mat<zz_p> random_SL_matrix(int64_t n);

// returns the trace of M
ZZ trace(const Mat<ZZ> M);

// computes result = transpose( A * transpose(B) )
template<class T>
void transpose_mul_transpose(Mat<T> &result,const Mat<T> &A, const Mat<T> &B)
{
    int64_t rowsA = A.NumRows();
    int64_t colsA = A.NumCols();
    int64_t colsB = B.NumRows();
    int64_t rowsB = B.NumCols();
    assert(colsA == rowsB);
    result.SetDims(colsB,rowsA);
    int64_t i, j, k;
    T acc, tmp;
    for(j = 1; j <= colsB; j++)
    {
        for(i = 1; i <= rowsA; ++i )
        {
            clear(acc);
            for(k = 1; k <= rowsB; k++)
            {
                mul(tmp,A(i,k),B(j,k));
                add(acc,acc,tmp);
            }
            result(j,i) = acc;
        }
    }
}

template <> void transpose_mul_transpose<ZZ_p>(Mat<ZZ_p> &result,const Mat<ZZ_p> &A, const Mat<ZZ_p> &B);
extern template void transpose_mul_transpose<ZZ_p>(Mat<ZZ_p> &result,const Mat<ZZ_p> &A, const Mat<ZZ_p> &B);




#endif
