// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include <assert.h>
#include <gmp.h>
#include <flint/fmpz.h>
#include <flint/fmpz_mat.h>
#include <NTL/matrix.h>
#include <NTL/mat_ZZ.h>
#include <NTL/ZZ.h>

#include <cstdint>

using namespace NTL;

//Some routines for "Mat"<mpz_t> and "Vec"<mpz_t>
// result = A * B
void mul(mpz_t* result, mpz_t* A, int64_t rowsA, int64_t colsA, mpz_t* B, int64_t rowsB, int64_t colsB)
{
    assert(colsA == rowsB);
    int64_t i,j,k;
    mpz_t tmp, acc;
    mpz_inits(tmp, acc, NULL);
    for(i = 0; i < rowsA; ++i)
    {
        for(j = 0; j < colsB; j++)
        {
            mpz_set_ui(acc,0);
            for(k = 0; k < rowsB; k++)
            {
                mpz_mul(tmp, A[i*rowsA + k], B[k*rowsB + j]);
                mpz_add(acc,acc,tmp);
            }
            mpz_set(result[i*rowsA + j],acc);
        }
    }
    mpz_clears(tmp, acc, NULL);
}


//result = a * A;
void mul(mpz_t* result, mpz_t a, mpz_t* A, int64_t rowsA, int64_t colsA)
{
    int64_t i, j;
    for(i = 0; i < rowsA; ++i)
    {
        for(j = 0; j < colsA; j++)
        {
            mpz_mul(result[i*rowsA+j], a, A[i*rowsA+j]);
        }
    }
}

// result = A * v
void mul(mpz_t* result, mpz_t* A, int64_t rowsA, int64_t colsA, mpz_t* v, int64_t v_size)
{
    assert(v_size == colsA);
    int64_t i,j;
    mpz_t tmp, acc;
    mpz_inits(tmp,acc,NULL);
    for(i = 0; i < rowsA; ++i)
    {
        mpz_set_ui(acc,0);
        for(j = 0; j < colsA; j++)
        {
            mpz_mul(tmp,A[i*rowsA +j],v[j]);
            mpz_add(acc,acc,tmp);
        }
        mpz_set(result[i],acc);
    }
}
// result = a * v;
void mul(mpz_t* result, mpz_t a, mpz_t* v, int64_t v_length)
{
    int64_t i;
    for(i = 0; i < v_length; ++i)
    {
        mpz_mul(result[i],a,v[i]);
    }
}

// result = a + b;
// suuports a = a + b
void mul(mpz_t* result, mpz_t* a, int64_t a_length, mpz_t * b, int64_t b_length)
{
    assert(a_length == b_length);
    int64_t i;
    for(i = 0; i < a_length; ++i)
        mpz_add(result[i],a[i],b[i]);
}
//result = <u,v>
void mul(mpz_t &result, mpz_t* u, int64_t u_length, mpz_t* v, int64_t v_length)
{
    assert(u_length == v_length);
    mpz_set_ui(result,0);
    mpz_t tmp;
    mpz_init(tmp);
    int64_t i;
    for(i = 0; i < u_length; ++i)
    {
        mpz_mul(tmp,u[i],v[i]);
        mpz_add(result, result, tmp);
    }
}

//Some routines for fmpz_{mat,vec} and nmod_{mat,vec}
//result = A * v
void mul(fmpz * result, const fmpz_mat_t A, const fmpz * v)
{
    int64_t ar, ac, i, j;
    ar = fmpz_mat_nrows(A);
    ac = fmpz_mat_ncols(A);
    for(i = 0; i < ar; ++i)
    {
        fmpz_mul(result + i,fmpz_mat_entry(A, i, 0), v);
        for(j = 1; j < ac; j++)
            fmpz_addmul(result + i, fmpz_mat_entry(A, i, j), v + j);
    }
}
//result = <v, w>;
void mul(fmpz_t result, const fmpz * v, const fmpz * w, const int64_t len)
{
    int64_t i;
    fmpz_mul(result, v, w);
    for(i = 1; i < len; ++i)
        fmpz_addmul(result, v + i, w + i);
}

