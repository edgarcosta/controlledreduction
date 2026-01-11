// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
// solve_system.cpp: routines to solve systems
  


#include "solve_system.h"
#include "conv.h"
#include "NTL/mat_ZZ.h"

using namespace NTL;

void solve_system_zz_p_flint(Vec<int64_t> &B, nmod_mat_t U, const nmod_mat_t T)
{
    /*
     * input T: m * n matrix over Fp, representing a Fp-linear map from Fp^n to Fp^m
     *
     * output: B, U 
     *
     * This function computes a subset B of {0, 1, ..., m-1}, such that the basis elements e_i of Fp^m, for i in B, descend to a basis of Fp^m / im(T).
     *
     * Let J be the matrix of the corresponding inclusion Fp^N -> Fp^m
     * Consider the map T + J from Fp^n \oplus Fp^B to Fp^m,
     * By definition of B, this map is surjective
     *
     * This function also computes a right inverse of T+J.
     * i.e. a map U from Fp^m to Fp^n \oplus Fp^B, such that (T + J) U = identity on Fp^m.
     *
     * In other words,  this function shows how to write every element of Fp^m
     * as a linear combination of the columns of T, plus possibly some basis
     * elements of Fp^m not hit by T.
     *
     * Here B is a vector containing the indices corresponding to the basis elements, in increasing order.
     * 
     */
    int64_t nrows, ncols, rank, position;
    mp_limb_t modulus = zz_p::modulus();
    int64_t i, j;
    nmod_mat_t T_transpose;
    nmod_mat_t T_copy;
    nmod_mat_t Y; // Y = pivot cols of T, plus J
    nmod_mat_t Z; // inverse of Y
    Vec<int64_t> pivot_cols;

    nrows = nmod_mat_nrows(T);
    ncols = nmod_mat_ncols(T);
    if(ncols == 0)
    {
        B.SetLength(nrows);
        nmod_mat_init(U, nrows, nrows, modulus);
        for(i = 0; i < nrows; i++)
        {
            B[i] = i;
            nmod_mat_entry(U, i, i) = 1;
        }
        return;
    }

    nmod_mat_init_set(T_copy, T);
    nmod_mat_init(T_transpose, ncols, nrows, modulus);
    
    nmod_mat_transpose(T_transpose, T);

    //Compute the pivot cols of T
    long * P;
    P = (long *)flint_malloc(sizeof(long) * nrows);
    rank = nmod_mat_lu(P, T_copy, false);
    flint_free(P);
    pivot_cols.SetLength(rank);
    j = 0;
    for( i = 0; i < rank; i++)
    {
        while(nmod_mat_entry(T_copy, i, j) == 0 and j < ncols)
        {
            j++;
        }
        pivot_cols[i] = j;
        j++;
    }

    nmod_mat_clear(T_copy);

    //Compute part of Y
    nmod_mat_init(Z, nrows, nrows, modulus);
    nmod_mat_init(Y, nrows, nrows,modulus);
    for( i = 0; i < rank; i++)
    {
        position = pivot_cols[i];
        for( j = 0; j < nrows; j++)
        {
            nmod_mat_entry(Y, j, i) = nmod_mat_entry(T_transpose, position, j);
        }
    }
    //Compute B
    
    //T_transpose becomes its row echelon form
    P = (long *)flint_malloc(sizeof(long) * ncols);
    rank = nmod_mat_lu(P, T_transpose ,false);
    flint_free(P);
    assert(rank == (int64_t)pivot_cols.length());
    // B = set of nonpivot cols of Tranpose(T)
    // = set of nonpivot rows of T if we don't row swapping 
    B.SetLength(nrows - rank);
    
    j = 0;
    position = 0;
    for( i = 0; i < ncols; i++)
    {
        while(j < nrows and nmod_mat_entry(T_transpose, i, j) == 0)
        {
            B[position] = j;
            position++;
            j++;
        }
        j++;
    }
    for( ; j < nrows ; j++)
    {
        B[position] = j;
        position++;
    }

     //Finishing Y
    for( i = 0; i < nrows - rank; i++)
    {
        nmod_mat_entry(Y, B[i], rank + i) = 1;

    }
    
    nmod_mat_clear(T_transpose);
    nmod_mat_inv(Z, Y);
    nmod_mat_clear(Y);

    //Recall
    // m = nrows
    // n = ncols
    // |B| = nrows - rank
    //U from Fp^m to Fp^n \oplus Fp^B, such that (T + J) U = identity on Fp^m
    
    nmod_mat_init(U, ncols+nrows-rank, nrows, modulus);

    for(i = 0; i < rank; i++)
    {
        position = pivot_cols[i];
        for(j = 0; j < nrows; j++)
        {
            nmod_mat_entry(U, position, j) = nmod_mat_entry(Z, i, j);
        }
    }

    for(i = rank; i < nrows; i++)
    {
        position = ncols - rank + i;
        for(j = 0; j < nrows; j++ )
        {
            nmod_mat_entry(U, position, j) = nmod_mat_entry(Z, i, j);
        }
    }
    nmod_mat_clear(Z);
}



    






 




void solve_system_zz_p_NTL(Vec<int64_t> &B, Mat<zz_p> &U, const Mat<zz_p> &T)
{
    /*
     * input T: m * n matrix over Fp, representing a Fp-linear map from Fp^n to Fp^m
     *
     * output: B, U 
     *
     * This function computes a subset B of {0, 1, ..., m-1}, such that the basis elements e_i of Fp^m, for i in B, descend to a basis of Fp^m / im(T).
     *
     * Let J be the matrix of the corresponding inclusion Fp^N -> Fp^m
     * Consider the map T + J from Fp^n \oplus Fp^B to Fp^m,
     * By definition of B, this map is surjective
     *
     * This function also computes a right inverse of T+J.
     * i.e. a map U from Fp^m to Fp^n \oplus Fp^B, such that (T + J) U = identity on Fp^m.
     *
     * In other words,  this function shows how to write every element of Fp^m
     * as a linear combination of the columns of T, plus possibly some basis
     * elements of Fp^m not hit by T.
     *
     * Here B is a vector containing the indices corresponding to the basis elements, in increasing order.
     * 
     */

    int64_t nrows, ncols, rank, position;
    int64_t i, j;
    Mat<zz_p> T_transpose;
    Mat<zz_p> T_copy;
    Mat<zz_p> Y; // Y = pivot cols of T, plus J
    Mat<zz_p> Z; // inverse of Y
    Vec<int64_t> pivot_cols;

    nrows = T.NumRows();
    ncols = T.NumCols();
    if(ncols == 0)
    {
        B.SetLength(nrows);
        clear(U);
        U.SetDims(nrows,nrows);
        for(i = 0; i < nrows; i++)
        {
            B[i] = i;
            NTL::set(U[i][i]);
        }
        return;
    }
    T_copy = T;
    
    transpose(T_transpose, T);
    
    //Compute the pivot cols of T
    rank = gauss(T_copy);
    pivot_cols.SetLength(rank);
    j = 0;
    for( i = 0; i < rank; i++)
    {
        while(T_copy[i][j] == 0 and j < ncols)
        {
            j++;
        }
        pivot_cols[i] = j;
        j++;
    }

    //Compute part of Y
    Y.SetDims(nrows, nrows);
    for( i = 0; i < rank; i++)
    {
        position = pivot_cols[i];
        for( j = 0; j < nrows; j++)
        {
            Y[j][i] = T_transpose[ position ][j];
        }
    }



    //Compute B
    
    //T_transpose becomes its row echelon form
    rank = gauss( T_transpose );
    // B = set of nonpivot cols of Tranpose(T)
    // = set of nonpivot rows of T if we don't row swapping 
    B.SetLength(nrows - rank);
    
    j = 0;
    position = 0;
    for( i = 0; i < ncols; i++)
    {
        while(j < nrows and T_transpose[i][j] == 0)
        {
            B[position] = j;
            position++;
            j++;
        }
        j++;
    }
    for( ; j < nrows ; j++)
    {
        B[position] = j;
        position++;
    }
    //Finishing Y
    for( i = 0; i < nrows - rank; i++)
    {
        NTL::set(Y[ B[i] ][rank + i ]);

    }
    inv(Z, Y);
    //Recall
    // m = nrows
    // n = ncols
    // |B| = nrows - rank
    //U from Fp^m to Fp^n \oplus Fp^B, such that (T + J) U = identity on Fp^m
    clear(U);
    U.SetDims(ncols + nrows  - rank, nrows);
    
    for(i = 0; i < rank; i++)
    {
        position = pivot_cols[i];
        for(j = 0; j < nrows; j++)
        {
            U[position][j] = Z[i][j];
        }
    }

    for(i = rank; i < nrows; i++)
    {
        position = ncols - rank + i;
        for(j = 0; j < nrows; j++ )
        {
            U[position][j] = Z[i][j];
        }
    }
}

void solve_system_zz_p(Vec<int64_t> &B,  Mat<zz_p> &U, const Mat<zz_p> &T, bool flint)
{
    if(!flint)
        solve_system_zz_p_NTL(B,U,T);
    else
    {
        nmod_mat_t U_flint;
        nmod_mat_t T_flint;
        nmod_mat_init(T_flint, T.NumRows(), T.NumCols(), zz_p::modulus());
        conv(T_flint, T);
        solve_system_zz_p_flint(B,U_flint,T_flint);
        conv(U,U_flint);
        nmod_mat_clear(T_flint);
        nmod_mat_clear(U_flint);
    }
}



void solve_system_padic_flint(Vec<int64_t> &B, Mat<ZZ_p> &U, const Mat<ZZ_p> &T, int64_t precision)
{
    int64_t nrows, ncols, Blength;
    int64_t i, j;
    mp_limb_t modulus = conv<ulong>(ZZ_p::modulus());

    nmod_mat_t T_Fp;
    nmod_mat_t U_Fp;

    nrows = T.NumRows();
    ncols = T.NumCols();
    if(ncols == 0)
    {
        B.SetLength(nrows);
        U.SetDims(nrows,nrows);
        for(i = 0; i < nrows; i++)
        {
            B[i] = i;
            NTL::set(U[i][i]);
        }
        return;
    }
    nmod_mat_init(T_Fp, nrows, ncols, zz_p::modulus());
    conv(T_Fp, T);
    solve_system_zz_p_flint(B, U_Fp, T_Fp);

    Blength = B.length();
   
    
    if((ulong) modulus == ZZ_p::modulus())
    {
        nmod_mat_t X; // X = T + J
        //Compute X
        nmod_mat_init(X, nrows, ncols + Blength, modulus);

        for(i = 0; i < nrows; i++)
        {
            for(j = 0; j < ncols; j++)
            {
                nmod_mat_entry(X, i, j) = conv<ulong>(T[i][j]);
            }
        }
        for(i = 0; i < Blength; i++)
        {
            nmod_mat_entry(X, B[i], i + ncols) = 1;
        }

        nmod_mat_t U_flint;
        nmod_mat_init(U_flint, ncols + Blength, nrows, modulus);

    
        // lift U_Fp
        for(i = 0; i < ncols + Blength; i++)
        {
            for(j = 0; j < nrows; j++)
            {
                nmod_mat_entry(U_flint, i, j) = nmod_mat_entry(U_Fp, i, j);
            }
        }


        nmod_mat_t XU, UXU;

        nmod_mat_init(XU, nrows, nrows, modulus);
        nmod_mat_init(UXU, ncols + Blength, nrows, modulus);
        
        j = 1;
        while( j < precision)
        {
            // U = 2U - UXU
            nmod_mat_mul(XU, X, U_flint);
            nmod_mat_mul(UXU, U_flint, XU);
            nmod_mat_scalar_mul(U_flint, U_flint, 2);
            nmod_mat_sub(U_flint, U_flint, UXU);
            j *=2;
        }
        conv(U,U_flint);
        nmod_mat_clear(X);
        nmod_mat_clear(U_flint);
        nmod_mat_clear(XU);
        nmod_mat_clear(UXU);

    }
    else
    {
        fmpz_t fmpz_modulus;
        fmpz_init(fmpz_modulus);
        fmpz_set_ui(fmpz_modulus, zz_p::modulus() );

        fmpz_mat_t X; // X = T + J
        //Compute X
        fmpz_mat_init(X, nrows, ncols + Blength);

        for(i = 0; i < nrows; i++)
        {
            for(j = 0; j < ncols; j++)
            {
                conv(fmpz_mat_entry(X, i, j), rep(T[i][j]));
            }
        }
        for(i = 0; i < Blength; i++)
        {
            fmpz_set_ui(fmpz_mat_entry(X, B[i], i + ncols), 1);
        }

        fmpz_mat_t U_flint;
        fmpz_mat_init(U_flint, ncols + Blength, nrows);

        // lift U_Fp
        for(i = 0; i < ncols + Blength; i++)
        {
            for(j = 0; j < nrows; j++)
            {
                fmpz_set_ui(fmpz_mat_entry(U_flint, i, j), nmod_mat_entry(U_Fp, i, j));
            }
        }


        fmpz_mat_t XU, UXU;

        fmpz_mat_init(XU, nrows, nrows);
        fmpz_mat_init(UXU, ncols + Blength, nrows);
        
        int64_t l, k;
        j = 1;
        while( j < precision)
        {
            // U = 2U - UXU
            fmpz_mat_mul(XU, X, U_flint); 
            fmpz_mat_mul(UXU, U_flint, XU);
            fmpz_mat_scalar_mul_ui(U_flint, U_flint, 2);
            fmpz_mat_sub(U_flint, U_flint, UXU);
            
            j *=2;
            if(j < precision)
            {
                fmpz_mul(fmpz_modulus,fmpz_modulus,fmpz_modulus);
                for(l = 0; l <  ncols + Blength; l++)
                    for(k = 0; k < nrows; k++)
                        fmpz_mod(fmpz_mat_entry(U_flint, l, k), fmpz_mat_entry(U_flint, l, k), fmpz_modulus);
            }
        }
        conv(U,U_flint);
        fmpz_mat_clear(X);
        fmpz_mat_clear(U_flint);
        fmpz_mat_clear(XU);
        fmpz_mat_clear(UXU);
        fmpz_clear(fmpz_modulus);

    }
    nmod_mat_clear(T_Fp);
    nmod_mat_clear(U_Fp);
}


void solve_system_padic_NTL(Vec<int64_t> &B, Mat<ZZ_p> &U, const Mat<ZZ_p> &T, int64_t precision)
{
    int64_t nrows, ncols, Blength;
    int64_t i, j;

    Mat<zz_p> T_Fp;
    Mat<zz_p> U_Fp;
    Mat<ZZ_p> X; // X = T + J
    Mat<ZZ_p> Utranspose;
    Mat<ZZ_p> tmp0, tmp1;


    nrows = T.NumRows();
    ncols = T.NumCols();
    if(ncols == 0)
    {
        B.SetLength(nrows);
        U.SetDims(nrows,nrows);
        for(i = 0; i < nrows; i++)
        {
            B[i] = i;
            NTL::set(U[i][i]);
        }
        return;
    }
    T_Fp.SetDims(nrows, ncols);
 

    for(i = 0; i < nrows; i++)
    {
        for(j = 0; j < ncols; j++)
        {
            T_Fp[i][j] = conv<zz_p>(rep(T[i][j]));
        }
    }

    solve_system_zz_p(B, U_Fp, T_Fp);
    Blength = B.length();

    //Compute X
    X.SetDims(nrows, ncols + Blength );
    for(i = 0; i < nrows; i++)
    {
        for(j = 0; j < ncols; j++)
        {
            X[i][j] = T[i][j];
        }
    }
    for(i = 0; i < Blength; i++)
    {
        NTL::set(X[B[i]][i+ncols]);
    }
    
    U.SetDims(ncols + Blength, nrows);
    // lift U_Fp
    for(i = 0; i < ncols + Blength; i++)
    {
        for(j = 0; j < nrows; j++)
        {
            U[i][j] = conv<ZZ_p>(rep(U_Fp[i][j]));
        }
    }
    j = 1;
    while( j < precision )
    {
        // lift from mod p^j to mod p^(2j)
        // U = 2U - UXU
        // U = 2U - U*(X*U^T)^T
        //printf("U = %lu x %lu, X = %lu x %lu \n", U.NumRows(), U.NumCols(), X.NumRows(), X.NumCols());
        //U = 2*U - U*(X*U);
        mul(tmp0, X, U);
        mul(tmp1, U, tmp0);
        mul(U, 2, U);
        sub(U, U, tmp1);
        /*
        transpose(Utranspose,U);
        transpose_mul_transpose(tmp0,X,Utranspose);
        mul_transpose(tmp1,U,tmp0);
        mul(tmp0,2,U);
        sub(U,tmp0,tmp1);*/
        j *= 2;
    }
}

void solve_system_padic(Vec<int64_t> &B, Mat<ZZ_p> &U, const Mat<ZZ_p> &T, int64_t precision, bool flint)
{
    if(flint)
        solve_system_padic_flint(B,U,T,precision);
    else
        solve_system_padic_NTL(B,U,T,precision);
}


