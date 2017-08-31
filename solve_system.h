// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
// solve_system.h: header file for the routines to solve systems
   

#ifndef SOLVE_SYSTEM_H
#define SOLVE_SYSTEM_H

#include <NTL/mat_ZZ_p.h>
#include <NTL/mat_lzz_p.h>
#include "tools.h"
#include <stdint.h>
#include <assert.h>

using namespace std;
using namespace NTL;
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
void solve_system_zz_p_NTL(Vec<int64_t> &B, Mat<zz_p> &U, const Mat<zz_p> &T);
void solve_system_zz_p_flint(Vec<int64_t> &B, nmod_mat_t U, const nmod_mat_t T);
void solve_system_zz_p(Vec<int64_t> &B,  Mat<zz_p> &U, const Mat<zz_p> &T, bool flint = true);


/* 
* Same as above but over the Z/p^precision
*/

void solve_system_padic(Vec<int64_t> &B, Mat<ZZ_p> &U, const Mat<ZZ_p> &T, int64_t precision, bool flint = true);
void solve_system_padic_flint(Vec<int64_t> &B, Mat<ZZ_p> &U, const Mat<ZZ_p> &T, int64_t precision);
void solve_system_padic_NTL(Vec<int64_t> &B, Mat<ZZ_p> &U, const Mat<ZZ_p> &T, int64_t precision);

#endif
