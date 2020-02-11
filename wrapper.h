// Copyright 2017 Edgar Costa
// See LICENSE file for license details.
// wrapper.h: header file for the wrapper that computes the Frobenius matrix and the zeta function

#ifndef WRAPPER_H_
#define WRAPPER_H_

#include "vec_int64.h"
#include <cstdint>
#include <map>
#include <vector>
#include <NTL/vector.h>
#include <NTL/ZZX.h>
#include <NTL/mat_ZZ.h>
#include <NTL/vec_ZZ.h>

void zeta_function(
    NTL::ZZX &zeta,  // output: the zeta function
    NTL::Mat<NTL::ZZ> &Frob_ZZ, // output: the Frobenius matrix
    const std::map< NTL::Vec<int64_t>, NTL::ZZ, vi64less> &f, // f as vector
    const int64_t &p, // the prime p
    const bool &verbose = false, //enable/disable verbose mode
    const int &threads = 1, //number of threads
    const int &min_abs_precision = 0, // in case we want Frob correct mod p^min_abs_precision,
    const bool &find_better_model = true // if one should try to find a non-degenerate model, this usually speeds up the overall computation
    );

void zeta_function(
    NTL::ZZX &zeta,  // output: the zeta function
    NTL::Mat<NTL::ZZ> &Frob_ZZ, // output: the Frobenius matrix
    // f as vector monomials + coeffs
    const std::vector< std::vector<int64_t> > &monomials,
    const std::vector<int64_t> &coef,
    const int64_t &p, // the prime p
    const bool &verbose = false, //enable/disable verbose mode
    const int &threads = 1, //number of threads
    const int &min_abs_precision = 0, // in case we want Frob correct mod p^min_abs_precision,
    const bool &find_better_model = true // if one should try to find a non-degenerate model, this usually speeds up the overall computation
    );

void zeta_function(
    NTL::ZZX &zeta,  // output: the zeta function
    NTL::Mat<NTL::ZZ> &Frob_ZZ, // output: the Frobenius matrix
    // f as vector monomials + coeffs
    const std::vector< std::vector<int64_t> > &monomials,
    const std::vector<NTL::ZZ> &coef,
    const int64_t &p, // the prime p
    const bool &verbose = false, //enable/disable verbose mode
    const int &threads = 1, //number of threads
    const int &min_abs_precision = 0, // in case we want Frob correct mod p^min_abs_precision,
    const bool &find_better_model = true // if one should try to find a non-degenerate model, this usually speeds up the overall computation
    );


/*
 * same as above, but the input is given through string in the following format:
 *      p
 *      f.keys()
 *      f.values()
 */
void zeta_function(
    NTL::ZZX &zeta,  // output: the zeta function
    NTL::Mat<NTL::ZZ> &Frob_ZZ, // output: the Frobenius matrix
    const char* input,
    const bool &verbose = false, //enable/disable verbose mode
    const int &threads = 1, //number of threads
    const int &min_abs_precision = 0, // in case we want Frob correct mod p^min_abs_precision,
    const bool &find_better_model = true // if one should try to find a non-degenerate model, this usually speeds up the overall computation
    );

#endif // WRAPPER_H_



