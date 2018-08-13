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

void zeta_function(NTL::ZZX &zeta, const std::map< NTL::Vec<int64_t>, NTL::ZZ, vi64less> &f,const int64_t &p, const bool &verbose = false, const int threads = 1);

void zeta_function(NTL::ZZX &zeta, const std::vector< std::vector<int64_t> > &monomials, const std::vector<int64_t> &coef, const int64_t &p,  const bool &verbose = false, const int threads = 1);
void zeta_function(NTL::ZZX &zeta, const std::vector< std::vector<int64_t> > &monomials, const std::vector<NTL::ZZ> &coef, const int64_t &p,  const bool &verbose = false, const int threads = 1);


/*
 * same as above, but the input is given through string in the following format:
 *      p
 *      f.keys()
 *      f.values()
 */
void zeta_function(NTL::ZZX &zeta, const char* input, const bool &verbose = false, const int threads = 1);

#endif // WRAPPER_H_



