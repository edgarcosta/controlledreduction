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

void zeta_function(NTL::ZZX &zeta, const std::map< NTL::Vec<int64_t>, int64_t, vi64less> &f,const int64_t &p, bool verbose = false);

void zeta_function(NTL::ZZX &zeta, std::vector< std::vector<int64_t> > &monomials, std::vector<int64_t> &coef, int64_t p,  bool verbose = false);

#endif // WRAPPER_H_



