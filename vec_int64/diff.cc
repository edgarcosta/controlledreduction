// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include <NTL/vector.h>

#include <cassert>
#include <cstdint>


/*
 * returns v-ei
 */
NTL::Vec<int64_t> diff(const NTL::Vec<int64_t> v, const int64_t i) {
    NTL::Vec<int64_t> r(v);
    assert(r[i] > 0);
    --r[i];
    return r;
}


