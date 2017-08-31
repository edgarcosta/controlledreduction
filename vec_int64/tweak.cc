// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include <NTL/vector.h>

#include <cstdint>

#include "vec_int64.h"


NTL::Vec<int64_t> tweak(const NTL::Vec<int64_t> v, const int64_t r) {
    if ( r == 0 ) {
        NTL::Vec<int64_t> r;
        r = v;
        return r;
    } else {
        return tweak(tweak_step(v), r - 1);
    }
}

