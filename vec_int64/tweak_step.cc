// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.


#include <NTL/vector.h>

#include <cassert>
#include <cstdint>

#include "vec_int64.h"

NTL::Vec<int64_t> tweak_step(const NTL::Vec<int64_t> v) {
    NTL::Vec<int64_t> r(v);
    assert((int64_t) r.length() == v.length());
    for (int64_t i = 0 ; i < r.length(); ++i) {
        if (r[i] > 0) {
            --r[i];
            break;
        }
    }
    return r;
}
