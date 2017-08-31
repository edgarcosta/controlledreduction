// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "vec_int64.h"

#include <cstdint>
#include <assert.h>


#include <NTL/vector.h>

/*
 * returns v-ei
 */
Vec<int64_t> diff(const Vec<int64_t> v, const int64_t i)
{
    Vec<int64_t> r(v);
    assert(r[i]>0);
    --r[i];
    return r;
}


