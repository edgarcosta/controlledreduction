// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "vec_int64.h"

#include <cstdint>
#include <assert.h>


#include <NTL/vector.h>


Vec<int64_t> tweak( const Vec<int64_t> v, const int64_t r)
{
    if(r == 0)
    {
        Vec<int64_t> r;
        r = v;
        return r;
    }
    else
    {
        return tweak(tweak_step(v), r - 1);
    }
}

