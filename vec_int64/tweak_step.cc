// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "vec_int64.h"

#include <cstdint>
#include <assert.h>


#include <NTL/vector.h>

Vec<int64_t> tweak_step( const Vec<int64_t> v)
{
    Vec<int64_t> r(v);
    assert((int64_t) r.length() == v.length());
    for(int64_t i = 0 ; i < r.length() ; ++i)
    {
        if(r[i] > 0)
        {
            --r[i];
            break;
        }
    }
    return r;
}
