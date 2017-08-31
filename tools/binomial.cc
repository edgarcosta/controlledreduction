// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "tools.h"

#include <assert.h>
#include <stdint.h>

#include <NTL/ZZ.h>
using namespace NTL;

int64_t binomial(const int64_t n, int64_t k)
{
    assert(n >= k);
    // n!/k! = (k+1) ... n
    ZZ num(1); 
    // (n-k)!
    ZZ den(1);
    for(int64_t i = k + 1; i <= n; ++i)
    {
        num *= i;
    }
    for(int64_t i = 1; i <= n - k; ++i)
    {
        den *=i;
    }
    num = num/den;
    assert( NumBytes(num) <= 8 );
    return int64_t{conv<long>(num)};
}

