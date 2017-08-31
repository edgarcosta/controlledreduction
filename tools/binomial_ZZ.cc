// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "tools.h"

#include <assert.h>
#include <stdint.h>

#include <NTL/ZZ.h>
using namespace NTL;

ZZ binomial_ZZ(const int64_t n, const int64_t k)
{
    assert(n>=k);
    ZZ num(1) ; // n!/k! = (k+1) ... n
    ZZ den(1); // (n-k)!
    for(int64_t i = k + 1; i <= n; ++i)
    {
        num *= i;
    }
    for(int64_t i = 1; i <= n - k; ++i)
    {
        den *=i;
    }
    return num / den;
}
