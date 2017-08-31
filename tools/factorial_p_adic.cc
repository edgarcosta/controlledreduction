// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "tools.h"

#include <assert.h>
#include <stdint.h>

#include <NTL/ZZ_p.h>

using namespace NTL;

int64_t factorial_p_adic(ZZ_p &result, const int64_t n, const int64_t p, const int64_t start)
{
    int64_t val = 0;
    result = 1;
    for (int64_t i = start; i <=n ; i++) {
        if (i % p != 0)
            result *= i;
        else {
            int64_t tmp = i;
            while ( tmp % p == 0) {
                tmp = tmp / p;
                ++val;
            }
            result *= tmp;
        }
    }
    if( start == 1)
        assert(val ==  valuation_of_factorial(n, p));

    return val;
}
