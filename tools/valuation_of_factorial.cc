// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include <cstdint>

int64_t valuation_of_factorial(const int64_t n, const int64_t p)
{
    int64_t sum = 0;
    int64_t tmp = n;
    while ( tmp != 0 ) {
        sum += tmp%p;
        tmp = tmp/p;
    }
    return (n - sum)/(p - 1);
}
