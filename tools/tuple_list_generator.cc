// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "tools.h"

#include <stdint.h>

#include <NTL/vector.h>
using namespace NTL;

void tuple_list_generator(Vec< Vec<int64_t> > &result, const int64_t d, const int64_t n)
{
    result.SetLength(binomial(d+n-1,n-1));
    int64_t * tuple = new int64_t[n];
    for(int64_t i=0; i < n - 1; ++i)
    {
        tuple[i] = 0;
    }
    tuple[n - 1] = d;
    int64_t j = 0;
    conv(result[j], tuple, n);

    while(tuple[0] < d)
    {
        for(int64_t i = n - 1; i > 0; --i)
        {
            if(tuple[i] > 0)
            {
                tuple[i-1]++;
                tuple[i]--;
                int64_t sum = 0;
                for(; i < n-1; i++)
                {
                    sum += tuple[i];
                    tuple[i] = 0;
                }
                tuple[n-1] += sum;
                break;
            }
                
        }
        j++;
        conv(result[j], tuple, n);
    }
    assert( j+1 == binomial(d+n-1,n-1) );
    delete[] tuple;
}
