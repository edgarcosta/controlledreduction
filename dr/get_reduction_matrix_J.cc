// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
   
#include "dr.h"

/*
 * returns the reduction matrix from V_{u+v} to V_u
 */
Mat<ZZ_p> de_Rham_local::get_reduction_matrix_J(const Vec<int64_t> u, const Vec<int64_t> v)
{
    int64_t i , sum;
    Mat<ZZ_p> result;
    sum = 0;
    for(i = 0; i <= n; i++)
    {
        sum += v[i];
        assert(u[i]==0 || v[i] > 0 );
    }
    assert( sum == d);

    Mat<ZZ_p> temp;
    map< Vec<int64_t> ,  Vec<Mat<ZZ_p> >, vi64less>::const_iterator it;
    it = compute_reduction_matrix_J(v);
    result = it->second[0];
    for( i = 0; i <= n ; i++)
    {
        mul(temp, it->second[i+1], u[i]);
        add(result,result,temp);
    }
    return result;
}


