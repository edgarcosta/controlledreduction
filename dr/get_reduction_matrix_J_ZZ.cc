// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
   
#include "dr.h"

/*
 * returns the reduction matrix (over ZZ) from V_{u+v} to V_u
 * that will be correct in ZZ_p
 */
Mat<ZZ> de_Rham_local::get_reduction_matrix_J_ZZ(const Vec<int64_t> u, const Vec<int64_t> v)
{
    int64_t i , sum;
    Mat<ZZ> result;
    sum = 0;
    for(i = 0; i <= n; i++)
    {
        sum += v[i];
        assert(u[i]==0 || v[i] > 0 );
    }
    assert( sum == d);

    map< Vec<int64_t>, Vec<Mat<ZZ> >, vi64less>::const_iterator it;
    it =  compute_reduction_matrix_J_ZZ(v);
    result = it->second[0];
    for( i = 0; i <= n ; i++)
    {
        result += it->second[i+1]*u[i];
    }
    return result;
}


