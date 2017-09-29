// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
   
#include "dr.h"

// reduces a vector from from V_{u+v} to V_u
// temporarily working over ZZ
void de_Rham_local::reduce_vector_J_ZZ(Vec<ZZ> &result, const Vec<int64_t> u, const Vec<int64_t> v, const Vec<ZZ> G)
{
    int64_t i, sum;
    sum = 0;
    for(i = 0; i <= n; i++)
    {
        sum += v[i];
        assert(u[i]==0 || v[i] > 0 );
        //u_zz_p[i] = conv<ZZ_p>(u[i]);
    }
    assert( sum == d);

        
    map< Vec<int64_t> ,  Vec<Mat<ZZ> >, vi64less>::const_iterator it;
    it = compute_reduction_matrix_J_ZZ(v);
    
    result = it->second[0] * G;
    for( i = 0; i <= n ; i++)
    {
        result += u[i] * (it->second[i+1] * G);
    }
}

