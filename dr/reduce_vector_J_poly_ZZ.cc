// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
   
#include "dr.h"

/*
 *  reduces a vector from V_{u+kv) to V_u
 *  just like reduce_vector_J_poly, but working temporarily over ZZ
 */
void de_Rham_local::reduce_vector_J_poly_ZZ(Vec<ZZ> &result, const Vec<int64_t> u, const Vec<int64_t> v, const int64_t iterations, const Vec<ZZ> G)
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

    Mat<ZZ> M0, M1;
    Vec<ZZ> Gin, Gout, Gout0, Gout1;
    int64_t x;

    M0 = it->second[0];
    M1.SetDims( G.length(), G.length() );
    for( i = 0; i <= n ; i++)
    {
        M0 += u[i] * it->second[i+1];
        M1 += v[i] * it->second[i+1];
    }

    Gin = G;
    for(x = iterations - 1; x != (int64_t)-1; x--)//counting on overflow
    {
        mul(Gout0, M0, Gin);
        mul(Gout1, M1, Gin);
        mul(Gout1,x,Gout1);
        add(Gout,Gout0,Gout1);
        for(i = 0; i < (int64_t) G.length(); i++)
            rem(Gin[i], Gout[i], ZZ_p::modulus() );
    }
    result = Gin;
}

