// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
   
#include "dr.h"

/*
 *  reduces a vector from V_{u+kv) to V_u
 */
void de_Rham_local::reduce_vector_J_poly(Vec<ZZ_p> &result, const Vec<int64_t> u, const Vec<int64_t> v, const int64_t iterations, const Vec<ZZ_p> G)
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
    map< Vec<int64_t> ,  Vec<Mat<ZZ_p> >, vi64less>::const_iterator it;
    it = compute_reduction_matrix_J(v);

    Mat<ZZ_p>  M0, M1;
    Vec<ZZ_p> Gin,  Gout0, Gout1;
    //Vec<ZZ_p> Gout;
    int64_t x;

    M0 = it->second[0];
    M1.SetDims( G.length(), G.length() );
    for( i = 0; i <= n ; i++)
    {
        M0 += u[i] * it->second[i+1];
        M1 += v[i] * it->second[i+1];
    }

    Gin = G;
    for(x = iterations - 1; x != (int64_t)-1 ; x--) //couting on overflow
    {
        //Gout = M0 * Gin;
        //Gout += x * (M1 * Gin);
        mul(Gout0, M0, Gin);
        mul(Gout1, M1, Gin);
        mul(Gout1,x,Gout1);
        add(Gin,Gout0,Gout1);
        //swap(Gin,Gout);
    }
    result = Gin;
}

