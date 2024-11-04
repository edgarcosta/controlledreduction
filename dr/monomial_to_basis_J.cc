// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
   
#include "dr.h"

/*
 * Computes the coordinates of (m-1)! x^u \Omega / x0 ... xn f^m
 * in the cohomology basis
 */

Vec<ZZ_p> de_Rham_local::monomial_to_basis_J(const Vec<int64_t> u)
{
    assert( (int64_t) u.length() == n+1 );
    int64_t i, m, sum;

    sum = 0;
    for( i = 0; i <= n ; i++)
    {
        sum += u[i];
        assert( u[i] > 0 );
    }
    assert( sum%d == 0);

    m = sum / d;
    
    if( m > n)
    {
        int64_t e;
        Vec<int64_t> s, v, w;
        Mat<ZZ_p> M_red;
        Mat<ZZ_p>* M;
        Vec<ZZ_p> G, G_new;
        G.SetLength( tuple_list[d*n-n].length() );
        s.SetLength(n+1);
        v.SetLength(n+1);
        w = tweak( u, d * n - n );
        sub(s, u, w);
        // u = w + s
        NTL::set( G[ tuple_dict[d*n - n][s] ] );

        for( e = m; e > n; e--)
        {
            sum = 0;
            for( i = 0; i <= n; i++)
            {
                if(w[i] > 0)
                {
                    v[i] = 1;
                    sum++;
                }
                else
                {
                    v[i] = 0;
                }
            }
            while( sum < d )
            {
                for( i = 0; i <= n; i++)
                {
                    if( w[i] > 0 && v[i] < w[i] )
                    {
                        v[i]++;
                        sum++;
                        break;
                    }
                }
            }
            //w = w - v;
            sub(w,w,v);
            M_red = get_reduction_matrix_J(w, v);
            mul(G, M_red, G);
            //mul(G_new, M_red, G);
            //swap(G_new, G);
            
        }
        sum = 0;
        for(i = 0; i <= n; i++)
            sum += w[i];

        assert(sum == n);
        // Now we have:
        // (m-1)! x^u \Omega / x0...xn ~ (n-1)! x^w \Omega / x0 ... xn f^n
        M = get_inclusion_matrix_J(w);

        //mul(G_new, *M, G);
        mul(G, *M, G);


        M = get_final_reduction_matrix_J(n);
        //mul(G, *M, G_new);
        mul(G, *M, G);

        return G;
    }
    else
    {
        Vec<ZZ_p> G, G_new;
        Vec<int64_t> w;
        Mat<ZZ_p>* M;

        G.SetLength( tuple_list[ sum - (n+1)].length() );
        w = u;
        for( i = 0; i <= n ; i++)
            w[i]--;

        NTL::set( G[ tuple_dict[ sum - (n+1)][w] ] );
        M = get_final_reduction_matrix_J(m);
        mul(G_new, *M, G);

        return G_new;
    }
}

