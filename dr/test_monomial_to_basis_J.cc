// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
   
#include "dr.h"

/*
* checks that a basis element x^w \Omega / f^k ~ x^w f^l \Omega / f^(l+k)
* for l <= N
*/

bool de_Rham_local::test_monomial_to_basis_J(int64_t N)
{
    if(verbose)
        cout <<"de_Rham_local::test_monomial_to_basis_J("<<N<<")"<<endl;

    Vec< map< Vec<int64_t>, ZZ_p, vi64less > > fpow;
    map< Vec<int64_t>, ZZ_p, vi64less>::iterator it1, it2;
    Vec<int64_t> u;
    Vec<ZZ_p> v;
    int64_t i, j, k, m, fact, sum;
    
    v.SetLength( coKernels_J_basis.length() );

    fpow.SetLength(N+1);

    u.SetLength(n+1);
    for( i = 0; i <= n; i++)
        u[i] = 0;

    NTL::set( fpow[0][u] );

    for( i = 0; i < N; i++)
    {
        for( it1 = fpow[i].begin() ; it1 != fpow[i].end() ; it1++ )
        {
            for( it2 = f.begin() ; it2 != f.end() ; it2++ )
            {
                u = it1->first + it2->first;
                fpow[i+1][u] += it1->second * it2->second;
            }
        }
    }
    for( i = 0; i <= N ; i++)
    {
        if(verbose)
            cout<<"i = "<<i<<endl;
        for( j = 0; j < (int64_t) coKernels_J_basis.length() ; j++ )
        {
            sum = 0;
            for(k = 0 ; k <= n; k++)
                sum += coKernels_J_basis[j][k];
            sum += (n + 1);
            m = sum/d + i;
            fact = 1;
            for( k = 1; k < m; k++)
            {
                fact *= k;
            }

            for( k = 0; k < (int64_t) coKernels_J_basis.length() ; k++ )
                v[k] = 0;

            for( it1 = fpow[i].begin(); it1 != fpow[i].end(); it1++ )
            {
                u = it1->first + coKernels_J_basis[j];
                for( k = 0 ; k <=n; k++ )
                    u[k]++;
                v += it1->second * monomial_to_basis_J(u);
            }
            for( k = 0 ; k < (int64_t) coKernels_J_basis.length() ; k++ )
            {
                if ( k == j and v[k] != to_ZZ_p( fact ) )
                    return false;
                if ( k !=j and not IsZero(v[k]) )
                    return false;
            }
        }
    }
    return true;
}
