// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
   
#include "dr.h"


/*
* tests the path independence while performing reductions
*/
bool de_Rham_local::test_paths_J(int64_t trials, int64_t paths)
{
    if(verbose)
        cout << "de_Rham_local::test_paths_J("<<trials<<", "<<paths<<")\n";
    assert(trials > 0);
    assert(paths > 0);

    int64_t attempt, pathi, i, k, sum, sum_u;
    Vec<int64_t> u_orig, u, v;
    Vec< Mat<ZZ_p> > M_saved;
    Mat<ZZ_p> M;
    Mat<ZZ_p>* Mfinalred;
    Mat<ZZ_p>* Minclusion;
    Mat<ZZ_p> Mred;
    u_orig.SetLength(n+1);
    v.SetLength(n+1);
    
    for(attempt = 0; attempt < trials; attempt++)
    {
        if(verbose)
            cout << "attempt = "<<attempt<<endl;
        M_saved.kill();
        M_saved.SetLength(paths);
        sum = 0;
        for( i = 0 ; i <=n ; i++)
        {
            u_orig[i] = RandomBnd(20);
            sum += u_orig[i];
        }
        
        while( sum % d != n )
        {
            k = RandomBnd(n+1);
            u_orig[k]++;
            sum++;
        }


        for( pathi = 0 ; pathi < paths ; pathi++)
        {
            if(verbose)
                cout << "path = "<<pathi<<endl;
            M.kill();
            M.SetDims( tuple_list[d*n-n].length(), tuple_list[d*n-n].length());
            for( i = 0; i < (int64_t) tuple_list[d*n-n].length(); i++)
                NTL::set( M[i][i] );
            sum_u = 0;
            u = u_orig;
            for( i = 0; i < n+1; i++)
                sum_u += u[i];
            while( sum_u > n)
            {
                sum = 0;
                for( i = 0; i <= n; i++)
                {
                    if( u[i] > 0 )
                    {
                        v[i] = 1;
                        sum++;
                    }
                    else
                        v[i] = 0;
                }
                while( sum < d)
                {
                    k = RandomBnd(n + 1);
                    if( u[k] > 0 and v[k]<u[k])
                    {
                        v[k]++;
                        sum++;
                    }
                }
                u = u - v;
                sum_u -= d;

                Mred = get_reduction_matrix_J(u,v);
                mul(M, Mred, M);
            }
            Minclusion = get_inclusion_matrix_J(u);
            mul(M, *Minclusion, M);

            Mfinalred = get_final_reduction_matrix_J(n);
            mul(M, *Mfinalred, M);
            M_saved[pathi] = M;
            if(pathi > 0)
                if (M_saved[pathi-1] != M_saved[pathi])
                    return false;
        }
    }
    return true;
}
