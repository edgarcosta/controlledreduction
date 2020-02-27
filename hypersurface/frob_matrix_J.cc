// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "hypersurface.h"
#ifdef _OPENMP
# include <omp.h>
#endif


/*
 * returns the matrix coordinates of the p-adic approximation of frob(e_i)
 * using N_k terms, where (k+1) is the pole order of e_i
 */
Mat<ZZ_p> hypersurface::frob_matrix_J(Vec<int64_t> N)
{
    assert( n == (int64_t) N.length() );
    Mat<ZZ_p> F;
    int64_t i;
    #ifdef _OPENMP
    if(omp_get_max_threads() > 1) {
      dR->compute_everything_J();
      compute_fpow(max(N) - 1);
    }
    ZZ_pContext context;
    context.save();
    #endif
    F.SetDims( dR->coKernels_J_basis.length(), dR->coKernels_J_basis.length() );

    #pragma omp parallel for schedule(dynamic)
    for( i = 0; i < (int64_t) dR->coKernels_J_basis.length(); i++)
    {
        #ifdef _OPENMP
        context.restore();
        #endif
        int64_t j, m, sum;
        sum = 0;
        for( j = 0; j<=n; j++)
            sum += dR->coKernels_J_basis[i][j];

        m = (sum + n + 1)/d;

        if(verbose)
            cout<<"Computing F("<<dR->coKernels_J_basis[i]<<") m = "<<m<<" N = "<<N[m-1]<<endl;
        if(N[m-1]>0)
            F[i] = frob_J(i, N[m-1]);
    }
    return transpose(F);
}


