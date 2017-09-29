// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
   
#include "dr.h"

/*
 * returns the final reduction matrix for
 * (k -1)! G \Omega / f^k to the cohomology basis elements
 */
Mat<ZZ_p>* de_Rham_local::get_final_reduction_matrix_J(int64_t k)
{
    assert(k <= n);
    assert(k > 0);
    int64_t i;
    if( k > (int64_t)final_reduction_matrix_J_dict.length() )
    {
        i = final_reduction_matrix_J_dict.length() + 1;
        final_reduction_matrix_J_dict.SetLength(k);
        for(  ; i <= k ; i++)
        {
            if( verbose )
                cout <<"Computing the final reduction matrix J for pole = "<<i<<endl;
            compute_final_reduction_matrix_J(i);
        }
    }
    return &(final_reduction_matrix_J_dict[k-1]);
}
