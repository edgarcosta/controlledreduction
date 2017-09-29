// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
   
#include "dr.h"

/*
 * returns a pointer to the inclusion_matrix_J(u)
 * which is the matrix that maps
 * V_u - > {(n-1)! H Omega/f^n : deg H = d*n - n - 1}
 * and automatically adds it to the dict
 */


Mat<ZZ_p>* de_Rham_local::get_inclusion_matrix_J(Vec<int64_t> u)
{
    int64_t i, sum;
    map< Vec<int64_t>, Mat<ZZ_p>, vi64less>::iterator it;

    sum = 0;
    for( i = 0; i <=n ; i++)
        sum += u[i];
    assert(sum == n);

    it = inclusion_matrix_J_dict.find(u);
    if( it == inclusion_matrix_J_dict.end() )
    {
        if( verbose )
            cout << "Computing the inclusion matrix J for u = "<<u<<endl;
        compute_inclusion_matrix_J(u);
        it = inclusion_matrix_J_dict.find(u);
        assert( it !=  inclusion_matrix_J_dict.end() ); 
    }
    return &(it->second);
}

