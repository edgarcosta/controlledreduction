// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "dr.h"

/*
 * computes inclusion_matrix_J(u) is the matrix that maps
 * V_u - > {(n-1)! H Omega/f^n : deg H = d*n - n - 1}
 * and automatically adds it to the dict
 */
void de_Rham_local::compute_inclusion_matrix_J(Vec<int64_t> u)
{
    assert( inclusion_matrix_J_dict.find(u) == inclusion_matrix_J_dict.end() );
    Vec< Vec<int64_t> >* list_G;
    map< Vec<int64_t>, int64_t, vi64less >* dict_H;
    Mat<ZZ_p>* M;
    Vec<int64_t> w;
    int64_t i, coordinate_of_monomial;
    bool boolean;

    list_G = &tuple_list[ d * n - n ];
    dict_H = &tuple_dict[ d * n - n - 1 ];

    M = &(inclusion_matrix_J_dict[u]);
    M->SetDims( dict_H->size(), list_G->length());
    for(coordinate_of_monomial = 0; coordinate_of_monomial < (int64_t) list_G->length(); coordinate_of_monomial++)
    {
        w = (*list_G)[coordinate_of_monomial] + u;
        boolean = true;
        for( i = 0; i <= n ; i++)
        {
            if(w[i] > 0)
                w[i]--;
            else
            {
                boolean = false;
                break;
            }
        } 
        if(boolean)
            NTL::set( (*M)[  (*dict_H)[w] ][ coordinate_of_monomial ] );
    }
}

