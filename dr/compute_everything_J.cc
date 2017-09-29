// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
   
#include "dr.h"

/*
 * computes all the reduction and inclusion matrices
 */
void de_Rham_local::compute_everything_J()
{
    int64_t i;
    Vec<int64_t> u;
    for(i = 0; i < (int64_t) tuple_list[d].length(); i++)
    {
        u = tuple_list[d][i];
        compute_reduction_matrix_J( u );
        compute_reduction_matrix_J_ZZ( u );
    }

    for(i = 0; i  < (int64_t) tuple_list[n].length(); i++)
    {
        u = tuple_list[n][i];
        get_inclusion_matrix_J(u);
    }
    get_final_reduction_matrix_J(n); 
}
