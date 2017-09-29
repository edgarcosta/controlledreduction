// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
   
#include "dr.h"
/*
* computes the matrix map
*  (H0,...,Hn) -> H0 * df0 + ... + Hn * dfn
*  Where Hi have degree = l - (d-1).
*/

void de_Rham_local::matrix_J(Mat<ZZ_p> &result, int64_t l)
{
    //free storage and make 0 x 0
    result.kill();

    if( d - 1 > l)
    {
        result.SetDims( tuple_list[l].length(),0);
        return;
    }

    int64_t i, j;
    int64_t  len_list_lminusd1;
    

    Vec< Vec<int64_t> > *tuple_list_l;
    Vec< Vec<int64_t> > *tuple_list_lminusd1;

    map< Vec<int64_t>, int64_t, vi64less> *tuple_dict_l;
    map< Vec<int64_t>, ZZ_p, vi64less>::const_iterator itf;
   
    tuple_list_l = &tuple_list[l];
    tuple_list_lminusd1 = &tuple_list[l - (d - 1)];
    
    tuple_dict_l = &tuple_dict[l];
    
    len_list_lminusd1 = tuple_list_lminusd1->length();

    
    result.SetDims( tuple_list_l->length() , (n+1) * len_list_lminusd1);
    
    for(i = 0 ; i <= n ; i++)
    {
        for(j = 0; j < len_list_lminusd1; j++)
        {
            for(itf = f.begin(); itf != f.end(); itf++)
            {
                if( (itf->first)[i] > 0 )
                {
                    result[ (*tuple_dict_l)[ diff(itf->first, i) + (*tuple_list_lminusd1)[j] ]][ i*len_list_lminusd1 + j ] = (itf->first)[i] * itf->second;
                }
            }
         }
    }
}
