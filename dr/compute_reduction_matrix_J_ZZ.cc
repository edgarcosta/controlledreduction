// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
   
#include "dr.h"


/*
 * returns the iterator for the reduction matrix from V_{u+v} to V_u as a polynomial in u0,...,un
 * if not computed adds it to the map
 */
map< Vec<int64_t>, Vec<Mat<ZZ> >, vi64less>::const_iterator de_Rham_local::compute_reduction_matrix_J_ZZ(const Vec<int64_t> v)
{
    map< Vec<int64_t> ,  Vec<Mat<ZZ> >, vi64less>::const_iterator it_ZZ;
    it_ZZ = reduction_matrix_J_ZZ_dict.find(v);
    if(it_ZZ != reduction_matrix_J_ZZ_dict.end())
    {
        return it_ZZ;
    }
    else
    {
        map< Vec<int64_t> ,  Vec<Mat<ZZ_p> >, vi64less>::const_iterator it;
        it = reduction_matrix_J_dict.find(v);
        if( it == reduction_matrix_J_dict.end() )
        {
            it = compute_reduction_matrix_J(v);
            reduction_matrix_J_ZZ_dict[v] = conv< Vec< Mat<ZZ> > >(it->second);
            if( save_memory )
                //intel compiler dislikes
                //reduction_matrix_J_dict.erase( it );
                reduction_matrix_J_dict.erase( v );
        }
        else
            reduction_matrix_J_ZZ_dict[v] = conv< Vec< Mat<ZZ> > >(it->second);

        return reduction_matrix_J_ZZ_dict.find(v);
    }
}
