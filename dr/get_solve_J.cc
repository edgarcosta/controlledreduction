// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
   
#include "dr.h"

/*
 * computes the pair solve_Jl and adds it to the dictionary (if not compute already)
 * returns the pointer for the pair
 */
pair< Vec<int64_t>, Mat<ZZ_p> >* de_Rham_local::get_solve_J(int64_t level)
{
    map< int64_t, pair< Vec<int64_t>, Mat<ZZ_p> > >::iterator it;
    it = solve_J.find(level);
    if( it != solve_J.end() )
        return &(it->second);


    Mat<ZZ_p> MJ;
    matrix_J(MJ, level);

    if(verbose)
        cout << "Computing and solving matrix of relations at degree = "<<level<<" ( "<<MJ.NumRows()<<"x"<<MJ.NumCols()<<" )."<<endl;
    solve_system_padic(solve_J[level].first, solve_J[level].second, MJ, precision);
    return &(solve_J[ level ]);
}
