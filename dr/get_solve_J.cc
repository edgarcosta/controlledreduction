// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
#include "dr.h"
#include "timing.h"

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

    timestamp_type wtime1, wtime2;
    double wall_time, user_time;

    Mat<ZZ_p> MJ;
    matrix_J(MJ, level);
    if(verbose) {
        user_time = get_cpu_time();
        get_timestamp(&wtime1);
        cout << "Computing and solving matrix of relations at degree = "<<level<<" ( "<<MJ.NumRows()<<"x"<<MJ.NumCols()<<" )."<<endl;
    }

    solve_system_padic(solve_J[level].first, solve_J[level].second, MJ, precision);
    if (verbose) {
        get_timestamp(&wtime2);
        wall_time = timestamp_diff_in_seconds(wtime1,wtime2);
        user_time = get_cpu_time() - user_time;
        printf("Time: CPU %.2f s, Wall: %.2f s\n", user_time, wall_time );
    }
    return &(solve_J[ level ]);
}
