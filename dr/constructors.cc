// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
// dr.cpp: routines for the class de_Rham_local
   
#include "dr.h"
//constructors

de_Rham_local::de_Rham_local(int64_t p, int64_t precision, map< Vec<int64_t>, zz_p, vi64less> fbar, bool verbose, bool save_memory)
{
    init(p, precision, fbar, verbose, save_memory);
}

de_Rham_local::de_Rham_local(int64_t p, int64_t precision, int64_t n, int64_t d, Vec<zz_p> fbar_vector, bool verbose, bool save_memory)
{
    Vec< Vec<int64_t> > tuple_d;
    map< Vec<int64_t>, zz_p, vi64less> fbar_map;
    int64_t i;
    tuple_list_generator(tuple_d, d, (int64_t) (n + 1) );
    assert(tuple_d.length() == fbar_vector.length() );
    for( i = 0; i < (int64_t) fbar_vector.length() ; i++)
        fbar_map[ tuple_d[i] ] = fbar_vector[i];

    init(p, precision, fbar_map, verbose, save_memory);

}

de_Rham_local::de_Rham_local(int64_t p, int64_t precision, int64_t n, int64_t d, bool verbose, bool save_memory)
{
    int64_t i ,j;
    Vec<int64_t> B;
    zz_p x;
    map< Vec<int64_t>, zz_p, vi64less>::iterator it;

    bool boolean = true;
    this->precision = 1;
    this->verbose = verbose;
    this->n = n;
    this->d = d;

    tuple_list.SetLength( (n + 1) * d + 1);
    tuple_dict.SetLength( (n + 1) * d + 1);
    for(i = 0 ; i < (n + 1) * d + 1 ; i++)
    {
        tuple_list_generator( tuple_list[i], i, (int64_t) (n+1));
        for( j = 0; j < (int64_t)tuple_list[i].length(); j++)
        {
            tuple_dict[i][ ((tuple_list[i])[j]) ] = j;
        }
    }

    while(boolean)
    {
        solve_J.clear();
        for( i = 0; i < (int64_t) tuple_list[d].length(); i++ )
        {
            x = random_zz_p();
            fbar[ tuple_list[d][i] ] = x;
            f[ tuple_list[d][i] ] = rep( x );
        
        }
        B = (get_solve_J((d-2)*(n+1)+1))->first;
        boolean = bool( B.length() != 0 );
    }

    solve_J.clear();
    this->precision = precision;
    
    init(p, precision, fbar, verbose, save_memory);

    
}

de_Rham_local::de_Rham_local(const char* filename)
{
    int64_t i, j;
    ifstream file;
    Vec<zz_p> fbar_vector;
    pair< Vec<int64_t>, Mat<ZZ_p> > *solve_pair;
    file.open(filename);

    if(file.is_open())
    {
        file >> p;
        assert( p == (int64_t)  zz_p::modulus() );
        file >> precision;
        assert(  power_ZZ(p, precision) == ZZ_p::modulus() );
        file >> verbose;
        file >> save_memory;
        file >> n;
        file >> d;
        file >> Hilbert_J;
        file >> fbar_vector;
        
        tuple_list.SetLength( (n + 1) * d + 1);
        tuple_dict.SetLength( (n + 1) * d + 1);
        for(i = 0 ; i < (n + 1) * d + 1 ; i++)
        {
            tuple_list_generator( tuple_list[i], i, (int64_t) (n+1));
            for( j = 0; j < (int64_t)tuple_list[i].length(); j++)
            {
                tuple_dict[i][ ((tuple_list[i])[j]) ] = j;
            }
        }
        assert( tuple_list[d].length() == fbar_vector.length() );
        for(i = 0 ; i < (int64_t) fbar_vector.length() ; i++)
        {
            if( ! IsZero(fbar_vector[i]) )
            {
                f[ tuple_list[d][i] ] = rep( fbar_vector[i] );
                fbar[ tuple_list[d][i] ] = fbar_vector[i];
            }
        }

        for(i = 0; i < n; i++)
        {
            solve_pair = &( solve_J[(i+1)*d - (n+1)]);
            file >> solve_pair->first;
            file >> solve_pair->second;
            assert( Hilbert_J[(i+1)*d-(n+1)] == (int64_t) (solve_pair->first).length() );
            for( j = 0; j < (int64_t)(solve_pair->first).length(); j++)
            {
                append(coKernels_J_basis, tuple_list[ (i + 1) * d - (n + 1) ][ (solve_pair->first)[j] ]);
            }
        }
        solve_pair = &( solve_J[(d - 2) * (n + 1) + 1]);
        file >> solve_pair->first;
        file >> solve_pair->second;

        for( i = 0; i < (int64_t)coKernels_J_basis.length(); i++ )
        {
            coKernels_J_basis_dict[ coKernels_J_basis[i] ] = i;
        }
        
    }
    else
    {
        cout << "Could not open \'"<<filename<<"\'"<<endl;
        abort();
    }
}






















