// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.

#include "dr.h"

//basically a constructor
void de_Rham_local::init(int64_t p, int64_t precision, map< Vec<int64_t>, zz_p, vi64less> fbar, bool verbose, bool save_memory)
{
    int64_t i, j;
    ZZX pol, H;
    pair< Vec<int64_t>, Mat<ZZ_p> >* solve_pair;

    map< Vec<int64_t>, zz_p, vi64less>::iterator it;


    this->p = p;
    this->precision = precision;
    this->verbose = verbose;
    this->save_memory = save_memory;
    this->fbar = fbar;
    for(it = fbar.begin(); it != fbar.end(); it++)
    {
        f[it->first] = rep(it->second);
    }

    it = fbar.begin();
    n = (int64_t) (it->first).length() - 1;
    d = 0;
    for(i = 0; i <= n; i++)
        d += (it->first)[i];

    if(verbose)
    {
        cout<<"n = "<<n;
        cout<<" d = "<<d;
        cout<<" p = "<<p;
        cout<<" precision = "<<precision;
        cout<<endl;
        cout <<"fbar = \n";
        cout <<= fbar;
        cout <<endl;
    }


    assert(d>n);

    assert(d%p != 0);

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

    // Hilbert_J = (1 + x + x^2 + ... + x^(d-2) ) ^(n+1)
    Hilbert_J.SetLength((d - 2) * (n + 1) + 1);

    for(i = 0; i < d -1; i++)
        SetCoeff(pol, i);

    H = 1;

    for(i = 0; i< n + 1; i++)
        mul(H, H, pol);

    for(i = 0; i <= (d - 2)*(n + 1) ; i++)
        Hilbert_J[i] = (int64_t) to_ulong(coeff(H, i));

    if(verbose)
    {
        cout<<"Hilbert_J = "<<Hilbert_J<<endl;
    }

    if(verbose)
        cout << "Asserting that f is smooth.\n";

    solve_pair = get_solve_J((d - 2)*(n + 1) + 1);
    assert((solve_pair->first).length() == 0);

    for( i = 0 ; i < n; i++)
    {
        solve_pair = get_solve_J((i+1)*d-(n+1));
        assert( Hilbert_J[(i+1)*d-(n+1)] == (int64_t) (solve_pair->first).length() );
        for( j = 0; j < (int64_t)(solve_pair->first).length(); j++)
        {
            append(coKernels_J_basis, tuple_list[ (i + 1) * d - (n + 1) ][ (solve_pair->first)[j] ]);
        }

    }
    for( i = 0; i < (int64_t)coKernels_J_basis.length(); i++ )
    {
        coKernels_J_basis_dict[ coKernels_J_basis[i] ] = i;
    }

}

