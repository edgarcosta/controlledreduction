// Copyright 2013-2017 Edgar Costa
// See LICENSE file for license details.
   
#include "dr.h"

bool isSmooth(const map< Vec<int64_t>, zz_p, vi64less> &f )
{
    int64_t i, j;
    de_Rham_local D;
    D = de_Rham_local();
    D.precision = 1;
    D.verbose = false;
    map< Vec<int64_t>, zz_p, vi64less>::const_iterator it;
    it = f.begin();
    D.n = (it->first).length() - 1;
    D.d = 0;
    for(i = 0; i <=D.n ; i++)
        D.d += (it->first)[i];

    for(it = f.begin(); it != f.end() ; it++ )
        D.f[it->first] = rep( it->second);


    D.tuple_list.SetLength( (D.n + 1) * D.d + 1);
    D.tuple_dict.SetLength( (D.n + 1) * D.d + 1);
    for(i = 0 ; i < (D.n + 1) * D.d + 1 ; i++)
    {
        tuple_list_generator( D.tuple_list[i], i, (int64_t) (D.n+1));
        for( j = 0; j < (int64_t)D.tuple_list[i].length(); j++)
        {
            D.tuple_dict[i][ ((D.tuple_list[i])[j]) ] = j;
        }
    }

    Mat<ZZ_p> M;
    D.matrix_J(M, (D.d-2)*(D.n+1) + 1 );
    Mat<zz_p> M_zz_p;
    M_zz_p.SetDims(M.NumRows(), M.NumCols());
    for(i = 0; i < (int64_t) M.NumRows(); i++)
    {
        for(j = 0; j < (int64_t) M.NumCols(); j++)
        {
            M_zz_p[i][j] = conv<zz_p>(rep(M[i][j]));
        }
    }
    nmod_mat_t M_flint;
    nmod_mat_init(M_flint,M.NumRows(), M.NumCols(), zz_p::modulus());
    long * P;
    bool result;
    P = (long *)flint_malloc(sizeof(long) * M.NumRows());
    conv(M_flint, M_zz_p);
    result = bool (M.NumRows() == nmod_mat_lu(P, M_flint, false));
    flint_free(P);
    nmod_mat_clear(M_flint);

    return result;
}
